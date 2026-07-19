from typing import List, Optional
import bisect
import subprocess
import tempfile
import os
import re
from collections import namedtuple

# Try importing the Python anarci package. If unavailable, fall back to anarci.exe.
try:
    import anarci as _anarci_pkg
    _HAS_ANARCI_PKG = True
except ImportError:
    _anarci_pkg = None
    _HAS_ANARCI_PKG = False

CDRS = {
    'chothia': {
        'L': {
            1: slice(24, 35),
            2: slice(50, 57),
            3: slice(89, 98),
        },
        'H': {
            1: slice(26, 33),
            2: slice(52, 57),
            3: slice(95, 103),
        }
    }
}

# insertion_code defaults to ' ' to be consistent with ANARCI.
class ResidueCode(namedtuple('ResidueCode', ['index', 'insertion_code'], defaults=[' '])):
    def __init__(self, index, insertion_code=' '):
        if not isinstance(index, int):
            raise TypeError('index must be an integer')

class Chain:
    """
    Provides indexing functionality for a subset of the query residues.

    Stores a subset of the residues of the original query, as well as its own
    start index (offset)
    """
    def __init__(self, residues, offset):
        self.residues = residues
        self.offset = offset
        return

    def _integer_index(self, key: Optional[ResidueCode]):
        if key is None:
            return None
        return bisect.bisect_left(self.residues, key)

    @staticmethod
    def regularize_key(key):
        if isinstance(key, slice):
            start = None if key.start is None else Chain.to_residue_code(key.start)
            stop = None if key.stop is None else Chain.to_residue_code(key.stop)
            return slice(start, stop, key.step)
        else:
            return slice(Chain.to_residue_code(key), Chain.to_residue_code(key) + ('',))

    @staticmethod
    def to_residue_code(key):
        try:
            return ResidueCode(*key)
        except TypeError:  # not iterable
            return ResidueCode(key)

    def __getitem__(self, key):
        key = self.regularize_key(key)
        # "and" evaluates to 2nd operand if the 1st operand is False-like.
        # In this case, the start and stop elements are None or ResidueCode.
        start = self._integer_index(key.start)
        stop = self._integer_index(key.stop)
        max_range = range(self.offset, self.offset + len(self.residues))
        return max_range[start:stop:key.step]

    def __repr__(self):
        return f'{self.__class__.__name__}({self.residues}, offset={self.offset})'


class Annotation:
    def __init__(self, seq, run_name='ab', scheme='chothia'):
        self.seq = seq
        self.run_name = run_name
        self.scheme = scheme
        
        if _HAS_ANARCI_PKG:
            self._init_from_python_pkg(seq, run_name, scheme)
        else:
            self._init_from_subprocess(seq, run_name, scheme)

    def _init_from_python_pkg(self, seq, run_name, scheme):
        """Use the Python anarci package (WSL/Linux)."""
        numbering, alignment_details, hit_tables = _anarci_pkg.anarci([(run_name, seq)], scheme=scheme)
        self._chains = {}
        for chain, details in zip(numbering[0], alignment_details[0]):
            assert details['chain_type'] not in self._chains
            residues, start, stop = chain
            res_keys = [r[0] for r in residues if r[1] != '-']
            assert is_sorted(res_keys)
            assert stop - start == len(res_keys) - 1, f"{start=}, {len(res_keys)=}, {stop=}"
            self._chains[details['chain_type']] = Chain(res_keys, start)

    def _init_from_subprocess(self, seq, run_name, scheme):
        """Use anarci.exe via subprocess (Windows fallback)."""
        from ..platform_config import get_config
        cfg = get_config()
        anarci_exe = cfg.anarci
        
        if anarci_exe is None or not os.path.exists(anarci_exe):
            raise RuntimeError(
                "ANARCI is not available. Install the Python 'anarci' package, "
                "or place anarci.exe in Tools/ANARCI/.")
        
        # Write sequence to temp FASTA file
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, f"{run_name}.fasta")
            with open(fasta_path, 'w') as f:
                f.write(f">{run_name}\n{seq}\n")
            
            result = subprocess.run(
                [anarci_exe, '--sequence', fasta_path, '--scheme', scheme],
                capture_output=True, text=True, timeout=120,
            )
            
            if result.returncode != 0:
                raise RuntimeError(f"ANARCI failed: {result.stderr}")
        
        self._chains = _parse_anarci_text(result.stdout, scheme)

    @classmethod
    def from_traj(cls, traj, **kwargs):
        seq = get_sequence(traj.top)
        return cls(seq, **kwargs)

    def chain(self, name) -> Chain:
        try:
            return self._chains[name]
        except KeyError:
            if name == 'light':
                return self.light_chain()
            elif name == 'heavy':
                return self.heavy_chain()
            raise
    
    def light_chain(self) -> Chain:
        if 'K' in self._chains:
            return self._chains['K']
        return self._chains['L']
    
    def heavy_chain(self) -> Chain:
        return self._chains['H']

    def cdr_indices(self, cdr_def=None) -> List[int]:
        """Indices in the input sequence that correspond to the 6 CDR loops."""
        if cdr_def is None:
            cdr_def = CDRS[self.scheme]
        out = []
        heavy = self.chain('heavy')
        for loop in cdr_def['H'].values():
            for aa in heavy[loop]:
                out.append(aa)
        light = self.chain('light')
        for loop in cdr_def['L'].values():
            for aa in light[loop]:
                out.append(aa)
        return out

    def __repr__(self):
        return f'{self.__class__.__name__}({self.seq}, {self.run_name}, {self.scheme})'


def is_oneletter(r):
    if r is None:
        return False
    return len(r) == 1


def get_sequence(top):
    seq = []
    for r in top.residues:
        seq.append(r.code if is_oneletter(r.code) else 'X')
    assert len("".join(seq)) == top.n_residues
    return "".join(seq)

def is_sorted(lst):
    for i, element in enumerate(lst[1:]):
        if lst[i] > element:
            return False
    return True


def _parse_anarci_text(output: str, scheme: str = 'chothia'):
    """Parse ANARCI text output into Chain objects.
    
    The output format is:
        H 1       E
        H 2       V
        H 52    A P    (insertion code A in Chothia numbering)
        ...
        //
    
    Returns dict of chain_type -> Chain objects compatible with Annotation._chains.
    """
    import re
    
    chains_data = {}  # chain_type -> list of (chothia_number, insertion_code, seq_position)
    current_chain = None
    seq_position = 0
    
    for line in output.split('\n'):
        line = line.rstrip()
        if not line or line.startswith('#'):
            continue
        if line.strip() == '//':
            current_chain = None
            continue
        
        # Parse: H 26      G   or   H 52    A P
        match = re.match(r'^([HKL])\s+(\d+)\s*([A-Z]?)\s+([A-Z])', line)
        if not match:
            continue
        
        chain_type = match.group(1)
        chothia_num = int(match.group(2))
        insertion = match.group(3).strip() if match.group(3) else ''
        aa = match.group(4)
        
        if aa == '-':
            continue  # skip gaps
        
        if chain_type != current_chain:
            current_chain = chain_type
            seq_position = 0
            if chain_type not in chains_data:
                chains_data[chain_type] = []
        
        chains_data[chain_type].append((chothia_num, insertion, seq_position))
        seq_position += 1
    
    # Convert to Chain objects compatible with existing API
    # The Chain class expects: residues = sorted list of ResidueCode, offset = min chothia number
    
    chains = {}
    for chain_type, data in chains_data.items():
        # Map chothia_number -> seq_position; also handle insertion codes
        res_keys = []
        for chothia_num, insertion, seq_pos in data:
            ic = insertion if insertion else ' '
            res_keys.append(ResidueCode(chothia_num, ic))
        
        res_keys.sort()
        offset = min(r.index for r in res_keys)
        chains[chain_type] = Chain(res_keys, offset)
    
    # Map chain types: ANARCI uses H/K/L, internal code also supports 'heavy'/'light'
    if 'H' not in chains and 'heavy' in chains:
        chains['H'] = chains['heavy']
    if 'K' in chains and 'L' not in chains:
        chains['L'] = chains['K']
    
    return chains
