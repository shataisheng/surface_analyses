"""
Platform auto-detection and tool path configuration for PEP-Patch.

Automatically detects Windows / WSL / native Linux and provides
paths to external binary tools (MSMS, APBS, pdb2pqr).

Usage:
    from surface_analyses.platform_config import get_config
    cfg = get_config()
    print(cfg.msms)       # path to msms binary
    print(cfg.apbs_dir)   # default APBS working directory
    cfg.setup_path()      # add tools to os.environ['PATH']
"""

import os
import sys
import shutil
import pathlib
from dataclasses import dataclass, field


def _is_wsl() -> bool:
    """Detect if running under Windows Subsystem for Linux."""
    if sys.platform != "linux":
        return False
    try:
        with open("/proc/version", "r") as f:
            return "microsoft" in f.read().lower() or "wsl" in f.read().lower()
    except (FileNotFoundError, PermissionError):
        pass
    # Fallback: check for WSL-specific env vars
    return bool(os.environ.get("WSL_DISTRO_NAME"))


def _is_windows() -> bool:
    return sys.platform == "win32"


@dataclass
class PlatformConfig:
    """Holds platform-specific paths to external tools."""
    
    platform: str  # "windows", "wsl", "linux"
    project_root: pathlib.Path
    tools_dir: pathlib.Path
    
    # Tool paths (relative to tools_dir or absolute system paths)
    msms_dir: pathlib.Path = field(init=False)
    apbs_dir: pathlib.Path = field(init=False)
    pdb2pqr_dir: pathlib.Path = field(init=False)
    anarci_dir: pathlib.Path = field(init=False)
    
    # Default APBS working directory (where intermediate files go)
    default_apbs_work_dir: str = field(init=False)
    
    def __post_init__(self):
        if self.platform == "windows":
            self._init_windows()
        elif self.platform == "wsl":
            self._init_wsl()
        else:
            self._init_linux()
    
    def _init_windows(self):
        self.msms_dir = self.tools_dir / "msms"
        self.apbs_dir = self.tools_dir / "APBS-3.4.1.Windows" / "bin"
        self.pdb2pqr_dir = self.tools_dir / "pdb2pqr-portable"
        self.anarci_dir = self.tools_dir / "ANARCI"
        self.default_apbs_work_dir = str(self.project_root / "Tools" / "apbs_work")
    
    def _init_wsl(self):
        # On WSL, tools are typically installed system-wide.
        # We still look in the project tools dir first, then fall back to system.
        self.msms_dir = self.tools_dir / "msms"
        self.apbs_dir = self.tools_dir / "APBS-3.4.1.Windows" / "bin"
        self.pdb2pqr_dir = self.tools_dir / "pdb2pqr-portable"
        self.anarci_dir = self.tools_dir / "ANARCI"
        self.default_apbs_work_dir = os.environ.get(
            "PEP_PATCH_APBS_DIR",
            str(self.project_root / "Tools" / "apbs_work")
        )
    
    def _init_linux(self):
        self.msms_dir = self.tools_dir / "msms"
        self.apbs_dir = self.tools_dir / "APBS-3.4.1.Windows" / "bin"
        self.pdb2pqr_dir = self.tools_dir / "pdb2pqr-portable"
        self.anarci_dir = self.tools_dir / "ANARCI"
        self.default_apbs_work_dir = os.environ.get(
            "PEP_PATCH_APBS_DIR",
            str(self.project_root / "Tools" / "apbs_work")
        )
    
    def _resolve(self, name: str, bundled: pathlib.Path) -> str:
        """Resolve external tool path: env override > bundled Tools/ > PATH > name.

        Bundled binary is preferred so a missing uv-installed git wrapper
        (e.g. msms-wrapper) does not break tool discovery.
        """
        env_val = os.environ.get(f"PEP_PATCH_{name.upper()}")
        if env_val and os.path.exists(env_val):
            return env_val
        if bundled.exists():
            return str(bundled)
        on_path = shutil.which(name)
        if on_path:
            return on_path
        return name

    @property
    def msms(self) -> str:
        """Return the msms executable path (bundled first)."""
        ext = "msms.exe" if self.platform == "windows" else "msms"
        return self._resolve("msms", self.msms_dir / ext)
    
    @property
    def apbs(self) -> str:
        """Return the apbs executable path (bundled first)."""
        ext = "apbs.exe" if self.platform == "windows" else "apbs"
        return self._resolve("apbs", self.apbs_dir / ext)
    
    @property
    def pdb2pqr(self) -> str:
        """Return the pdb2pqr executable path (bundled first)."""
        ext = "pdb2pqr.exe" if self.platform == "windows" else "pdb2pqr"
        return self._resolve("pdb2pqr", self.pdb2pqr_dir / ext)
    
    @property
    def anarci(self) -> str | None:
        """Return the anarci executable path, or None if not available."""
        env_val = os.environ.get("PEP_PATCH_ANARCI")
        if env_val and os.path.exists(env_val):
            return env_val
        ext = "anarci.exe" if self.platform == "windows" else "anarci"
        bundled = self.anarci_dir / ext
        if bundled.exists():
            return str(bundled)
        return shutil.which("anarci")
    
    def setup_path(self) -> None:
        """Add tool directories to os.environ['PATH']."""
        additions = []
        for d in [str(self.msms_dir), str(self.apbs_dir), str(self.pdb2pqr_dir), str(self.anarci_dir)]:
            if os.path.isdir(d) and d not in os.environ.get("PATH", ""):
                additions.append(d)
        if additions:
            os.environ["PATH"] = os.pathsep.join(additions + [os.environ.get("PATH", "")])
    
    def verify(self) -> dict:
        """Check that all tools are accessible. Returns dict of tool->bool."""
        def _available(attr: str) -> bool:
            val = getattr(self, attr)
            if val is None:
                return False
            return os.path.exists(val) or shutil.which(val) is not None

        return {
            "msms": _available("msms"),
            "apbs": _available("apbs"),
            "pdb2pqr": _available("pdb2pqr"),
            "anarci": self.anarci is not None,
            "platform": self.platform,
            "is_wsl": _is_wsl(),
            "is_windows": _is_windows(),
        }
    
    def __repr__(self) -> str:
        return (
            f"PlatformConfig(platform='{self.platform}', "
            f"msms='{self.msms}', apbs='{self.apbs}', pdb2pqr='{self.pdb2pqr}')"
        )


def _find_project_root() -> pathlib.Path:
    """Find the project root directory (where pyproject.toml lives)."""
    # Start from this file's location
    current = pathlib.Path(__file__).resolve().parent.parent
    # Walk up to find pyproject.toml
    for _ in range(5):
        if (current / "pyproject.toml").exists():
            return current
        current = current.parent
    # Fallback: current working directory
    return pathlib.Path.cwd()


# Singleton cache
_config_cache: PlatformConfig | None = None


def get_config() -> PlatformConfig:
    """Get the platform configuration (cached singleton)."""
    global _config_cache
    if _config_cache is None:
        if _is_windows():
            platform = "windows"
        elif _is_wsl():
            platform = "wsl"
        else:
            platform = "linux"
        
        root = _find_project_root()
        tools_dir = root / "Tools"
        
        _config_cache = PlatformConfig(
            platform=platform,
            project_root=root,
            tools_dir=tools_dir,
        )
    
    return _config_cache


def detect_platform() -> dict:
    """Quick platform detection, useful for shell scripts."""
    cfg = get_config()
    return cfg.verify()


# Allow running as a script: python -m surface_analyses.platform_config
if __name__ == "__main__":
    cfg = get_config()
    print(f"Platform: {cfg.platform}")
    print(f"Project root: {cfg.project_root}")
    print(f"msms: {cfg.msms}")
    print(f"apbs: {cfg.apbs}")
    print(f"pdb2pqr: {cfg.pdb2pqr}")
    print(f"default_apbs_work_dir: {cfg.default_apbs_work_dir}")
    print()
    print("Verification:")
    for tool, available in cfg.verify().items():
        status = "✓" if available else "✗"
        print(f"  {status} {tool}")
