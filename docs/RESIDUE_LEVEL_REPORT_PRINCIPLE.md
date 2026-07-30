# PEP-Patch 残基级报告原理

## 概述

PEP-Patch（`surface_analyses`）的核心计算在**几何顶点级别**操作：它生成蛋白质的溶剂可及分子表面（三角形网格），将静电势或疏水性数值映射到每个顶点，然后按连通性分割为离散的 surface patch。但抗体工程和可开发性分析需要的是**残基级别**的信息——哪些残基参与了哪个 patch、贡献了多少面积。

本文档详细说明如何从顶点级数据转换为残基级报告。

---

## 1. 分子表面与 Patch 分割

### 1.1 表面生成

使用 MSMS（Michel Sanner's Molecular Surface）或 Marching Cubes 算法生成溶剂排除表面（SES），输出为三角形网格：

```
顶点: (x₁, y₁, z₁), (x₂, y₂, z₂), ..., (xₙ, yₙ, zₙ)
三角面: (vᵢ, vⱼ, vₖ), ...
```

### 1.2 属性映射

每个顶点被赋予一个标量值：

- **静电势（ES）**：APBS 求解 Poisson-Boltzmann 方程 → 三维电位网格 → 插值到表面顶点（单位：kT/e）
- **疏水性（HB）**：Crippen logP 原子倾向性 → Heiden 距离加权映射到表面顶点（无量纲）

### 1.3 Patch 分割

按连通性将顶点分组。两个相邻顶点属于同一 patch 当且仅当：
- 它们通过三角面的边直接相邻
- 它们的标量值符号相同（ES）或跨越同一分段（HB）

每个顶点被分配一个 patch ID（0, 1, 2, ...），顶点颜色在 PLY 文件中编码此 ID。

### 1.4 顶点面积计算

对每个三角面 $(v_a, v_b, v_c)$：

$$\text{Area}_{\triangle} = \frac{1}{2} \| \vec{ab} \times \vec{ac} \|$$

$$\text{Area}_{v_i} = \frac{1}{3} \text{Area}_{\triangle}$$

每个顶点的面积为共享该顶点的所有三角面贡献之和。NPZ 模式下使用 `np.add.at` 累加；PLY 模式下直接读取面索引数组。

---

## 2. 顶点 → 残基映射

这是将几何表面数据转换为残基级信息的核心步骤。

### 2.1 PDB 解析（ES 模式）

```
PDB 文件
  │
  ├─ mdtraj 加载 → 原子拓扑 + 原子坐标 (Å)
  │
  └─ Bio.PDB 解析 → (chain_id, res_name, res_seq) 三元组
```

合并两者构建 `atom_df`：

| atom_idx | atom_name | chain_id | res_name | res_seq | res_id | seq_nr | seq_res_id |
|----------|-----------|----------|----------|---------|--------|--------|------------|
| 0 | N | H | GLY | 26 | GLY26 | 1 | GLY1 |
| 1 | CA | H | GLY | 26 | GLY26 | 1 | GLY1 |
| ... | ... | ... | ... | ... | ... | ... | ... |

`seq_nr`（顺序编号）将 PDB 中不连续的残基编号（含插入码）映射为连续的整数序列，便于后续输出。

### 2.2 最近邻映射（ES 模式）

使用 **KD-Tree**（k-d 树）空间索引：

```python
from scipy.spatial import KDTree
tree = KDTree(atom_coords)  # 对所有原子坐标建索引

# 对每个 patch 顶点，查询最近原子
dists, atom_idxs = tree.query(patch_vertices, k=1)
```

**物理原理：** patch 顶点位于溶剂可及表面上，距离最近的原子通常是支配该表面区域化学性质的残基原子。顶点面积归入该原子所属残基。

### 2.3 内建映射（HB 模式）

HB 的 NPZ 文件在生成时已经记录了每个表面顶点对应的原子索引：

```python
atom_arr = npz[f"{basename}:data:atom"]  # 每个顶点 → 原子索引
```

然后通过 PDB 解析的 `atom_map` 将原子索引转换为 `(chain_id, res_name, res_seq)`：

```python
chain, rname, rseq, gseq, rid = atom_map[atom_index]
```

### 2.4 按残基聚合

将同一 patch 内归属同一残基的顶点面积求和：

```
patch_id=3, patch_type="positive", patch_total_area=2279.22 Å²

顶点 i₁ → atom_idx 145 → HIS-H-52 → 1.2 Å²
顶点 i₂ → atom_idx 146 → HIS-H-52 → 0.8 Å²
顶点 i₃ → atom_idx 147 → HIS-H-52 → 2.1 Å²
...                    → HIS-H-52 → 总计 4.1 Å², frac=0.0018

顶点 iⱼ → atom_idx 200 → ARG-H-98 → 15.3 Å²
...                    → ARG-H-98 → 总计 89.7 Å², frac=0.0393
```

---

## 3. 输出格式

### 3.1 `{stem}_residues_detailed.csv`

每行 = 一个残基在一个 patch 中的占比信息：

| 列 | 含义 | 单位 |
|-----|------|------|
| `patch_nr` | Patch 编号 | — |
| `patch_type` | 类型（positive / negative / hydrophobic / hydrophilic） | — |
| `patch_total_area_A2` | 该 patch 总面积 | Å² |
| `chain_id` | 链标识 | — |
| `res_name` | 残基三字母名 | — |
| `res_seq` | PDB 残基编号 | — |
| `res_id` | 残基标识（如 `HIS52`） | — |
| `seq_nr` | 顺序编号（连续） | — |
| `seq_res_id` | 顺序残基标识 | — |
| `n_vertices` | 该残基在该 patch 中的顶点数 | — |
| `area_A2` | 该残基在该 patch 中的贡献面积 | Å² |
| `frac_of_patch` | 占该 patch 总面积的比例 | — |
| `mean_dist_A` | 该 patch 顶点到质心的平均距离 | Å |

### 3.2 `{stem}_patch_summary.csv`

每个 patch 一行，包含总面积、覆盖残基数、最大贡献残基等汇总信息。

---

## 4. 精度与局限性

### 4.1 精度

- **顶点面积**：精确（三角形几何直接计算）
- **顶点 → 残基映射**：最近邻匹配精度受限于网格分辨率（默认 0.05 nm = 0.5 Å）。在 patch 内部，顶点通常距其所属原子的 van der Waals 表面仅 1.4 Å（溶剂探针半径），最近邻可靠
- **patch 边界**：边界顶点可能被映射到相邻残基，引入微小的面积归属误差（< 2%）

### 4.2 局限性

1. **不区分主链/侧链**：顶点映射到原子，聚合到残基，但不区分原子类型
2. **溶剂暴露不归一化**：`area_A2` 是绝对面积，未除以该残基类型的参考 SASA
3. **依赖网格分辨率**：更细的网格（`--grid_spacing 0.03`）产生更精确的映射，但计算量和文件大小增加
4. **对称性误差**：完全对称的残基对（如两个对称的 ARG）可能因浮点精度造成微小面积差异

---

## 5. 数据流总览

```
                  PDB 文件
                     │
         ┌───────────┴───────────┐
         ▼                       ▼
   ES 分析 (APBS)          HB 分析 (Crippen)
         │                       │
   静电势网格                Crippen logP 原子倾向性
         │                       │
   MSMS 分子表面            Marching Cubes 分子表面
         │                       │
   顶点势能映射              Heiden 距离加权映射
         │                       │
   连通分量分割               patch 分割
         │                       │
   _es_patches.csv          _hb_out.npz
   _es-pos.ply              (顶点 + patch_id + atom_idx)
   _es-neg.ply
         │                       │
         └───────────┬───────────┘
                     ▼
           unified_analyzer.py
                     │
         ┌───────────┴───────────┐
         ▼                       ▼
   KD-Tree 最近邻映射       atom_map 直接查找
   (PLY 顶点 → PDB 原子)     (NPZ 原子索引 → 残基)
         │                       │
         └───────────┬───────────┘
                     ▼
              按残基聚合面积
                     │
                     ▼
         {stem}_residues_detailed.csv
         {stem}_patch_summary.csv
```
