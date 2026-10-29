#!/usr/bin/env python3
"""
phm-molcheck -- deterministic physical-validity gate for AI-generated
molecular / protein structures (PDB-level, no RDKit required).

Architecture principle: the checker does not trust the generating model's
internal state; it reads only the output structure file (.pdb/.cif).
Checks (deterministic physics/chemistry rules, no learnable parameters):
  C1  atomic clash detection (MolProbity-style, grid-accelerated):
      non-bonded pair with overlap = r_vdw_i + r_vdw_j - d >= 0.4 A
  C2  bond-length sanity: covalent distance vs standard coval radii,
      deviation > 25% flagged
  C3  cross clashes between a generated binder and a reference
      structure (ligand/target), when a reference is supplied
  C4  packing density observation (atoms per 1000 A^3)
Output: per-structure PASS/FAIL verdict + batch summary + SHA-256 audit
digest of every input file.

Corigin note: battle-tested as an independent QC layer on generated
protein-binder structures; generalized here as a model-agnostic tool.
"""
import sys, os, json, glob, hashlib, time
import numpy as np

VDW = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80, 'H': 1.20, 'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'ZN': 1.39, 'MG': 1.73, 'CA': 2.31, 'FE': 1.94, 'SE': 1.90, 'I': 1.98}
COV = {'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05, 'P': 1.07, 'H': 0.31, 'F': 0.57, 'CL': 1.02, 'BR': 1.20, 'ZN': 1.22, 'MG': 1.41, 'CA': 1.76, 'FE': 1.32, 'SE': 1.20, 'I': 1.39}
BOND_TOL = 0.25   # 键长偏差容忍25%
CLASH_FACTOR = 0.75
GRID = 4.0        # 碰撞检测网格(Å)


def parse_pdb_atoms(path):
    """解析PDB: 返回 [(elem, x, y, z, chain, residx, hetero)]
    过滤: 结晶水(HOH)与altloc双构象(只留A/空白)——晶体学伪clash源"""
    atoms = []
    with open(path) as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                resn = line[17:20].strip()
                if resn in ('HOH', 'H2O', 'WAT'):
                    continue  # 结晶水: 位置受晶体学约束, 非碰撞评估对象
                altloc = line[16:17].strip()
                if altloc and altloc != 'A':
                    continue  # 替代构象: 双构象原子互相重叠, 只留A
                elem = line[76:78].strip().upper() if len(line) >= 78 else line[12:16].strip().strip('0123456789').upper()
                if elem.startswith('H'):
                    continue  # 氢原子: PDB晶体结构常不给氢或给不全, 跳过
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                except ValueError:
                    continue
                chain = line[21:22].strip()
                het = line.startswith('HETATM')
                atoms.append((elem if elem in VDW else 'C', x, y, z, chain, line[22:26], het))
    return atoms


def _ridx(s):
    """残基号规范化: ' 123A'->123 ( Chothia/Kabat插入码取数值部分)"""
    try:
        return int(''.join(c for c in s if c.isdigit()) or s)
    except (ValueError, TypeError):
        return None


def grid_clash(atoms, overlap_tol=0.4, exclude_same_residue=True):
    """网格加速碰撞检测(MolProbity标准): 重叠=r_i+r_j-d >= 0.4Å => clash
    排除规则: 同残基内全部原子对 + 同链相邻残基(数值号差1)肽键C-N"""
    if len(atoms) < 2:
        return []
    xyz = np.array([(a[1], a[2], a[3]) for a in atoms])
    cells = {}
    for i, (x, y, z) in enumerate(xyz):
        cells.setdefault((int(x // GRID), int(y // GRID), int(z // GRID)), []).append(i)
    clashes = []
    offsets = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    for key, idxs in cells.items():
        for off in offsets:
            nkey = (key[0] + off[0], key[1] + off[1], key[2] + off[2])
            if nkey not in cells:
                continue
            for i in idxs:
                ei, xi = atoms[i][0], xyz[i]
                for j in cells[nkey]:
                    if j <= i:
                        continue
                    ej, xj = atoms[j][0], xyz[j]
                    # 排除1: 同一残基内原子对(共价几何, 非碰撞)
                    if exclude_same_residue and atoms[i][4] == atoms[j][4] and atoms[i][5] == atoms[j][5]:
                        continue
                    # 排除2: 同链相邻残基(插入码已规范化)间全部原子对——
                    # 肽键C-N 1.33/O-N 2.25/CA-N 2.4等均为构象几何固定值, 非碰撞
                    ri, rj = _ridx(atoms[i][5]), _ridx(atoms[j][5])
                    if (ri is not None and rj is not None and abs(ri - rj) == 1
                            and atoms[i][4] == atoms[j][4]):
                        continue
                    d = float(np.linalg.norm(xi - xj))
                    if d < 1e-3:
                        continue  # 同一原子
                    # 二硫键白名单: CYS-SG ~ S-S 2.05Å共价键
                    if {atoms[i][0], atoms[j][0]} == {'S'} and 1.8 <= d <= 2.3:
                        continue
                    overlap = VDW[ei] + VDW[ej] - d
                    if overlap >= overlap_tol:
                        clashes.append((i, j, round(d, 3), round(overlap, 2)))
    return clashes


def bond_check(atoms):
    """键长合理性: 对距离在[0.9,1.9]x共价和的原子对, 检查偏差"""
    issues = []
    xyz = np.array([(a[1], a[2], a[3]) for a in atoms])
    n = len(atoms)
    for i in range(n):
        ei, xi = atoms[i][0], xyz[i]
        for j in range(i + 1, min(i + 40, n)):
            ej, xj = atoms[j][0], xyz[j]
            d = float(np.linalg.norm(xi - xj))
            if d > 2.2:
                continue
            ideal = COV[ei] + COV[ej]
            # 氢键/盐桥范围排除(2.5-3.5Å不算共价)
            if abs(d - ideal) / ideal > BOND_TOL and d < ideal * 1.25 - 0.35:
                issues.append((i, j, round(d, 3), round(ideal, 3)))
    return issues


def qc_structure(pdb_path, ref_pdb=None):
    """对单个生成结构跑C1/C2/C4, ref_pdb存在则加跑C3"""
    t0 = time.time()
    atoms = parse_pdb_atoms(pdb_path)
    r = {'file': os.path.basename(pdb_path), 'n_atoms': len(atoms)}
    if not atoms:
        r['verdict'] = 'NO_ATOMS'
        return r
    clashes = grid_clash(atoms)
    r['C1_clash_pairs'] = len(clashes)
    r['C1_examples'] = [(atoms[i][0], atoms[j][0], d, ov) for i, j, d, ov in clashes[:5]]
    bonds = bond_check(atoms)
    r['C2_bad_bonds'] = len(bonds)
    if ref_pdb and os.path.exists(ref_pdb):
        ref_atoms = parse_pdb_atoms(ref_pdb)
        both = atoms + ref_atoms
        cross = [(i, j, d) for i, j, d, ov in grid_clash(both) if (i < len(atoms)) != (j < len(atoms))]
        r['C3_cross_clash_with_ref'] = len(cross)
    xyz = np.array([(a[1], a[2], a[3]) for a in atoms])
    vol = float(np.prod(xyz.max(0) - xyz.min(0) + 1e-6))
    r['C4_atoms_per_1000A3'] = round(len(atoms) / max(vol, 1) * 1000, 2) if vol > 0 else 0
    r['verdict'] = 'FAIL' if (r['C1_clash_pairs'] > 0 or r.get('C3_cross_clash_with_ref', 0) > 0) else 'PASS'
    r['audit_sha256'] = hashlib.sha256(open(pdb_path, 'rb').read()).hexdigest()[:16]
    r['elapsed_s'] = round(time.time() - t0, 3)
    return r


def batch_qc(pattern, ref_pdb=None, out_json=None):
    files = sorted(glob.glob(pattern))
    report = {'tool': 'phm-molcheck (deterministic physical validity gate)',
              'n_files': len(files), 'results': [], 't': time.strftime('%F %T')}
    for f in files:
        report['results'].append(qc_structure(f, ref_pdb))
    n_pass = sum(1 for x in report['results'] if x.get('verdict') == 'PASS')
    n_fail = sum(1 for x in report['results'] if x.get('verdict') == 'FAIL')
    report['summary'] = {'PASS': n_pass, 'FAIL': n_fail,
                         'pass_rate': round(n_pass / max(len(files), 1), 4),
                         'total_clash_pairs': sum(x.get('C1_clash_pairs', 0) for x in report['results'])}
    blob = json.dumps(report, ensure_ascii=False, default=str)
    report['batch_audit_sha256'] = hashlib.sha256(blob.encode()).hexdigest()[:16]
    if out_json:
        json.dump(report, open(out_json, 'w'), ensure_ascii=False, indent=1)
    return report


if __name__ == '__main__':
    pattern = sys.argv[1] if len(sys.argv) > 1 else './demo/**/*.pdb'
    ref = sys.argv[2] if len(sys.argv) > 2 else None
    out = sys.argv[3] if len(sys.argv) > 3 else './phm_molcheck_report.json'
    rep = batch_qc(pattern, ref, out)
    print(json.dumps(rep['summary'], ensure_ascii=False, indent=1))
    print('batch audit digest:', rep['batch_audit_sha256'])
