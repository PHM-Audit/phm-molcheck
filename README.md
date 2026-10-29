# phm-molcheck

Deterministic **physical-validity gate for AI-generated molecular and
protein structures** (PDB-level). Pure checks, no learnable parameters,
no RDKit dependency — just `numpy`.

The checker does not trust the generating model's internal state; it
reads only the output structure file and answers one question: *is this
structure physically plausible?*

Battle-tested as an independent QC layer on AI-generated protein-binder
structures; generalized here as a model-agnostic tool.

## Checks (all deterministic)

| ID | Check | Rule |
|----|-------|------|
| C1 | atomic clashes | non-bonded pair with overlap `r_vdw_i + r_vdw_j - d >= 0.4 A` (MolProbity-style, grid-accelerated) |
| C2 | bond-length sanity | covalent distance vs standard covalent radii, deviation > 25% flagged |
| C3 | cross clashes | generated structure vs a reference (ligand/target), when supplied |
| C4 | packing density | atoms per 1000 A^3 (observational) |

Parser hygiene: crystal waters (HOH/WAT), altloc B+ conformations, and
hydrogens are filtered — the classic false-clash sources in PDB files.

Every input file gets a SHA-256 digest; every batch report carries an
audit digest. Same input → same report, always.

## Install

```bash
pip install numpy
pip install git+https://github.com/PHM-Audit/phm-molcheck
```

## Use as a library

```python
from phm_molcheck import qc_structure, batch_qc

r = qc_structure("generated_binder.pdb", ref_pdb="target_complex.pdb")
print(r["verdict"])           # PASS / FAIL / NO_ATOMS
print(r["C1_clash_pairs"])    # number of severe non-bonded overlaps
print(r["C2_bad_bonds"])      # unreasonable covalent distances

rep = batch_qc("./out/**/*.pdb", ref_pdb="target_complex.pdb",
               out_json="report.json")
print(rep["summary"])         # {'PASS': .., 'FAIL': .., 'pass_rate': ..}
```

## Use as a CLI

```bash
python -m phm_molcheck.checker "./generated/*.pdb" target.pdb report.json
```

## Scope and limits

- This package contains deterministic PDB-level geometric checks only.
  It does not perform force-field minimization, learned quality
  prediction, or chemistry-aware graph validation (that would require
  RDKit or heavier machinery).
- Verdicts are mechanical: C1/C3 clashes drive FAIL; C4 is observational.
- Hydrogens are ignored (crystal structures typically omit them).

## Tests

```bash
pip install pytest numpy
pytest tests/ -q
```

## License

MIT
