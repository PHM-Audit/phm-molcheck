# phm-molcheck

Deterministic **physical-validity gate for AI-generated molecular and
protein structures** (PDB-level). Pure checks, no learnable parameters,
no RDKit dependency — just `numpy`.

The checker does not trust the generating model's internal state; it
reads only the output structure file and answers one question: *does
this structure contain obvious physical implausibilities?*

**These checks are necessary, not sufficient: a PASS verdict does not
imply chemical correctness and does not replace experimental validation.**

Originally built as an independent QC layer for AI-generated
protein-binder structures; generalized here as a model-agnostic tool.
This package contains only the deterministic PDB-level subset of that
QC layer.

## Checks (all deterministic)

| ID | Check | Rule |
|----|-------|------|
| C1 | atomic clashes | non-bonded pair with overlap `r_vdw_i + r_vdw_j - d >= 0.4 A` — inspired by the MolProbity clashscore concept (Davis et al., 2007); thresholds are our own; not affiliated with or endorsed by the MolProbity project |
| C2 | bond-length sanity | covalent distance vs standard covalent radii, deviation > 25% flagged |
| C3 | cross clashes | generated structure vs a reference (ligand/target), when supplied |
| C4 | packing density | atoms per 1000 A^3 (observational) |

Parser hygiene: crystal waters (HOH/WAT), altloc B+ conformations, and
hydrogens are filtered — the classic false-clash sources in PDB files.

Every input file gets a SHA-256 digest; every batch report carries an
audit digest. Same input → same report, always.

## Install

**Option A — from the GitHub repo** (no package index needed):

```bash
pip install git+https://github.com/PHM-Audit/phm-molcheck
```

**Option B — from a Release wheel** (faster, no git clone):

```bash
pip install https://github.com/PHM-Audit/phm-molcheck/releases/download/v1.0.0/phm_molcheck-1.0.0-py3-none-any.whl
```

PyPI listing is planned; both options above work today.

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

## Industry context (2026-11)

AI-generated protein designs are now reaching wet-lab validation at
frontier scale: OpenAI's November 2026 preprint ("GPT-6 as a Universal
Foundation Model for Protein Biology") reports 35 of 322 GPT-6-designed
minibinders showing binding evidence in wet-lab tests across five
targets, with sequence design handled by specialized tools. This is
real progress — and also a concrete case for independent,
structure-level QC of AI-generated designs: wet-lab binding success
(10.9% overall) coexists with imperfect structural predictions, and
downstream users of such structures need deterministic, replayable
physical-validity checks of the actual artifact files, independent of
the generating model's self-report.

`phm-molcheck` is exactly that kind of check layer — model-agnostic,
deterministic, audit-logged. It does not evaluate binding likelihood
and does not replace experimental validation; it gates the artifact
file itself.

Source: OpenAI preprint & blog, Nov 24 2026 (35/322 minibinders,
5 targets, PD-L1/IL-7R-alpha/TrkA/VEGF-A/IL-17A). Recorded here as
industry context only; OpenAI is not affiliated with this project.
