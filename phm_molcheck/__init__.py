# -*- coding: utf-8 -*-
"""phm-molcheck public API.

Deterministic physical-validity checks for AI-generated molecular /
protein structures (PDB-level):
    parse_pdb_atoms(path)   -- PDB parser (filters crystal waters, altloc B+, hydrogens)
    grid_clash(atoms)       -- C1 clash detection (MolProbity-style, grid-accelerated)
    bond_check(atoms)       -- C2 bond-length sanity
    qc_structure(path, ref) -- single-structure verdict (C1/C2/C4, +C3 with ref)
    batch_qc(pattern, ...)  -- batch report with SHA-256 audit digests
    VDW, COV, BOND_TOL, CLASH_FACTOR, GRID
"""

from .checker import (parse_pdb_atoms, grid_clash, bond_check,
                      qc_structure, batch_qc,
                      VDW, COV, BOND_TOL, CLASH_FACTOR, GRID)

__version__ = "1.0.0"
__all__ = ["parse_pdb_atoms", "grid_clash", "bond_check", "qc_structure",
           "batch_qc", "VDW", "COV", "BOND_TOL", "CLASH_FACTOR", "GRID",
           "__version__"]
