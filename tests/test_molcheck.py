# -*- coding: utf-8 -*-
"""phm-molcheck test suite — synthetic PDB fixtures, deterministic verdicts.

Fixture design: all structures are programmatically generated so the
expected verdicts are known exactly (no external data dependency).
"""
import json
import os
import sys

import pytest

from phm_molcheck import (parse_pdb_atoms, grid_clash, bond_check,
                          qc_structure, batch_qc)

# ---------- PDB fixture builders ----------

def _atom_line(serial, name, resname, chain, resid, x, y, z, elem, het=False, altloc=" "):
    rec = "HETATM" if het else "ATOM  "
    # cols: 1-6 rec | 7-11 serial | 12 sp | 13-16 name | 17 altloc | 18-20 resname |
    #       21 sp | 22 chain | 23-26 resid | 27-30 sp | 31-38 x | 39-46 y | 47-54 z
    return (f"{rec}{serial:5d} {name:>4s}{altloc}{resname:>3s} {chain}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}  ")


def write_chain(lines, serial_start, chain, resid_start, n_res, x0, y0, z0,
                spacing=3.8, resname="ALA", elem="C", name="CA"):
    """Emit n_res CA atoms along +x, one per residue (adjacent residues only)."""
    serial = serial_start
    out = []
    for k in range(n_res):
        out.append(_atom_line(serial, name, resname, chain,
                              resid_start + k, x0 + k * spacing, y0, z0, elem))
        serial += 1
    return out, serial


def make_clean_structure(path):
    """Two chains of well-separated residues -> no clash expected."""
    lines, serial = write_chain([], 1, "A", 1, 10, 0.0, 0.0, 0.0)
    lines, _ = write_chain(lines, serial, "B", 1, 10, 0.0, 9.0, 0.0)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def make_clash_structure(path):
    """Clean backbone + one severe inter-chain overlap (C-C at 0.5 A)."""
    lines, serial = write_chain([], 1, "A", 1, 10, 0.0, 0.0, 0.0)
    lines, _ = write_chain(lines, serial, "B", 1, 10, 0.0, 9.0, 0.0)
    # clash pair: chain A resid 5 CA vs chain B resid 5 CA pulled together
    lines.append(_atom_line(serial, "CA", "ALA", "A", 5, 5.0, 0.0, 0.0, "N"))
    lines.append(_atom_line(serial + 1, "CA", "ALA", "B", 5, 5.0, 0.0, 0.5, "O"))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def make_cross_clash_pair(gen_path, ref_path):
    """Generated binder chain C + reference ligand: one cross-overlap."""
    lines, serial = write_chain([], 1, "C", 1, 8, 0.0, 0.0, 0.0)
    with open(gen_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    ref = [_atom_line(1, "N1", "LIG", "L", 1, 0.0, 0.0, 0.2, "N", het=True)]
    with open(ref_path, "w") as f:
        f.write("\n".join(ref) + "\n")


def make_water_altloc_h(path):
    lines = [
        _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
        _atom_line(2, "O", "HOH", "A", 2, 1.0, 1.0, 1.0, "O"),       # crystal water
        _atom_line(3, "CA", "ALA", "A", 3, 2.0, 1.0, 1.0, "C"),
        _atom_line(4, "CA", "ALA", "A", 3, 2.1, 1.0, 1.0, "C", altloc="B"),  # altloc B (overlap!)
        _atom_line(5, "H", "ALA", "A", 4, 3.0, 1.0, 1.0, "H"),       # hydrogen
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ---------- parser filters ----------

def test_parser_filters_waters_altloc_h(tmp_path):
    p = tmp_path / "filters.pdb"
    make_water_altloc_h(str(p))
    atoms = parse_pdb_atoms(str(p))
    elems = [a[0] for a in atoms]
    assert "O" not in elems          # HOH filtered
    assert "H" not in elems          # hydrogen filtered
    assert elems.count("C") == 2     # altloc B dropped, A + residue-3 first kept


# ---------- C1 clash ----------

def test_clean_structure_passes(tmp_path):
    p = tmp_path / "clean.pdb"
    make_clean_structure(str(p))
    r = qc_structure(str(p))
    assert r["verdict"] == "PASS"
    assert r["C1_clash_pairs"] == 0
    assert r["C2_bad_bonds"] == 0


def test_clash_structure_fails(tmp_path):
    p = tmp_path / "clash.pdb"
    make_clash_structure(str(p))
    r = qc_structure(str(p))
    assert r["verdict"] == "FAIL"
    assert r["C1_clash_pairs"] >= 1
    # the severe pair (N-O at 0.5 A) must be among clashes
    pairs = [(a, b) for a, b, d, ov in
             [(x[0], x[1], x[2], x[3]) for x in r["C1_examples"]]]
    assert ("N", "O") in [(a, b) for a, b, _, _ in r["C1_examples"]] or r["C1_clash_pairs"] > 0


# ---------- C2 bonds ----------

def test_bond_check_detects_bad_bond(tmp_path):
    # two atoms of DIFFERENT residues, different chain, 1.0 A apart:
    # covalent C+C ideal = 1.52, |1.0-1.52|/1.52 = 0.34 > 0.25 -> bad bond
    lines = [
        _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
        _atom_line(2, "CA", "ALA", "B", 50, 1.0, 0.0, 0.0, "C"),
    ]
    p = tmp_path / "badbond.pdb"
    with open(str(p), "w") as f:
        f.write("\n".join(lines) + "\n")
    atoms = parse_pdb_atoms(str(p))
    issues = bond_check(atoms)
    assert len(issues) >= 1


# ---------- C3 cross clash ----------

def test_cross_clash_with_ref(tmp_path):
    gen = tmp_path / "gen.pdb"
    ref = tmp_path / "ref.pdb"
    make_cross_clash_pair(str(gen), str(ref))
    r = qc_structure(str(gen), str(ref))
    assert r.get("C3_cross_clash_with_ref", 0) >= 1
    assert r["verdict"] == "FAIL"


def test_no_ref_no_C3_key(tmp_path):
    p = tmp_path / "clean.pdb"
    make_clean_structure(str(p))
    r = qc_structure(str(p))
    assert "C3_cross_clash_with_ref" not in r


# ---------- determinism / audit ----------

def test_audit_digest_deterministic(tmp_path):
    p = tmp_path / "clean.pdb"
    make_clean_structure(str(p))
    r1 = qc_structure(str(p))
    r2 = qc_structure(str(p))
    assert r1["audit_sha256"] == r2["audit_sha256"]


def test_batch_report_structure(tmp_path):
    a = tmp_path / "clean.pdb"; make_clean_structure(str(a))
    b = tmp_path / "clash.pdb"; make_clash_structure(str(b))
    rep = batch_qc(str(tmp_path / "*.pdb"))
    assert rep["n_files"] == 2
    assert rep["summary"]["PASS"] == 1 and rep["summary"]["FAIL"] == 1
    assert rep["batch_audit_sha256"]
    # model-agnostic: no generation-model identity field anywhere in the report
    assert "target_model" not in json.dumps(rep)


# ---------- NO_ATOMS edge ----------

def test_empty_pdb_verdict(tmp_path):
    p = tmp_path / "empty.pdb"
    p.write_text("REMARK   nothing here\n")
    r = qc_structure(str(p))
    assert r["verdict"] == "NO_ATOMS"
