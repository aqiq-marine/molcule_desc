import pytest

from molcule_desc import h8_binol, tetralol, to_smiles


def test_tetralol_fragment_formula():
    molecule = tetralol.finalize()
    assert molecule.formula_string() == "C10H12O"
    assert molecule.ports == {}


def test_h8_binol_formula_and_ports():
    molecule = h8_binol.finalize()
    assert molecule.formula_string() == "C20H22O2"
    assert molecule.ports == {}


def test_h8_binol_converts_to_smiles():
    pytest.importorskip("rdkit")
    assert to_smiles(h8_binol)
