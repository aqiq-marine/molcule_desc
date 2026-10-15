import pytest

from molcule_desc import to_smiles
from examples.jacobsen_salen import build_jacobsen_salen


def test_jacobsen_salen_connectivity_formula():
    ligand = build_jacobsen_salen()
    assert ligand.finalize().formula_string() == "C36H54N2O2"


def test_jacobsen_salen_smiles_is_valid():
    pytest.importorskip("rdkit")
    smiles = to_smiles(build_jacobsen_salen())
    assert "N" in smiles
    assert "O" in smiles
    assert smiles.count("C") > 0
