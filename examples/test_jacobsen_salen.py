import pytest

from molcule_desc import to_smiles
from examples.jacobsen_salen import build_jacobsen_salen


REFERENCE_SMILES = (
    "CC(C)(C)C1=CC(C(C)(C)C)=C(O)C(/C=N/C2CCCCC2/N=C/"
    "C3=C(O)C(C(C)(C)C)=CC(C(C)(C)C)=C3)=C1"
)


def test_jacobsen_salen_connectivity_formula():
    ligand = build_jacobsen_salen()
    assert ligand.finalize().formula_string() == "C36H54N2O2"


def test_jacobsen_salen_matches_reference_smiles_connectivity():
    rdkit = pytest.importorskip("rdkit")
    from rdkit import Chem

    generated = Chem.MolFromSmiles(to_smiles(build_jacobsen_salen()))
    reference = Chem.MolFromSmiles(REFERENCE_SMILES)
    assert generated is not None
    assert reference is not None

    # The current graph model does not encode E/Z or R/S stereochemistry.
    # Compare canonical non-isomeric SMILES to validate connectivity only.
    generated_connectivity = Chem.MolToSmiles(
        generated, canonical=True, isomericSmiles=False
    )
    reference_connectivity = Chem.MolToSmiles(
        reference, canonical=True, isomericSmiles=False
    )
    assert generated_connectivity == reference_connectivity
