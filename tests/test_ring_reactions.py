import pytest

from molcule_desc import to_smiles
from examples import ring_reactions as examples


@pytest.mark.parametrize(
    ("name", "formula", "ring_size"),
    [
        ("beta_propiolactone", "C3H4O2", 4),
        ("gamma_butyrolactone", "C4H6O2", 5),
        ("delta_valerolactone", "C5H8O2", 6),
        ("epsilon_caprolactone", "C6H10O2", 7),
        ("seven_membered_lactone", "C7H12O2", 8),
    ],
)
def test_lactone_examples_have_expected_formulas_and_ring_sizes(name, formula, ring_size):
    molecule = getattr(examples, name).finalize()
    assert molecule.formula_string() == formula
    assert molecule.ports == {}
    assert len(molecule.atoms) == ring_size + 1  # includes carbonyl oxygen


def test_iodoetherification_example_has_expected_formula():
    molecule = examples.iodoetherification_5_exo.finalize()
    assert molecule.formula_string() == "C5H9IO"
    assert molecule.ports == {}


def test_ring_reaction_examples_convert_to_smiles():
    pytest.importorskip("rdkit")
    for name in (
        "beta_propiolactone",
        "gamma_butyrolactone",
        "delta_valerolactone",
        "epsilon_caprolactone",
        "seven_membered_lactone",
        "iodoetherification_5_exo",
    ):
        assert to_smiles(getattr(examples, name))
