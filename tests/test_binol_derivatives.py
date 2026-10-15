import pytest

from molcule_desc import to_smiles
from examples import binol_derivatives as derivatives


def test_example_binol_derivatives_have_expected_formulas():
    expected = {
        "binol_3_3_dibromo": "C20H12Br2O2",
        "binol_3_3_diiodo": "C20H12I2O2",
        "binol_6_6_dibromo": "C20H12Br2O2",
        "binol_3_3_dimethyl": "C22H18O2",
        "binol_2_2_diformyl": "C22H14O2",
        "binol_bis_2_pyridyl_imine": "C32H22N4",
        "binol_bis_2_pyrimidyl_imine": "C30H20N6",
        "binol_6_6_di_tert_butyl": "C28H30O2",
        "binol_6_6_diphenyl": "C32H22O2",
    }
    for name, formula in expected.items():
        molecule = getattr(derivatives, name).finalize()
        assert molecule.formula_string() == formula
        assert molecule.ports == {}


def test_example_binol_derivatives_convert_to_valid_smiles():
    pytest.importorskip("rdkit")
    for name in (
        "binol_3_3_dibromo", "binol_3_3_diiodo", "binol_6_6_dibromo",
        "binol_3_3_dimethyl", "binol_2_2_diformyl",
        "binol_bis_2_pyridyl_imine", "binol_bis_2_pyrimidyl_imine",
        "binol_6_6_di_tert_butyl", "binol_6_6_diphenyl",
    ):
        assert to_smiles(getattr(derivatives, name))
