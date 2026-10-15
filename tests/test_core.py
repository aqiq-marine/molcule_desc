import pytest

from molcule_desc import (
    Atom,
    Bond,
    BondType,
    PortRef,
    Structure,
    XS,
    XD,
    XT,
    XO,
    acetone,
    bond,
    carbonyl,
    double,
    ethylene,
    formaldehyde,
    methanimine,
    single,
    triple,
)


def test_dummy_templates_have_expected_ports():
    assert set(XS.ports) == {single(1)}
    assert set(XD.ports) == {single(1), single(2), double(1)}
    assert set(XT.ports) == {single(1), triple(1)}
    assert set(XO.ports) == {double(1)}


def test_ethylene_is_built_by_connecting_two_xd_fragments():
    molecule = ethylene.finalize()
    assert molecule.formula == {"C": 2, "H": 4}
    assert molecule.ports == {}
    assert any(edge.order is BondType.DOUBLE for edge in molecule.bonds)


def test_methanimine_is_built_by_connecting_c_and_n_double_ports():
    molecule = methanimine.finalize()
    assert molecule.formula == {"C": 1, "N": 1, "H": 3}
    assert any(edge.order is BondType.DOUBLE for edge in molecule.bonds)


def test_formaldehyde_is_built_by_connecting_xd_and_xo():
    molecule = formaldehyde.finalize()
    assert molecule.formula == {"C": 1, "H": 2, "O": 1}
    assert any(edge.order is BondType.DOUBLE for edge in molecule.bonds)


def test_acetone_is_built_from_carbonyl_and_methyl_fragments():
    molecule = acetone.finalize()
    assert molecule.formula == {"C": 3, "H": 6, "O": 1}
    assert any(edge.order is BondType.DOUBLE for edge in molecule.bonds)


def test_carbonyl_fragment_exposes_two_single_ports():
    assert set(carbonyl.ports) == {single(1), single(2)}


def test_mismatched_bond_types_are_rejected():
    left = Structure.atom("C", [single(1)])
    right = Structure.atom("O", [double(1)])
    with pytest.raises(ValueError, match="Bond types do not match"):
        bond(
            left, single(1), right, double(1),
            expose={},
        )


def test_consumed_port_cannot_be_published():
    left = Structure.atom("C", [single(1), single(2)])
    right = Structure.atom("C", [single(1)])
    with pytest.raises(ValueError, match="Cannot expose consumed port"):
        bond(
            left, single(1), right, single(1),
            expose={single(1): PortRef(0, single(1))},
        )


def test_xt_finalizes_to_ethyne_side_and_xs_to_hydrogen():
    assert XS.finalize().formula == {"H": 1}
    assert XT.finalize().formula == {"C": 2, "H": 2}
    assert XD.finalize().formula == {"C": 2, "H": 4}


def test_intramolecular_bond_closes_ports_without_duplicating_graph():
    chain = Structure(
        "chain",
        [Atom(0, "C"), Atom(1, "C")],
        [Bond(0, 1, BondType.SINGLE)],
        {single(1): 0, single(2): 1},
    )
    closed = bond(
        chain, single(1), chain, single(2),
        expose={},
        intramolecular=True,
    )
    assert len(closed.atoms) == 2
    assert len(closed.bonds) == 2
    assert closed.ports == {}
