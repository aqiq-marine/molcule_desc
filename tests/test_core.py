import pytest

from molcule_desc import (
    Atom,
    Bond,
    BondType,
    PortRef,
    Structure,
    hydrogen_single,
    carbon_double_two_single,
    carbon_triple_single,
    oxygen_double,
    acetone,
    benzene,
    atom,
    ported_atom,
    bond,
    ring,
    carbonyl,
    double,
    ethylene,
    formaldehyde,
    methanimine,
    single,
    triple,
    to_smiles,
)


def test_atom_topology_templates_have_expected_ports():
    assert set(hydrogen_single.ports) == {single(1)}
    assert set(carbon_double_two_single.ports) == {single(1), single(2), double(1)}
    assert set(carbon_triple_single.ports) == {single(1), triple(1)}
    assert set(oxygen_double.ports) == {double(1)}


def test_atom_declares_element_and_valence():
    carbon = atom("C", 4)
    assert carbon.element == "C"
    assert carbon.valence == 4


def test_ported_atom_wraps_an_atom_with_attachment_ports():
    carbon = atom("C", 4)
    template = ported_atom(carbon, [double(1), single(1), single(2)])
    assert template.atoms[0].element == "C"
    assert template.atoms[0].valence == 4
    assert set(template.ports) == {double(1), single(1), single(2)}


def test_ethylene_is_built_by_connecting_two_carbon_topologies():
    molecule = ethylene.finalize()
    assert molecule.formula == {"C": 2, "H": 4}
    assert molecule.ports == {}
    assert any(edge.order is BondType.DOUBLE for edge in molecule.bonds)


def test_methanimine_is_built_by_connecting_c_and_n_double_ports():
    molecule = methanimine.finalize()
    assert molecule.formula == {"C": 1, "N": 1, "H": 3}
    assert any(edge.order is BondType.DOUBLE for edge in molecule.bonds)


def test_formaldehyde_is_built_by_connecting_carbon_and_oxygen_topologies():
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
    left = ported_atom(atom("C", 4), [single(1)])
    right = ported_atom(atom("O", 2), [double(1)])
    with pytest.raises(ValueError, match="Bond types do not match"):
        bond(
            left, single(1), right, double(1),
            expose={},
        )


def test_consumed_port_cannot_be_published():
    left = ported_atom(atom("C", 4), [single(1), single(2)])
    right = ported_atom(atom("C", 4), [single(1)])
    with pytest.raises(ValueError, match="Cannot expose consumed port"):
        bond(
            left, single(1), right, single(1),
            expose={single(1): PortRef(0, single(1))},
        )


def test_terminal_templates_finalize_to_expected_formulas():
    assert hydrogen_single.finalize().formula == {"H": 1}
    assert carbon_triple_single.finalize().formula == {"C": 2, "H": 2}
    assert carbon_double_two_single.finalize().formula == {"C": 2, "H": 4}


def test_ring_closes_ports_without_duplicating_graph():
    chain = Structure(
        "chain",
        [Atom(0, "C", 4), Atom(1, "C", 4)],
        [Bond(0, 1, BondType.SINGLE)],
        {single(1): 0, single(2): 1},
    )
    closed = ring(
        chain, single(1), single(2),
        expose={},
    )
    assert len(closed.atoms) == 2
    assert len(closed.bonds) == 2
    assert closed.ports == {}



def test_ring_closes_two_ports_on_the_same_structure():
    chain = Structure(
        "chain",
        [Atom(0, "C", 4), Atom(1, "C", 4), Atom(2, "C", 4)],
        [
            Bond(0, 1, BondType.SINGLE),
            Bond(1, 2, BondType.SINGLE),
        ],
        {single(1): 0, single(2): 2},
    )
    cyclic = ring(chain, single(1), single(2), expose={})
    assert len(cyclic.atoms) == 3
    assert len(cyclic.bonds) == 3
    assert cyclic.ports == {}



def test_ring_exports_remaining_ports_without_port_refs():
    chain = Structure(
        "chain",
        [Atom(0, "C", 4), Atom(1, "C", 4), Atom(2, "C", 4)],
        [
            Bond(0, 1, BondType.SINGLE),
            Bond(1, 2, BondType.SINGLE),
        ],
        {single(1): 0, single(2): 2, single(3): 1},
    )
    cyclic = ring(
        chain,
        single(1),
        single(2),
        expose={single(1): single(3)},
    )
    assert cyclic.ports == {single(1): 1}



def test_rdkit_smiles_conversion_for_library_molecules():
    pytest.importorskip("rdkit")
    assert to_smiles(ethylene) == "C=C"
    assert to_smiles(formaldehyde) == "C=O"
    assert to_smiles(acetone) == "CC(C)=O"


def test_rdkit_smiles_conversion_requires_finalized_structure_when_requested():
    pytest.importorskip("rdkit")
    template = ported_atom(atom("C", 4), [single(1)])
    with pytest.raises(ValueError, match="exposed ports"):
        to_smiles(template, finalize=False)


def test_benzene_smiles_is_aromatic_benzene():
    pytest.importorskip("rdkit")
    assert benzene.formula == {"C": 6}
    assert len(benzene.bonds) == 6
    assert to_smiles(benzene) == "c1ccccc1"
