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
    assert len(benzene.ports) == 6
    assert set(benzene.ports.values()) == set(range(6))
    assert to_smiles(benzene) == "c1ccccc1"


def test_expanded_library_structures_finalize_with_expected_formulas():
    import molcule_desc as md

    expected = {
        "methane": "CH4",
        "ethane": "C2H6",
        "propane": "C3H8",
        "acetylene": "C2H2",
        "cyclohexane": "C6H12",
        "ethanol": "C2H6O",
        "pyridine": "C5H5N",
        "pyrrole": "C4H5N",
        "thiophene": "C4H4S",
        "imidazole": "C3H4N2",
        "naphthalene": "C10H8",
        "biphenyl": "C12H10",
        "diphenyl_ether": "C12H10O",
    }
    for name, formula in expected.items():
        assert getattr(md, name).finalize().formula_string() == formula


def test_expanded_library_smiles_for_common_structures():
    pytest.importorskip("rdkit")
    import molcule_desc as md

    expected = {
        "ethane": "CC",
        "propane": "CCC",
        "acetylene": "C#C",
        "cyclohexane": "C1CCCCC1",
        "pyridine": "c1ccncc1",
        "naphthalene": "c1ccc2ccccc2c1",
    }
    for name, smiles in expected.items():
        assert to_smiles(getattr(md, name)) == smiles


def test_library_fragments_expose_remaining_valences_as_ports():
    import molcule_desc as md

    assert len(md.methane.ports) == 4
    assert len(md.ethane.ports) == 6
    assert len(md.benzene.ports) == 6
    assert len(md.hydroxyl.ports) == 2
    assert len(md.amino.ports) == 3
    assert len(md.phosphine.ports) == 3
    assert len(md.phosphate_ester.ports) == 3



def test_tert_butyl_library_fragment_is_built_from_methane_units():
    import molcule_desc as md

    assert len(md.tert_butyl.ports) == 1
    # C4H9 is the tert-butyl substituent; its single port attaches to a parent.
    assert md.tert_butyl.formula == {"C": 4, "H": 9}



def test_substitute_replaces_multiple_host_ports_and_preserves_guest_ports():
    from molcule_desc import methane, substitute

    result = substitute(
        benzene,
        {
            single(1): (methane, single(1)),
            single(4): (methane, single(1)),
        },
    )
    assert len(result.atoms) == 8
    # Four benzene ports remain, plus three exposed ports from each methane.
    assert len(result.ports) == 10
    assert result.formula == {"C": 8}


def test_substitute_rejects_unknown_host_port():
    from molcule_desc import methane, substitute

    with pytest.raises(ValueError, match="Host has no exposed ports"):
        substitute(benzene, {single(99): (methane, single(1))})



def test_substitute_allows_explicit_expose_mapping():
    from molcule_desc import methane, substitute

    result = substitute(
        benzene,
        {
            single(1): (methane, single(1)),
            single(4): (methane, single(1)),
        },
        expose={
            single(1): (0, single(2)),
            single(2): (1, single(2)),
        },
    )
    assert set(result.ports) == {single(1), single(2)}
    assert len(result.ports) == 2



def test_library_supports_elements_with_multiple_allowed_valences():
    import molcule_desc as md

    # Explicitly selected higher-valence states retain the intended ports.
    assert len(md.phosphine.ports) == 3
    assert len(md.phosphine_oxide.ports) == 3
    assert len(md.sulfone.ports) == 2







def test_binol_is_composed_from_naphthol_and_exposes_remaining_ports():
    import molcule_desc as md

    assert md.binol.name == "BINOL"
    assert md.binol.finalize().formula_string() == "C20H14O2"
    # Two naphthol fragments each expose seven ring ports plus one O-H port;
    # one ring port on each fragment is consumed by the biaryl bond.
    assert len(md.binol.ports) == 12
    assert md.binol.finalize().ports == {}


def test_binol_converts_to_valid_smiles():
    pytest.importorskip("rdkit")
    import molcule_desc as md

    assert to_smiles(md.binol)
