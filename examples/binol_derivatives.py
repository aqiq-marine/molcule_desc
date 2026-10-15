"""Example BINOL derivatives built compositionally from reusable fragments.

These structures are examples, not part of the public built-in library. The
atom indices below refer to the current two-naphthol composition.
"""
from molcule_desc.core import Port, Structure, single, substitute
from molcule_desc.library import (
    bromo, carbonyl, imine, iodo, methane, naphthalene, phenyl, pyridine,
    _graph_structure,
)
from molcule_desc.library_advanced import binol, tert_butyl, _port_at_atom


def _replace_with_fragment(host: Structure, anchor_atom: int, guest: Structure, *, name: str) -> Structure:
    """Attach a guest fragment while preserving all other host ports."""
    host_port = _port_at_atom(host, anchor_atom)
    guest_port = next(iter(guest.ports))
    expose = {
        single(index): (0, port)
        for index, port in enumerate((p for p in host.ports if p != host_port), start=1)
    }
    return substitute(host, {host_port: (guest, guest_port)}, expose=expose, name=name)


def _derivative(name: str, substituents: dict[int, str]) -> Structure:
    guests = {"Br": bromo, "I": iodo, "Me": methane, "tBu": tert_butyl, "phenyl": phenyl}
    result = binol
    for atom_id, substituent in substituents.items():
        result = _replace_with_fragment(result, atom_id, guests[substituent], name=name)
    result.name = name
    return result


binol_3_3_dibromo = _derivative("3,3'-dibromo-BINOL", {3: "Br", 14: "Br"})
binol_3_3_diiodo = _derivative("3,3'-diiodo-BINOL", {3: "I", 14: "I"})
binol_6_6_dibromo = _derivative("6,6'-dibromo-BINOL", {7: "Br", 18: "Br"})
binol_3_3_dimethyl = _derivative("3,3'-dimethyl-BINOL", {3: "Me", 14: "Me"})
binol_6_6_di_tert_butyl = _derivative("6,6'-di-tert-butyl-BINOL", {7: "tBu", 18: "tBu"})
binol_6_6_diphenyl = _derivative("6,6'-diphenyl-BINOL", {7: "phenyl", 18: "phenyl"})


def _compose_binaphthyl_base() -> Structure:
    """Join two naphthalene fragments at their 1-position ports."""
    port = _port_at_atom(naphthalene, 1)
    return substitute(naphthalene, {port: (naphthalene, port)}, name="1,1'-binaphthyl")


def _attach_preserving_guest_port(
    host: Structure, host_anchor: int, guest: Structure,
    guest_anchor: int, keep_guest_anchor: int, *, name: str,
) -> tuple[Structure, int]:
    """Attach a fragment and preserve one guest port for a subsequent step."""
    host_port = _port_at_atom(host, host_anchor)
    guest_port = _port_at_atom(guest, guest_anchor)
    keep_port = _port_at_atom(guest, keep_guest_anchor)
    expose: dict[Port, tuple[int, Port]] = {}
    output_index = 1
    for port in host.ports:
        if port != host_port:
            expose[single(output_index)] = (0, port)
            output_index += 1
    expose[single(output_index)] = (1, keep_port)
    new_anchor = len(host.atoms) + keep_guest_anchor
    result = substitute(host, {host_port: (guest, guest_port)}, expose=expose, name=name)
    return result, new_anchor



def _bis_imine(amine: str) -> Structure:
    heteroaromatic = pyridine if amine == "2-aminopyridine" else pyrimidine
    result = _compose_binaphthyl_base()
    for position in (2, 13):
        result, imine_n = _attach_preserving_guest_port(
            result, position, imine, guest_anchor=0, keep_guest_anchor=1,
            name=f"binaphthyl bis({amine} imine)",
        )
        result = _replace_with_fragment(
            result, imine_n, heteroaromatic, name=f"binaphthyl bis({amine} imine)"
        )
    result.name = f"binaphthyl bis({amine} imine)"
    return result


def _diformyl() -> Structure:
    result = _compose_binaphthyl_base()
    result = _replace_with_fragment(result, 2, carbonyl, name="2-formyl-1,1'-binaphthyl")
    return _replace_with_fragment(result, 13, carbonyl, name="2,2'-diformyl-1,1'-binaphthyl")


binol_2_2_diformyl = _diformyl()
binol_bis_2_pyridyl_imine = _bis_imine("2-aminopyridine")
binol_bis_2_pyrimidyl_imine = _bis_imine("2-aminopyrimidine")
