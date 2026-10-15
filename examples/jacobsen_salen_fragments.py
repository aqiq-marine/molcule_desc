"""Fragment-composed example: Jacobsen-type salen ligand.

The molecule is assembled from reusable library fragments. This example
intentionally avoids manually enumerating atoms and bonds. Stereochemistry is
not represented by the current Structure model.
"""

from __future__ import annotations

from dataclasses import dataclass

from molcule_desc import (
    Port, PortRef, Structure, bond, methane, benzene, hydroxyl, methanimine,
    cyclohexane, single, to_smiles,
)


@dataclass
class Fragment:
    """A structure plus stable symbolic names for its exposed attachment ports."""

    structure: Structure
    ports: dict[str, Port]


def attach(
    base: Fragment,
    base_port_name: str,
    substituent: Fragment,
    substituent_port_name: str,
    *,
    name: str,
) -> Fragment:
    """Attach two fragments and preserve all unconsumed ports by symbolic name."""
    bp = base.ports[base_port_name]
    sp = substituent.ports[substituent_port_name]
    remaining: list[tuple[str, int, Port]] = []
    for label, port in base.ports.items():
        if port != bp:
            remaining.append((label, 0, port))
    for label, port in substituent.ports.items():
        if port != sp:
            remaining.append((label, 1, port))

    exports: dict[Port, PortRef] = {}
    new_names: dict[str, Port] = {}
    for index, (label, source, old_port) in enumerate(remaining, start=1):
        new_port = Port(old_port.bond, index)
        exports[new_port] = PortRef(source, old_port)
        new_names[label] = new_port

    combined = bond(
        base.structure,
        bp,
        substituent.structure,
        sp,
        expose=exports,
        name=name,
    )
    return Fragment(combined, new_names)


def carbon_fragment(name: str = "carbon") -> Fragment:
    """A carbon atom with four named single-bond attachment ports."""
    from molcule_desc import atom, ported_atom

    structure = ported_atom(
        atom("C", 4),
        [single(i) for i in range(1, 5)],
        name=name,
    )
    return Fragment(structure, {f"p{i}": single(i) for i in range(1, 5)})


def tert_butyl() -> Fragment:
    """Build tert-butyl from one central carbon and three methyl carbon fragments."""
    result = carbon_fragment("tert-butyl central carbon")
    # Three carbon fragments each consume one port of the central carbon;
    # the remaining ports become the three methyl groups' C-H valences.
    for i in range(3):
        methyl = carbon_fragment(f"tert-butyl methyl {i + 1}")
        result = attach(result, f"p{i + 1}", methyl, "p1", name="tert-butyl")
    result.ports = {"attachment": result.ports["p4"]}
    return result


def salicylidene_half() -> Fragment:
    """Build one 3,5-di-tert-butylsalicylidene half from ring and group fragments."""
    result = Fragment(
        benzene,
        {f"aryl{i}": single(i) for i in range(1, 7)},
    )
    # Position aliases are local to this chosen benzene port convention.
    result = attach(result, "aryl3", tert_butyl(), "attachment", name="aryl tert-butyl")
    # Resolve aliases from the returned mapping; labels for unconsumed ports survive.
    result = attach(result, "aryl5", tert_butyl(), "attachment", name="di-tert-butyl aryl")
    result = attach(result, "aryl2", Fragment(hydroxyl, {"attachment": single(1), "hydrogen": single(2)}),
                    "attachment", name="salicylidene phenol")
    result = attach(result, "aryl1", Fragment(methanimine, {
        "imine_carbon_a": single(1),
        "imine_carbon_b": single(2),
        "nitrogen": single(3),
    }), "imine_carbon_a", name="salicylidene imine")
    # Expose only the imine nitrogen for connection to the diamine fragment.
    nitrogen_port = result.ports["nitrogen"]
    return Fragment(result.structure, {"nitrogen": nitrogen_port})


def build_jacobsen_salen() -> Structure:
    """Compose two salicylidene halves with a 1,2-disubstituted cyclohexane."""
    left = salicylidene_half()
    right = salicylidene_half()
    diamine = Fragment(
        cyclohexane,
        {
            # In the built-in cyclohexane template, ports 1 and 3 are on
            # adjacent ring carbons (each carbon has two ports).
            "amine1": single(1),
            "amine2": single(3),
        },
    )
    with_left = attach(diamine, "amine1", left, "nitrogen", name="left imine attached")
    # The cyclohexane's second port alias survives the first attachment.
    return attach(with_left, "amine2", right, "nitrogen", name="Jacobsen salen").structure


if __name__ == "__main__":
    ligand = build_jacobsen_salen()
    print("Formula:", ligand.finalize().formula_string())
    print("SMILES:", to_smiles(ligand))
