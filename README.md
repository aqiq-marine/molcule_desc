# molcule_desc

**molcule_desc** is an experimental, fragment-based language for describing molecular structures compositionally.

Instead of writing every atom and bond directly, a user defines reusable structures and composes them by connecting their terminal bonds. A structure may be a small substituent, a large fragment, or a complete molecule: the language does not distinguish these as separate kinds of objects.

The project prioritizes a small core language, explicit composition rules, and predictable behavior. Human-friendly names and editor autocomplete can make common attachment positions convenient without requiring the language itself to infer chemical nomenclature.

> **Status:** Early prototype. The Python API in `molcule_desc/` implements an atom/bond graph, runtime port checks, explicit port exports, and terminal-port completion. The surface language remains a design proposal.

## Design principles

- **Composition over atom-by-atom notation.** Build larger structures from reusable structures.
- **One concept for fragments and molecules.** A structure is a structure, whether or not it has unconsumed terminal bonds.
- **Terminal bonds are attachment points.** Terminal single, double, and triple bonds are the primitive connection points.
- **Connections are explicit.** Definitions list how source ports map to the ports exposed by the resulting structure.
- **No inferred port mapping.** The language implementation must not guess which attachment points a composition uses or exports.
- **Implicit hydrogen completion.** Unused terminal valences are completed with hydrogen according to the rules defined in the specification.
- **Friendly surface syntax.** Names such as positional labels may be provided by libraries and editor tooling; internally, they resolve to explicit port identifiers.

## Conceptual examples

The following examples illustrate the intended model. They are **pseudocode**, not a promise of finalized syntax.

A structure can be composed from other structures:

```text
ethylene(phenyl, phenyl)
```

A substituent can be attached to a named position:

```text
benzene(4:OH)
```

A larger structure can be defined by composing structures and explicitly listing the ports that remain available:

```text
structure example = compose(A, B) {
    connect A.single[1] with B.single[2]

    export:
        single[1] = A.single[2]
        single[2] = B.single[1]
}
```

The exact grammar is still under design. See [SPECIFICATION.md](SPECIFICATION.md) for the current design proposal and unresolved questions.

## Intended workflow

1. Define or import reusable structures.
2. Compose them by supplying structures to attachment positions.
3. Explicitly connect terminal bonds when a composition requires a ring closure or another intramolecular connection.
4. Explicitly declare which unused terminal bonds are exposed by the resulting structure.
5. Resolve any remaining terminal valences to implicit hydrogens when producing a completed molecular structure.

## Scope

The language describes **molecular structure**, not synthetic procedures or reaction mechanisms. Operations such as protecting groups, reaction transformations, iteration, and general-purpose control flow are not part of the proposed core.

## Python prototype

The core API is implemented with runtime validation. Bond order is the only enum; ports use a bond kind and a structure-local integer identifier.

```python
from molcule_desc import (
    acetone, ethylene, methanimine,
    hydrogen_single, carbon_double_two_single, carbon_triple_single, oxygen_double,
)

print(ethylene.finalize().formula_string())    # C2H4
print(methanimine.finalize().formula_string())  # CH3N
print(acetone.finalize().formula_string())      # C3H6O

# atom() declares an atom and its valence; ported_atom() gives that atom ports.
from molcule_desc import atom, ported_atom, single, double

carbon = atom("C", 4)
carbon_template = ported_atom(
    carbon, [double(1), single(1), single(2)], name="carbon_double_two_single"
)

# Atom-topology templates, named by their available bond ports.
print(carbon_double_two_single.ports)
print(carbon_triple_single.ports)
print(oxygen_double.ports)
print(hydrogen_single.ports)
```

`atom(element, valence)` creates an atom specification with an explicit valence (for example, `atom("C", 4)` or `atom("O", 2)`). `ported_atom(atom, ports)` turns that atom into a one-atom `Structure` with explicit single, double, or triple bond ports. This separates the atom's element/valence definition from the attachment template.

Structures are composed with `bond(left, left_port, right, right_port, expose=...)` for connecting two structures, and `ring(structure, first_port, second_port, expose=...)` for connecting two ports on one structure.
The `expose` mapping is explicit: each result-local port points to a source port through `PortRef(0, ...)` (left input) or `PortRef(1, ...)` (right input). Compatible ports are consumed and a graph bond is added. Unexported, unconsumed ports are terminalized during composition; call `finalize()` to resolve the ports that remain exposed by the final structure. To connect two ports in one structure, use `ring(structure, port1, port2, expose=...)`. Its `expose` mapping maps each result-local port directly to an existing port on that structure, so `PortRef` is not needed.

This is a small prototype, not a full cheminformatics toolkit. It currently supports common valences and does not implement aromaticity or stereochemistry in its own graph model. Optional RDKit integration converts structures to SMILES:

```python
from molcule_desc import acetone, to_smiles

print(to_smiles(acetone))  # CC(C)=O
```

Install RDKit integration with `pip install -e '.[chemistry]'`, then run tests with `pip install -e . pytest` followed by `pytest`. The converter finalizes exposed ports by default, sanitizes the RDKit molecule, and emits canonical SMILES by default. Aromaticity and stereochemistry are not inferred by the converter.

## Built-in structure library

The prototype includes reusable hydrocarbon skeletons, functional groups, heterocycles,
halogen substituents, sulfur/phosphorus groups, and linkers. Library structures expose
their remaining valences as single-bond ports; when a port is left unused, `finalize()`
completes it with hydrogen.

```python
from molcule_desc import benzene, hydroxyl, to_smiles

# Every benzene position is an available attachment port.
assert len(benzene.ports) == 6
print(to_smiles(benzene))  # c1ccccc1
```

Available names include:

- **Hydrocarbons:** `methane`, `ethane`, `propane`, `ethylene`, `acetylene`,
  `cyclohexane`, `benzene`, `naphthalene`
- **Oxygen groups:** `hydroxyl`, `ether`, `carbonyl`, `carboxyl`, `ester`
- **Nitrogen groups:** `amino`, `amine`, `amide`, `nitrile`, `imine`
- **Heterocycles:** `pyridine`, `pyrrole`, `thiophene`, `imidazole`
- **Halogens:** `fluoro`, `chloro`, `bromo`, `iodo`
- **Sulfur/phosphorus:** `thiol`, `thioether`, `sulfone`, `phosphine`,
  `phosphine_oxide`, `phosphate_ester`
- **Linkers and examples:** `methylene`, `ethylene_linker`, `phenyl`,
  `phenylene`, `biphenyl`, `diphenyl_ether`, `ethanol`

These are graph templates, not a chemically exhaustive database. In particular,
aromatic rings are stored as Kekulé structures and aromaticity is perceived by RDKit
during SMILES conversion. Port positions currently have numeric identifiers; named
positions such as ortho/meta/para are not yet modeled.

## Documentation

- [Language specification (proposal)](SPECIFICATION.md)

## Name

The repository name is `molcule_desc`. The project explores a compositional description language for molecules.

## Contributing

Design feedback and examples of structures that are difficult to express predictably are welcome. Since the language is experimental, proposals that simplify the core semantics are especially useful.
