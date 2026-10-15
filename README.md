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
    XS, XD, XT,
)

print(ethylene.finalize().formula_string())    # C2H4
print(methanimine.finalize().formula_string())  # CH3N
print(acetone.finalize().formula_string())      # C3H6O

# atom() declares an atom and its valence; ported_atom() gives that atom ports.
from molcule_desc import atom, ported_atom, single, double

carbon = atom("C", 4)
carbon_template = ported_atom(
    carbon, [double(1), single(1), single(2)], name="XD"
)

# Terminal templates:
# XS = H with one single-bond port
# XD = C with one double-bond port and two single-bond (H) ports
# XT = C with one triple-bond port and one single-bond (H) port
print(XD.ports)
```

`atom(element, valence)` creates an atom specification with an explicit valence (for example, `atom("C", 4)` or `atom("O", 2)`). `ported_atom(atom, ports)` turns that atom into a one-atom `Structure` with explicit single, double, or triple bond ports. This separates the atom's element/valence definition from the attachment template.

Structures are composed with `bond(left, left_port, right, right_port, expose=...)`.
The `expose` mapping is explicit: each result-local port points to a source port through `PortRef(0, ...)` (left input) or `PortRef(1, ...)` (right input). Compatible ports are consumed and a graph bond is added. Unexported, unconsumed ports are terminalized during composition; call `finalize()` to resolve the ports that remain exposed by the final structure. To connect two ports in one structure, pass the same structure twice and set `intramolecular=True`.

This is a small prototype, not a full cheminformatics toolkit. It currently supports common neutral valences and does not implement aromaticity, stereochemistry, formal-charge-aware valence rules, or SMILES conversion. Install and test with `pip install -e . pytest` followed by `pytest`.

## Documentation

- [Language specification (proposal)](SPECIFICATION.md)

## Name

The repository name is `molcule_desc`. The project explores a compositional description language for molecules.

## Contributing

Design feedback and examples of structures that are difficult to express predictably are welcome. Since the language is experimental, proposals that simplify the core semantics are especially useful.
