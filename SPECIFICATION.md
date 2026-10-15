# molcule_desc Language Specification

**Status:** Design proposal  
**Version:** 0.1-draft

This document records the current design direction for molcule_desc. It is intentionally explicit about design principles while leaving syntax details provisional where they have not yet been settled.

## 1. Purpose

molcule_desc is a compositional language for describing molecular structures. It aims to make structures readable and reusable by treating a named structure as a callable building block whose attachment positions can receive other structures.

The language should have a small core and deterministic composition semantics. It should not require the parser to infer chemical intent from a structure's name or from an IUPAC numbering convention.

## 2. Core principles

1. **A single structure concept.** A substituent, fragment, intermediate structure, and complete molecule are all represented by the same kind of value. There is no separate language-level `Fragment` versus `Molecule` type.
2. **Composition is the default construction mechanism.** Larger structures are defined by composing existing structures rather than by manually specifying every atom and bond.
3. **Terminal bonds are attachment points.** The core recognizes terminal single, double, and triple bonds as connectable ports.
4. **Port mapping is explicit.** Every composition definition declares which source ports are consumed and which source ports are exposed by the result.
5. **No implicit port guessing.** The implementation must not infer a connection from atom order, structural similarity, or nomenclature unless an explicit user-facing alias has already resolved to a specific port.
6. **Hydrogen completion is implicit.** Unused terminal valences are completed according to the completion rules in Section 6.
7. **Surface syntax may be ergonomic.** Human-readable position names and autocomplete are tooling/library features that resolve to explicit ports. They do not replace explicit semantics.

## 3. Structures and definitions

A structure represents a molecular graph together with a set of exposed terminal-bond ports.

A structure definition may:
- refer to existing structures;
- compose structures by supplying them to attachment positions;
- connect compatible terminal-bond ports, including ports belonging to the same composed structure;
- explicitly declare the ports exposed by the resulting structure.

Definitions may be nested: a structure created from other structures can itself be reused as a structure in another definition.

The exact declaration grammar is not finalized. The following is illustrative pseudocode only:

```text
structure diphenylethane = compose(ethylene, phenyl, phenyl) {
    # The actual syntax for composition is provisional.
    # The resulting structure explicitly exports selected ports.
}
```

The implementation must preserve a deterministic mapping from each exposed result port to its source port or source structure location.

## 4. Terminal-bond ports

### 4.1 Port kinds

The initial core defines three port kinds:

- `single`: terminal single bond
- `double`: terminal double bond
- `triple`: terminal triple bond

A port is identified by its kind and an identifier within the structure, for example `single[1]` or `double[2]`. Identifiers are local to a structure definition; they are not inherently atom indices or IUPAC locants.

A port represents an exposed terminal bond, not a free-standing atom and not a general-purpose arbitrary bond operation.

### 4.2 Connection compatibility

In the initial model, ports may be connected only when they have the same kind:

- single + single
- double + double
- triple + triple

Connecting two ports consumes both ports and creates the corresponding bond between their attached atoms. The consumed ports are no longer available for further composition.

This is a syntactic compatibility rule, not a complete guarantee of chemical validity. The resulting structure must also pass valence and graph validation. More detailed treatment of charge, aromaticity, unusual valence states, and other chemical edge cases remains future work.

### 4.3 Explicit intramolecular connections

The same connection mechanism is used when both ports belong to the same composed structure. This supports ring closure and other intramolecular links without requiring a distinct core primitive. The Python prototype also provides `ring(...)` as a convenience wrapper around the same intramolecular connection operation; it does not introduce different semantics.

For example, a cyclic structure can be defined by constructing an open chain and explicitly connecting its two appropriate terminal ports. The particular port identifiers must be listed by the definition.

## 5. Composition and port export

When structures are composed:

1. Each input structure contributes its atoms, bonds, and exposed ports.
2. Every declared connection consumes exactly the two referenced ports.
3. Any port not consumed by a connection remains available internally.
4. The definition explicitly declares which remaining ports are exposed by the resulting structure, and assigns result-local identifiers to them.
5. Unexported ports are not available to callers. They are internal to the definition and are handled by the definition's completion policy.

The exact syntax for composition and export declarations is provisional. The required semantic property is not: callers must be able to determine the resulting structure's available ports from the definition alone.

Example pseudocode:

```text
structure linked = compose(A, B) {
    connect A.single[1] with B.single[2]

    export:
        single[1] = A.single[2]
        single[2] = B.single[1]
}
```

Here, `A.single[1]` and `B.single[2]` are consumed. The two listed ports are exposed by `linked`. This example demonstrates semantics, not finalized syntax.

A definition must not silently renumber or export ports based on implementation-dependent atom traversal order.

## 6. Implicit hydrogen completion

The intended user-facing rule is that terminal bonds left unused after composition are completed with hydrogen when a completed molecular structure is requested.

The initial design targets common neutral organic structures. Informally, a carbon with:
- one remaining single-bond valence is completed as `CH3` where appropriate;
- two remaining single-bond valences is completed as `CH2`;
- three remaining single-bond valences is completed as `CH`;
- four single-bond valences is completed as `C`.

These examples are explanatory rather than a complete valence algorithm. Actual completion must be based on the atom's existing bond orders, element, formal charge, and supported valence rules—not solely on the number of ports.

The following behavior is required:
- Completion must be deterministic.
- Explicitly connected ports are not hydrogen-capped.
- A completed structure must not retain unresolved terminal ports.
- If completion is ambiguous or violates the supported valence model, the implementation must report an error rather than silently choose a structure.

The exact point at which completion occurs—at definition time, composition time, or finalization time—remains to be specified. The preferred behavior is to preserve exposed ports while composing and perform hydrogen completion only when a structure is finalized.

## 7. Variables and reusable definitions

The language is intended to support named definitions and variables. A defined structure can be used as an argument in another definition, enabling arbitrary levels of composition.

The precise declaration syntax, parameter binding rules, and name-resolution rules are not finalized. These features must not alter the port semantics: passing a structure as an argument must preserve its graph and its exposed port identities.

## 8. Human-friendly position labels

The language may offer labels such as `4`, `4'`, `ortho`, `meta`, or `para` when a library defines those labels for a particular structure.

These labels are aliases, not universal rules of the core language. For example, a library may define a familiar position label for a benzene-derived structure and map it to a particular port identifier. The compiler or editor resolves the label to that explicit port before composition.

Requirements:
- Every alias must resolve unambiguously within its structure/library context.
- Aliases must not be assumed to follow one universal numbering scheme.
- Tooling should provide autocomplete and display the structural context of candidate ports.
- The internal representation must retain explicit port identities after alias resolution.

## 9. Errors and validation

At minimum, an implementation should diagnose:
- references to nonexistent ports;
- attempts to consume the same port more than once;
- connections between incompatible port kinds;
- exported ports that were consumed or do not exist;
- duplicate result-port identifiers;
- structures that cannot be completed under the supported valence rules;
- definitions whose port mappings cannot be resolved deterministically.

Validation of chemical plausibility is separate from parsing. A syntactically valid definition is not necessarily a chemically valid molecule.

## 10. Out of scope for the core

The current design does not require dedicated primitives for:
- arbitrary `bond()` construction;
- ring notation or ring-specific operations;
- protecting groups or reaction transformations;
- loops, `map`, or general-purpose iteration;
- general reaction mechanisms or synthetic procedures.

Rings are expressed using explicit connections between terminal ports. A host-language convenience helper such as Python's `ring(...)` may wrap this operation, but ring closure is not a separate language-level primitive. Common molecular building blocks belong in libraries rather than being hard-coded into the language.

Stereochemistry, isotopes, formal charges, aromaticity, coordination compounds, polymers, and nonstandard valence states are not yet fully specified.

## 11. Open design questions

The following decisions remain open and should be settled with examples and tests before the syntax is treated as stable:

1. **Surface grammar:** Exact syntax for structure definitions, composition, connection declarations, exports, and named arguments.
2. **Port identity:** Whether identifiers are declared directly in definitions, generated and then named, or both.
3. **Unexported ports:** Whether these are always hydrogen-completed at the end of a definition, or can remain private for later internal use.
4. **Hydrogen completion:** Supported elements, valence tables, formal charges, aromatic atoms, and error conditions.
5. **Aromaticity:** How aromatic systems are represented internally and how their ports interact with single/double/triple terminal-bond rules.
6. **Stereochemistry and isotopes:** Syntax and preservation through composition.
7. **Validation backend:** The boundary between language-level checks and validation delegated to a cheminformatics toolkit.
8. **Serialization:** Which formats (for example SMILES, molfile, or a molecular graph representation) are supported as input/output targets.

## 12. Compatibility principle

Until a versioned grammar and semantic test suite are published, examples in this document are illustrative. A future implementation may change surface syntax, but should preserve the core guarantees: composition of structures, explicit port connections, explicit port exports, deterministic port mapping, and predictable hydrogen completion.
