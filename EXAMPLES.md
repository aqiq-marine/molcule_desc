# Language Examples

This document illustrates the intended feel of the language. **The syntax is provisional**: these examples describe the design direction, not a finalized grammar or executable implementation.

The language treats every structure uniformly. A structure may be a small functional group, a fragment, or a complete molecule; the distinction is not part of the language's type system. New structures are composed from existing structures rather than created by arbitrary atom-by-atom editing.

## Conventions used below

- Names such as `ethane`, `benzene`, and `OH` refer to structures supplied by a library.
- A structure exposes named ports. Port names are local to that structure and need not follow IUPAC numbering.
- For example, the six replaceable hydrogens of `ethane` are exposed as ports `1` through `6`.
- `4` and `9` are just port names in the biphenyl example below; the language does not give them special meaning.
- A port includes a consumable port atom. When left unconnected in a finalized structure, it is evaluated using its terminal form: a single-bond port as hydrogen, a double-bond port as one side of ethene, and a triple-bond port as one side of ethyne.
- When compatible ports are connected, the corresponding port atoms are consumed according to the connection rule. The resulting structure explicitly declares which ports remain exposed.
- The examples use concise illustrative notation. Exact syntax, port compatibility, and validation rules remain to be specified.

## 1. Substitute one hydrogen of ethane

```text
nitroethane = ethane(1: nitro)
```

Conceptually, this starts with the library structure `ethane` and connects its port `1` to the single-bond port of `nitro`. The connected ports are consumed. The unused ethane ports are finalized as hydrogens.

The six ports of `ethane` are explicit; they are not inferred from the name `ethane` during composition.

## 2. Define a reusable structure from other structures

```text
phenol = benzene(4: OH)
```

This defines `phenol` by composing the existing `benzene` and `OH` structures. Port `4` is simply the selected port on `benzene`; it is not a reserved word or a special syntax for aromatic substitution.

The definition of `benzene` in the library determines its available ports and their internal connections.

## 3. Use arbitrary port names

```text
dihydroxybiphenyl = biphenyl(
    4: OH,
    9: OH
)
```

Here, `4` and `9` are names chosen for two ports exposed by `biphenyl`. They could instead be called `left_para` and `right_para`, provided the structure definition exposes those names.

The language itself does not need to know what “4,4′” means. A library or editor may offer familiar locants as convenient port names.

## 4. Compose structures recursively

```text
hydroxyphenyl = benzene(4: OH)
linked_structure = biphenyl(left: hydroxyphenyl)
```

This illustrates the compositional idea: a structure defined by the user can be used just like a library structure in another definition.

This example is schematic. A real `biphenyl` definition must specify its connection ports and how the ports of `hydroxyphenyl` map into the resulting structure. No port mapping should be guessed from atom order or structure names.

## 5. Declare which ports survive composition

A structure definition must explicitly describe the ports exposed by its result. For example, suppose `A` and `B` each expose two single-bond ports:

```text
AB = compose(A, B) {
    connect A.1 with B.2

    export:
        left  = A.2
        right = B.1
}
```

The notation means:

1. Connect `A` port `1` to `B` port `2`.
2. Consume those two connected ports.
3. Expose the remaining `A` port `2` as `AB.left`.
4. Expose the remaining `B` port `1` as `AB.right`.

The names `left` and `right` are newly assigned names. Port propagation is explicit, not based on ordering.

## 6. Close a chain into a ring

A chain structure can expose ports at both ends. Those two ports can then be connected to each other within the same composition:

```text
cyclohexane = compose(hexane) {
    connect hexane.left with hexane.right
}
```

This is an illustration of intramolecular port connection, not a special `ring()` operation. The library definition of `hexane` must expose suitable end ports, and the connection must consume them according to the same rules used for connections between separate structures.

## 7. Double- and triple-bond ports

The same composition model can use ports of different bond orders:

```text
alkene_based_structure = compose(A, B) {
    connect A.double_port with B.double_port
}

alkyne_based_structure = compose(C, D) {
    connect C.triple_port with D.triple_port
}
```

These names are illustrative. The library definitions determine which ports are double- or triple-bond ports. A finalized, unused double-bond port is represented by the appropriate half-ethene terminal structure; a finalized, unused triple-bond port is represented by the appropriate half-ethyne terminal structure.

The exact graph rewrite for consuming these port atoms must be defined precisely so that the resulting bond orders and valences are correct.

## Design principle

The core language should only need a small set of operations:

- refer to a named structure;
- define a structure by composing existing structures;
- connect compatible ports;
- explicitly export and rename remaining ports;
- finalize a structure by evaluating its unconnected ports.

Library structures and user-defined structures follow the same rules. The library supplies useful starting points; it does not define a separate kind of object.
