# Agent instructions

## Molecular structure authoring

- Prefer a **fragment-based approach** when writing molecular structures and examples. Compose structures from reusable fragments in `molcule_desc.library` using operations such as `bond()`, `ring()`, and `substitute()`, rather than manually constructing atom and bond graphs or defining one-off atoms when an existing library fragment is suitable.
- Build chains with `bond()` and close cycles with `ring()` where practical. Explicitly expose only the ports that need to remain available for subsequent composition; unexposed ports are capped by the composition operation's default behavior.
- Reuse library fragments for common atoms and groups (for example, `methylene`, `hydroxyl`, and `carbonyl`) instead of defining duplicate local fragments.

## Library fragment ports

- **Expose all chemically valid attachment ports by default** for fragments added to `molcule_desc.library`. This keeps fragments reusable in different molecular contexts and lets callers decide which ports to consume or expose during composition.
- Make an exception when the attachment point is part of the fragment's meaning and changing it would change the intended structure or make the fragment ambiguous. For example, a named substituent such as `tert_butyl` should expose the chemically intended attachment point(s), rather than indiscriminately exposing every atom valence as an interchangeable connection.
- Document any intentional restriction on exposed ports so users can understand why the fragment is not fully open for substitution.
