import pytest

from molcule_desc import (
    benzene,
    fuse,
    single,
    tetramethylene,
    to_smiles,
)


def test_fuse_creates_tetralin_skeleton_and_preserves_other_ports():
    fused = fuse(
        benzene,
        tetramethylene,
        cyclic_ports=(single(1), single(2)),
        chain_ports=(single(1), single(10)),
        name="tetralin",
    )

    # Benzene has six exposed ports and tetramethylene has ten; two ports
    # from each fragment are consumed by the two new fusion bonds.
    assert len(fused.atoms) == 10
    assert len(fused.bonds) == 11
    assert len(fused.ports) == 12

    finalized = fused.finalize()
    assert finalized.formula_string() == "C10H12"
    assert finalized.ports == {}
    assert to_smiles(fused)


def test_fuse_rejects_missing_ports():
    with pytest.raises(ValueError, match="no exposed port"):
        fuse(
            benzene,
            tetramethylene,
            cyclic_ports=(single(1), single(99)),
            chain_ports=(single(1), single(10)),
        )


def test_fuse_rejects_repeated_ports():
    with pytest.raises(ValueError, match="two distinct ports"):
        fuse(
            benzene,
            tetramethylene,
            cyclic_ports=(single(1), single(1)),
            chain_ports=(single(1), single(10)),
        )
