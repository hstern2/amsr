from rdkit import Chem
from rdkit.Chem import AllChem

import amsr
from amsr.tokens import ToTokens, _insert_implicit_carbon, _remove_implicit_carbon

caffeine_smi = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
taxol_smi = (
    "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]"
    "([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)N"
    "C(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C"
)


def test_version() -> None:
    assert amsr.__version__


def test_caffeine() -> None:
    assert amsr.CheckSmiles(caffeine_smi)


def test_taxol() -> None:
    assert amsr.CheckSmiles(taxol_smi)


def test_cage() -> None:
    assert amsr.CheckAMSR("CCccCccc6oC..CCCC6C6.6")


def test_no_stereo() -> None:
    assert amsr.CheckAMSR(amsr.FromSmiles(taxol_smi, useStereo=False))


def test_canonical() -> None:
    m = Chem.MolFromSmiles(caffeine_smi)
    s = amsr.FromSmiles(Chem.MolToSmiles(m, doRandom=True), canonical=True)
    for _ in range(20):
        assert amsr.FromSmiles(Chem.MolToSmiles(m, doRandom=True), canonical=True) == s


def _mol_with_3d(smi):
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    return Chem.RemoveHs(mol)


def test_tokenize_bracketed_group_with_digits() -> None:
    # bracketed groups may contain digits after letters (e.g., [NMe2], [NO2], [CF3])
    assert ToTokens("[NMe2]") == ["[NMe2]"]
    assert ToTokens("[NO2]") == ["[NO2]"]
    assert ToTokens("[CF3]") == ["[CF3]"]
    # isotope-style brackets (digits first) must still tokenize as one atom
    assert ToTokens("[12C]") == ["[12C]"]
    assert ToTokens("[13CH3]") == ["[13CH3]"]


def test_implicit_carbon_insert() -> None:
    assert _insert_implicit_carbon("^^__") == "^^C__"
    assert _insert_implicit_carbon("^^__>>") == "^^C__C>>"
    assert _insert_implicit_carbon("C^^__C") == "C^^C__C"
    assert _insert_implicit_carbon("C^^C__C") == "C^^C__C"
    # multi-char dihedrals
    assert _insert_implicit_carbon("^\\<\\") == "^\\C<\\"
    assert _insert_implicit_carbon("_/<\\") == "_/C<\\"


def test_implicit_carbon_remove() -> None:
    assert _remove_implicit_carbon(["C", "^^", "C", "__", "C"]) == ["C", "^^", "__", "C"]
    assert _remove_implicit_carbon(["C", "^^", "C", "__", "C", "^^", "C"]) == [
        "C",
        "^^",
        "__",
        "^^",
        "C",
    ]
    # modified carbons should NOT be removed
    assert _remove_implicit_carbon(["C", "^^", "C'", "__", "C"]) == [
        "C",
        "^^",
        "C'",
        "__",
        "C",
    ]
    assert _remove_implicit_carbon(["C", "^^", "c", "__", "C"]) == [
        "C",
        "^^",
        "c",
        "__",
        "C",
    ]


def test_implicit_carbon_roundtrip() -> None:
    for smi in ["CCCCC", "CCCCCC", "CCCCCCC", "CCCCCCCC"]:
        mol = _mol_with_3d(smi)
        s = amsr.FromMol(mol)
        # roundtrip check
        mol2 = amsr.ToMol(s)
        assert Chem.MolToInchi(mol, options="-FixedH") == Chem.MolToInchi(
            mol2, options="-FixedH"
        ), f"InChI mismatch for {smi}: {s}"
        # verify dihedral info survives roundtrip
        dihedral: dict[tuple[int, int, int, int], int] = {}
        amsr.ToMol(s, dihedral=dihedral)
        assert len(dihedral) > 0, f"No dihedrals decoded for {smi}: {s}"
