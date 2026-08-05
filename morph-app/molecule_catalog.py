"""Curated molecule catalog for quick local lookup.

SMILES values in this file were resolved from PubChem PUG REST and kept with
their PubChem CIDs for traceability.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class MoleculeReference:
    name: str
    smiles: str
    collection: str
    aliases: tuple[str, ...] = ()
    pubchem_cid: Optional[int] = None


def molecule(
    name: str,
    smiles: str,
    *,
    cid: Optional[int],
    collection: str,
    aliases: tuple[str, ...] = (),
) -> MoleculeReference:
    return MoleculeReference(
        name=name,
        smiles=smiles,
        collection=collection,
        aliases=aliases,
        pubchem_cid=cid,
    )


DEFAULT_SMILES_1 = (
    "CC[C@H]1C[C@@H]2C[C@@H]3[C@H]1N(C2)CCC4=C3NC5=C4C=C(C=C5)OC"  # Ibogaine, PubChem CID 197060
)
DEFAULT_SMILES_2 = "C1C[C@@H]2[C@H](C[C@H]1N2)C3=CN=C(C=C3)Cl"  # Epibatidine, PubChem CID 854023


KNOWN_MOLECULES: dict[str, MoleculeReference] = {
    "ibogaine": molecule(
        "Ibogaine",
        DEFAULT_SMILES_1,
        cid=197060,
        collection="Natural product",
    ),
    "aspirin": molecule(
        "Aspirin",
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        cid=2244,
        collection="FDA-approved drug",
    ),
    "caffeine": molecule(
        "Caffeine",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cid=2519,
        collection="Natural product",
    ),
    "acetaminophen": molecule(
        "Acetaminophen",
        "CC(=O)NC1=CC=C(C=C1)O",
        cid=1983,
        collection="FDA-approved drug",
        aliases=("paracetamol",),
    ),
    "ibuprofen": molecule(
        "Ibuprofen",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        cid=3672,
        collection="FDA-approved drug",
    ),
    "nicotine": molecule(
        "Nicotine",
        "CN1CCC[C@H]1C2=CN=CC=C2",
        cid=89594,
        collection="Natural product",
    ),
    "dopamine": molecule(
        "Dopamine",
        "C1=CC(=C(C=C1CCN)O)O",
        cid=681,
        collection="FDA-approved drug",
    ),
    "serotonin": molecule(
        "Serotonin",
        "C1=CC2=C(C=C1O)C(=CN2)CCN",
        cid=5202,
        collection="Natural product",
    ),
    "aconitine": molecule(
        "Aconitine",
        "CCN1C[C@@]2([C@@H](C[C@@H]([C@@]34[C@@H]2[C@H]([C@@H](C31)[C@@]5([C@@H]6[C@H]4C[C@@]([C@@H]6OC(=O)C7=CC=CC=C7)([C@H]([C@@H]5O)OC)O)OC(=O)C)OC)OC)O)COC",
        cid=245005,
        collection="Natural product",
    ),
    "artemisinin": molecule(
        "Artemisinin",
        "C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C",
        cid=68827,
        collection="Natural product",
    ),
    "atropine": molecule(
        "Atropine",
        "CN1[C@@H]2CC[C@H]1CC(C2)OC(=O)C(CO)C3=CC=CC=C3",
        cid=174174,
        collection="Natural product",
    ),
    "berberine": molecule(
        "Berberine",
        "COC1=C(C2=C[N+]3=C(C=C2C=C1)C4=CC5=C(C=C4CC3)OCO5)OC",
        cid=2353,
        collection="Natural product",
    ),
    "capsaicin": molecule(
        "Capsaicin",
        "CC(C)/C=C/CCCCC(=O)NCC1=CC(=C(C=C1)O)OC",
        cid=1548943,
        collection="Natural product",
    ),
    "camphor": molecule(
        "Camphor", "CC1(C2CCC1(C(=O)C2)C)C", cid=2537, collection="Natural product"
    ),
    "cannabidiol": molecule(
        "Cannabidiol",
        "CCCCCC1=CC(=C(C(=C1)O)[C@@H]2C=C(CC[C@H]2C(=C)C)C)O",
        cid=644019,
        collection="Natural product",
    ),
    "cannabinol": molecule(
        "Cannabinol",
        "CCCCCC1=CC(=C2C(=C1)OC(C3=C2C=C(C=C3)C)(C)C)O",
        cid=2543,
        collection="Natural product",
    ),
    "colchicine": molecule(
        "Colchicine",
        "CC(=O)N[C@H]1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC",
        cid=6167,
        collection="Natural product",
    ),
    "curcumin": molecule(
        "Curcumin",
        "COC1=C(C=CC(=C1)/C=C/C(=O)CC(=O)/C=C/C2=CC(=C(C=C2)O)OC)O",
        cid=969516,
        collection="Natural product",
    ),
    "digoxin": molecule(
        "Digoxin",
        "C[C@@H]1[C@H]([C@H](C[C@@H](O1)O[C@@H]2[C@H](O[C@H](C[C@@H]2O)O[C@@H]3[C@H](O[C@H](C[C@@H]3O)O[C@H]4CC[C@]5([C@@H](C4)CC[C@@H]6[C@@H]5C[C@H]([C@]7([C@@]6(CC[C@@H]7C8=CC(=O)OC8)O)C)O)C)C)C)O)O",
        cid=2724385,
        collection="Natural product",
    ),
    "doxorubicin": molecule(
        "Doxorubicin",
        "C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C(=CC=C5)OC)O)(C(=O)CO)O)N)O",
        cid=31703,
        collection="Natural product",
    ),
    "epibatidine": molecule(
        "Epibatidine",
        "C1C[C@@H]2[C@H](C[C@H]1N2)C3=CN=C(C=C3)Cl",
        cid=854023,
        collection="Natural product",
    ),
    "epigallocatechin_gallate": molecule(
        "Epigallocatechin gallate",
        "C1[C@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C(=C3)O)O)O)OC(=O)C4=CC(=C(C(=C4)O)O)O",
        cid=65064,
        collection="Natural product",
        aliases=("egcg",),
    ),
    "ergotamine": molecule(
        "Ergotamine",
        "C[C@@]1(C(=O)N2[C@H](C(=O)N3CCC[C@H]3[C@@]2(O1)O)CC4=CC=CC=C4)NC(=O)[C@H]5CN([C@@H]6CC7=CNC8=CC=CC(=C78)C6=C5)C",
        cid=8223,
        collection="Natural product",
    ),
    "erythromycin": molecule(
        "Erythromycin",
        "CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
        cid=12560,
        collection="Natural product",
    ),
    "galantamine": molecule(
        "Galantamine",
        "CN1CC[C@@]23C=C[C@@H](C[C@@H]2OC4=C(C=CC(=C34)C1)OC)O",
        cid=9651,
        collection="Natural product",
    ),
    "genistein": molecule(
        "Genistein",
        "C1=CC(=CC=C1C2=COC3=CC(=CC(=C3C2=O)O)O)O",
        cid=5280961,
        collection="Natural product",
    ),
    "harmine": molecule(
        "Harmine", "CC1=NC=CC2=C1NC3=C2C=CC(=C3)OC", cid=5280953, collection="Natural product"
    ),
    "harmaline": molecule(
        "Harmaline", "CC1=NCCC2=C1NC3=C2C=CC(=C3)OC", cid=3564, collection="Natural product"
    ),
    "huperzine_a": molecule(
        "Huperzine A",
        "C/C=C/1\\[C@@H]2CC3=C([C@]1(CC(=C2)C)N)C=CC(=O)N3",
        cid=854026,
        collection="Natural product",
    ),
    "indirubin": molecule(
        "Indirubin",
        "C1=CC=C2C(=C1)C(=C(N2)O)C3=NC4=CC=CC=C4C3=O",
        cid=10177,
        collection="Natural product",
    ),
    "kavain": molecule(
        "Kavain",
        "COC1=CC(=O)O[C@H](C1)/C=C/C2=CC=CC=C2",
        cid=5281565,
        collection="Natural product",
    ),
    "lovastatin": molecule(
        "Lovastatin",
        "CC[C@H](C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C",
        cid=53232,
        collection="Natural product",
    ),
    "menthol": molecule("Menthol", "CC1CCC(C(C1)O)C(C)C", cid=1254, collection="Natural product"),
    "mescaline": molecule(
        "Mescaline",
        "COC1=CC(=CC(=C1OC)OC)CCN",
        cid=4076,
        collection="Natural product",
    ),
    "morphine": molecule(
        "Morphine",
        "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O",
        cid=5288826,
        collection="Natural product",
    ),
    "mycophenolic_acid": molecule(
        "Mycophenolic acid",
        "CC1=C2COC(=O)C2=C(C(=C1OC)C/C=C(\\C)/CCC(=O)O)O",
        cid=446541,
        collection="Natural product",
    ),
    "noscapine": molecule(
        "Noscapine",
        "CN1CCC2=CC3=C(C(=C2[C@@H]1[C@@H]4C5=C(C(=C(C=C5)OC)OC)C(=O)O4)OC)OCO3",
        cid=275196,
        collection="Natural product",
    ),
    "paclitaxel": molecule(
        "Paclitaxel",
        "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
        cid=36314,
        collection="Natural product",
        aliases=("taxol",),
    ),
    "penicillin_g": molecule(
        "Penicillin G",
        "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
        cid=5904,
        collection="Natural product",
    ),
    "psilocybin": molecule(
        "Psilocybin",
        "CN(C)CCC1=CNC2=C1C(=CC=C2)OP(=O)(O)O",
        cid=10624,
        collection="Natural product",
    ),
    "psilocin": molecule(
        "Psilocin",
        "CN(C)CCC1=CNC2=C1C(=CC=C2)O",
        cid=4980,
        collection="Natural product",
    ),
    "quercetin": molecule(
        "Quercetin",
        "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
        cid=5280343,
        collection="Natural product",
    ),
    "quinine": molecule(
        "Quinine",
        "COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O",
        cid=3034034,
        collection="Natural product",
    ),
    "reserpine": molecule(
        "Reserpine",
        "CO[C@H]1[C@@H](C[C@@H]2CN3CCC4=C([C@H]3C[C@@H]2[C@@H]1C(=O)OC)NC5=C4C=CC(=C5)OC)OC(=O)C6=CC(=C(C(=C6)OC)OC)OC",
        cid=5770,
        collection="Natural product",
    ),
    "resveratrol": molecule(
        "Resveratrol",
        "C1=CC(=CC=C1/C=C/C2=CC(=CC(=C2)O)O)O",
        cid=445154,
        collection="Natural product",
    ),
    "salvinorin_a": molecule(
        "Salvinorin A",
        "CC(=O)O[C@H]1C[C@H]([C@@]2(CC[C@H]3C(=O)O[C@@H](C[C@@]3([C@H]2C1=O)C)C4=COC=C4)C)C(=O)OC",
        cid=128563,
        collection="Natural product",
    ),
    "scopolamine": molecule(
        "Scopolamine",
        "CN1[C@@H]2CC(C[C@H]1[C@H]3[C@@H]2O3)OC(=O)[C@H](CO)C4=CC=CC=C4",
        cid=3000322,
        collection="Natural product",
    ),
    "shikonin": molecule(
        "Shikonin",
        "CC(=CC[C@H](C1=CC(=O)C2=C(C=CC(=C2C1=O)O)O)O)C",
        cid=479503,
        collection="Natural product",
    ),
    "silybin": molecule(
        "Silybin",
        "COC1=C(C=CC(=C1)[C@@H]2[C@H](OC3=C(O2)C=C(C=C3)[C@@H]4[C@H](C(=O)C5=C(C=C(C=C5O4)O)O)O)CO)O",
        cid=31553,
        collection="Natural product",
        aliases=("silibinin",),
    ),
    "sirolimus": molecule(
        "Sirolimus",
        "C[C@@H]1CC[C@H]2C[C@@H](/C(=C/C=C/C=C/[C@H](C[C@H](C(=O)[C@@H]([C@@H](/C(=C/[C@H](C(=O)C[C@H](OC(=O)[C@@H]3CCCCN3C(=O)C(=O)[C@@]1(O2)O)[C@H](C)C[C@@H]4CC[C@H]([C@@H](C4)OC)O)C)/C)O)OC)C)C)/C)OC",
        cid=5284616,
        collection="Natural product",
    ),
    "tacrolimus": molecule(
        "Tacrolimus",
        "C[C@@H]1C[C@@H]([C@@H]2[C@H](C[C@H]([C@@](O2)(C(=O)C(=O)N3CCCC[C@H]3C(=O)O[C@@H]([C@@H]([C@H](CC(=O)[C@@H](/C=C(/C1)\\C)CC=C)O)C)/C(=C/[C@@H]4CC[C@H]([C@@H](C4)OC)O)/C)O)C)OC)OC",
        cid=445643,
        collection="Natural product",
    ),
    "taxifolin": molecule(
        "Taxifolin",
        "C1=CC(=C(C=C1[C@@H]2[C@H](C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
        cid=439533,
        collection="Natural product",
    ),
    "tetrandrine": molecule(
        "Tetrandrine",
        "CN1CCC2=CC(=C3C=C2[C@@H]1CC4=CC=C(C=C4)OC5=C(C=CC(=C5)C[C@H]6C7=C(O3)C(=C(C=C7CCN6C)OC)OC)OC)OC",
        cid=73078,
        collection="Natural product",
    ),
    "theobromine": molecule(
        "Theobromine",
        "CN1C=NC2=C1C(=O)NC(=O)N2C",
        cid=5429,
        collection="Natural product",
    ),
    "theophylline": molecule(
        "Theophylline",
        "CN1C2=C(C(=O)N(C1=O)C)NC=N2",
        cid=2153,
        collection="Natural product",
    ),
    "thymol": molecule("Thymol", "CC1=CC(=C(C=C1)C(C)C)O", cid=6989, collection="Natural product"),
    "ursolic_acid": molecule(
        "Ursolic acid",
        "C[C@@H]1CC[C@@]2(CC[C@@]3(C(=CC[C@H]4[C@]3(CC[C@@H]5[C@@]4(CC[C@@H](C5(C)C)O)C)C)[C@@H]2[C@H]1C)C)C(=O)O",
        cid=64945,
        collection="Natural product",
    ),
    "vinblastine": molecule(
        "Vinblastine",
        "CC[C@@]1(C[C@H]2C[C@@](C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)[C@]78CCN9[C@H]7[C@@](C=CC9)([C@H]([C@@]([C@@H]8N6C)(C(=O)OC)O)OC(=O)C)CC)OC)C(=O)OC)O",
        cid=13342,
        collection="Natural product",
    ),
    "vincristine": molecule(
        "Vincristine",
        "CC[C@@]1(C[C@@H]2C[C@@](C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)[C@]78CCN9[C@H]7[C@@](C=CC9)([C@H]([C@@]([C@@H]8N6C=O)(C(=O)OC)O)OC(=O)C)CC)OC)C(=O)OC)O",
        cid=5978,
        collection="Natural product",
    ),
    "withaferin_a": molecule(
        "Withaferin A",
        "CC1=C(C(=O)O[C@H](C1)[C@@H](C)[C@H]2CC[C@@H]3[C@@]2(CC[C@H]4[C@H]3C[C@@H]5[C@]6([C@@]4(C(=O)C=C[C@@H]6O)C)O5)C)CO",
        cid=265237,
        collection="Natural product",
    ),
    "yohimbine": molecule(
        "Yohimbine",
        "COC(=O)[C@H]1[C@H](CC[C@@H]2[C@@H]1C[C@H]3C4=C(CCN3C2)C5=CC=CC=C5N4)O",
        cid=8969,
        collection="Natural product",
    ),
    "abacavir": molecule(
        "Abacavir",
        "C1CC1NC2=C3C(=NC(=N2)N)N(C=N3)[C@@H]4C[C@@H](C=C4)CO",
        cid=441300,
        collection="FDA-approved drug",
    ),
    "acyclovir": molecule(
        "Acyclovir", "C1=NC2=C(N1COCCO)N=C(NC2=O)N", cid=135398513, collection="FDA-approved drug"
    ),
    "albuterol": molecule(
        "Albuterol",
        "CC(C)(C)NCC(C1=CC(=C(C=C1)O)CO)O",
        cid=2083,
        collection="FDA-approved drug",
        aliases=("salbutamol",),
    ),
    "allopurinol": molecule(
        "Allopurinol", "C1=NNC2=C1C(=O)NC=N2", cid=135401907, collection="FDA-approved drug"
    ),
    "alprazolam": molecule(
        "Alprazolam",
        "CC1=NN=C2N1C3=C(C=C(C=C3)Cl)C(=NC2)C4=CC=CC=C4",
        cid=2118,
        collection="FDA-approved drug",
    ),
    "amlodipine": molecule(
        "Amlodipine",
        "CCOC(=O)C1=C(NC(=C(C1C2=CC=CC=C2Cl)C(=O)OC)C)COCCN",
        cid=2162,
        collection="FDA-approved drug",
    ),
    "amoxicillin": molecule(
        "Amoxicillin",
        "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)[C@@H](C3=CC=C(C=C3)O)N)C(=O)O)C",
        cid=33613,
        collection="FDA-approved drug",
    ),
    "apixaban": molecule(
        "Apixaban",
        "COC1=CC=C(C=C1)N2C3=C(CCN(C3=O)C4=CC=C(C=C4)N5CCCCC5=O)C(=N2)C(=O)N",
        cid=10182969,
        collection="FDA-approved drug",
    ),
    "aripiprazole": molecule(
        "Aripiprazole",
        "C1CC(=O)NC2=C1C=CC(=C2)OCCCCN3CCN(CC3)C4=C(C(=CC=C4)Cl)Cl",
        cid=60795,
        collection="FDA-approved drug",
    ),
    "atorvastatin": molecule(
        "Atorvastatin",
        "CC(C)C1=C(C(=C(N1CC[C@H](C[C@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
        cid=60823,
        collection="FDA-approved drug",
    ),
    "azithromycin": molecule(
        "Azithromycin",
        "CC[C@@H]1[C@@]([C@@H]([C@H](N(C[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O",
        cid=447043,
        collection="FDA-approved drug",
    ),
    "bupropion": molecule(
        "Bupropion", "CC(C(=O)C1=CC(=CC=C1)Cl)NC(C)(C)C", cid=444, collection="FDA-approved drug"
    ),
    "buprenorphine": molecule(
        "Buprenorphine",
        "C[C@]([C@H]1C[C@@]23CC[C@@]1([C@H]4[C@@]25CCN([C@@H]3CC6=C5C(=C(C=C6)O)O4)CC7CC7)OC)(C(C)(C)C)O",
        cid=644073,
        collection="FDA-approved drug",
    ),
    "celecoxib": molecule(
        "Celecoxib",
        "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F",
        cid=2662,
        collection="FDA-approved drug",
    ),
    "cetirizine": molecule(
        "Cetirizine",
        "C1CN(CCN1CCOCC(=O)O)C(C2=CC=CC=C2)C3=CC=C(C=C3)Cl",
        cid=2678,
        collection="FDA-approved drug",
    ),
    "ciprofloxacin": molecule(
        "Ciprofloxacin",
        "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O",
        cid=2764,
        collection="FDA-approved drug",
    ),
    "clopidogrel": molecule(
        "Clopidogrel",
        "COC(=O)[C@H](C1=CC=CC=C1Cl)N2CCC3=C(C2)C=CS3",
        cid=60606,
        collection="FDA-approved drug",
    ),
    "dapagliflozin": molecule(
        "Dapagliflozin",
        "CCOC1=CC=C(C=C1)CC2=C(C=CC(=C2)[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)O)Cl",
        cid=9887712,
        collection="FDA-approved drug",
    ),
    "dasatinib": molecule(
        "Dasatinib",
        "CC1=C(C(=CC=C1)Cl)NC(=O)C2=CN=C(S2)NC3=CC(=NC(=N3)C)N4CCN(CC4)CCO",
        cid=3062316,
        collection="FDA-approved drug",
    ),
    "dexamethasone": molecule(
        "Dexamethasone",
        "C[C@@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@@]4([C@]3([C@H](C[C@@]2([C@]1(C(=O)CO)O)C)O)F)C",
        cid=5743,
        collection="FDA-approved drug",
    ),
    "diazepam": molecule(
        "Diazepam",
        "CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3",
        cid=3016,
        collection="FDA-approved drug",
    ),
    "diclofenac": molecule(
        "Diclofenac",
        "C1=CC=C(C(=C1)CC(=O)O)NC2=C(C=CC=C2Cl)Cl",
        cid=3033,
        collection="FDA-approved drug",
    ),
    "donepezil": molecule(
        "Donepezil",
        "COC1=C(C=C2C(=C1)CC(C2=O)CC3CCN(CC3)CC4=CC=CC=C4)OC",
        cid=3152,
        collection="FDA-approved drug",
    ),
    "duloxetine": molecule(
        "Duloxetine",
        "CNCC[C@@H](C1=CC=CS1)OC2=CC=CC3=CC=CC=C32",
        cid=60835,
        collection="FDA-approved drug",
    ),
    "empagliflozin": molecule(
        "Empagliflozin",
        "C1COC[C@H]1OC2=CC=C(C=C2)CC3=C(C=CC(=C3)[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)Cl",
        cid=11949646,
        collection="FDA-approved drug",
    ),
    "escitalopram": molecule(
        "Escitalopram",
        "CN(C)CCC[C@@]1(C2=C(CO1)C=C(C=C2)C#N)C3=CC=C(C=C3)F",
        cid=146570,
        collection="FDA-approved drug",
    ),
    "esomeprazole": molecule(
        "Esomeprazole",
        "CC1=CN=C(C(=C1OC)C)C[S@](=O)C2=NC3=C(N2)C=C(C=C3)OC",
        cid=9568614,
        collection="FDA-approved drug",
    ),
    "famotidine": molecule(
        "Famotidine",
        "C1=C(N=C(S1)N=C(N)N)CSCC/C(=N/S(=O)(=O)N)/N",
        cid=5702160,
        collection="FDA-approved drug",
    ),
    "finasteride": molecule(
        "Finasteride",
        "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2C(=O)NC(C)(C)C)CC[C@@H]4[C@@]3(C=CC(=O)N4)C",
        cid=57363,
        collection="FDA-approved drug",
    ),
    "fluconazole": molecule(
        "Fluconazole",
        "C1=CC(=C(C=C1F)F)C(CN2C=NC=N2)(CN3C=NC=N3)O",
        cid=3365,
        collection="FDA-approved drug",
    ),
    "fluoxetine": molecule(
        "Fluoxetine",
        "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F",
        cid=3386,
        collection="FDA-approved drug",
    ),
    "furosemide": molecule(
        "Furosemide",
        "C1=COC(=C1)CNC2=CC(=C(C=C2C(=O)O)S(=O)(=O)N)Cl",
        cid=3440,
        collection="FDA-approved drug",
    ),
    "gabapentin": molecule(
        "Gabapentin", "C1CCC(CC1)(CC(=O)O)CN", cid=3446, collection="FDA-approved drug"
    ),
    "hydrochlorothiazide": molecule(
        "Hydrochlorothiazide",
        "C1NC2=CC(=C(C=C2S(=O)(=O)N1)S(=O)(=O)N)Cl",
        cid=3639,
        collection="FDA-approved drug",
    ),
    "hydroxychloroquine": molecule(
        "Hydroxychloroquine",
        "CCN(CCCC(C)NC1=C2C=CC(=CC2=NC=C1)Cl)CCO",
        cid=3652,
        collection="FDA-approved drug",
    ),
    "imatinib": molecule(
        "Imatinib",
        "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
        cid=5291,
        collection="FDA-approved drug",
    ),
    "lamotrigine": molecule(
        "Lamotrigine",
        "C1=CC(=C(C(=C1)Cl)Cl)C2=C(N=C(N=N2)N)N",
        cid=3878,
        collection="FDA-approved drug",
    ),
    "levetiracetam": molecule(
        "Levetiracetam",
        "CC[C@@H](C(=O)N)N1CCCC1=O",
        cid=5284583,
        collection="FDA-approved drug",
    ),
    "levodopa": molecule(
        "Levodopa",
        "C1=CC(=C(C=C1C[C@@H](C(=O)O)N)O)O",
        cid=6047,
        collection="FDA-approved drug",
        aliases=("L-DOPA",),
    ),
    "lisinopril": molecule(
        "Lisinopril",
        "C1C[C@H](N(C1)C(=O)[C@H](CCCCN)N[C@@H](CCC2=CC=CC=C2)C(=O)O)C(=O)O",
        cid=5362119,
        collection="FDA-approved drug",
    ),
    "loratadine": molecule(
        "Loratadine",
        "CCOC(=O)N1CCC(=C2C3=C(CCC4=C2N=CC=C4)C=C(C=C3)Cl)CC1",
        cid=3957,
        collection="FDA-approved drug",
    ),
    "losartan": molecule(
        "Losartan",
        "CCCCC1=NC(=C(N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4)CO)Cl",
        cid=3961,
        collection="FDA-approved drug",
    ),
    "meloxicam": molecule(
        "Meloxicam",
        "CC1=CN=C(S1)NC(=O)C2=C(C3=CC=CC=C3S(=O)(=O)N2C)O",
        cid=54677470,
        collection="FDA-approved drug",
    ),
    "metformin": molecule(
        "Metformin", "CN(C)C(=N)N=C(N)N", cid=4091, collection="FDA-approved drug"
    ),
    "methadone": molecule(
        "Methadone",
        "CCC(=O)C(CC(C)N(C)C)(C1=CC=CC=C1)C2=CC=CC=C2",
        cid=4095,
        collection="FDA-approved drug",
    ),
    "methylphenidate": molecule(
        "Methylphenidate",
        "COC(=O)C(C1CCCCN1)C2=CC=CC=C2",
        cid=4158,
        collection="FDA-approved drug",
    ),
    "metoprolol": molecule(
        "Metoprolol",
        "CC(C)NCC(COC1=CC=C(C=C1)CCOC)O",
        cid=4171,
        collection="FDA-approved drug",
    ),
    "minoxidil": molecule(
        "Minoxidil", "C1CCN(CC1)C2=NC(=N)N(C(=C2)N)O", cid=4201, collection="FDA-approved drug"
    ),
    "montelukast": molecule(
        "Montelukast",
        "CC(C)(C1=CC=CC=C1CC[C@H](C2=CC=CC(=C2)/C=C/C3=NC4=C(C=CC(=C4)Cl)C=C3)SCC5(CC5)CC(=O)O)O",
        cid=5281040,
        collection="FDA-approved drug",
    ),
    "naloxone": molecule(
        "Naloxone",
        "C=CCN1CC[C@]23[C@@H]4C(=O)CC[C@]2([C@H]1CC5=C3C(=C(C=C5)O)O4)O",
        cid=5284596,
        collection="FDA-approved drug",
    ),
    "naproxen": molecule(
        "Naproxen",
        "C[C@@H](C1=CC2=C(C=C1)C=C(C=C2)OC)C(=O)O",
        cid=156391,
        collection="FDA-approved drug",
    ),
    "olanzapine": molecule(
        "Olanzapine",
        "CC1=CC2=C(S1)NC3=CC=CC=C3N=C2N4CCN(CC4)C",
        cid=135398745,
        collection="FDA-approved drug",
    ),
    "omeprazole": molecule(
        "Omeprazole",
        "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC",
        cid=4594,
        collection="FDA-approved drug",
    ),
    "ondansetron": molecule(
        "Ondansetron",
        "CC1=NC=CN1CC2CCC3=C(C2=O)C4=CC=CC=C4N3C",
        cid=4595,
        collection="FDA-approved drug",
    ),
    "oseltamivir": molecule(
        "Oseltamivir",
        "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)OCC",
        cid=65028,
        collection="FDA-approved drug",
    ),
    "oxycodone": molecule(
        "Oxycodone",
        "CN1CC[C@]23[C@@H]4C(=O)CC[C@]2([C@H]1CC5=C3C(=C(C=C5)OC)O4)O",
        cid=5284603,
        collection="FDA-approved drug",
    ),
    "pantoprazole": molecule(
        "Pantoprazole",
        "COC1=C(C(=NC=C1)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC(F)F)OC",
        cid=4679,
        collection="FDA-approved drug",
    ),
    "pregabalin": molecule(
        "Pregabalin",
        "CC(C)C[C@@H](CC(=O)O)CN",
        cid=5486971,
        collection="FDA-approved drug",
    ),
    "prednisone": molecule(
        "Prednisone",
        "C[C@]12CC(=O)[C@H]3[C@H]([C@@H]1CC[C@@]2(C(=O)CO)O)CCC4=CC(=O)C=C[C@]34C",
        cid=5865,
        collection="FDA-approved drug",
    ),
    "quetiapine": molecule(
        "Quetiapine",
        "C1CN(CCN1CCOCCO)C2=NC3=CC=CC=C3SC4=CC=CC=C42",
        cid=5002,
        collection="FDA-approved drug",
    ),
    "remdesivir": molecule(
        "Remdesivir",
        "CCC(CC)COC(=O)[C@H](C)N[P@](=O)(OC[C@@H]1[C@H]([C@H]([C@](O1)(C#N)C2=CC=C3N2N=CN=C3N)O)O)OC4=CC=CC=C4",
        cid=121304016,
        collection="FDA-approved drug",
    ),
    "rivaroxaban": molecule(
        "Rivaroxaban",
        "C1COCC(=O)N1C2=CC=C(C=C2)N3C[C@@H](OC3=O)CNC(=O)C4=CC=C(S4)Cl",
        cid=9875401,
        collection="FDA-approved drug",
    ),
    "rosuvastatin": molecule(
        "Rosuvastatin",
        "CC(C)C1=NC(=NC(=C1/C=C/[C@H](C[C@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C",
        cid=446157,
        collection="FDA-approved drug",
    ),
    "sertraline": molecule(
        "Sertraline",
        "CN[C@H]1CC[C@H](C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl",
        cid=68617,
        collection="FDA-approved drug",
    ),
    "sildenafil": molecule(
        "Sildenafil",
        "CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C",
        cid=135398744,
        collection="FDA-approved drug",
    ),
    "simvastatin": molecule(
        "Simvastatin",
        "CCC(C)(C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C",
        cid=54454,
        collection="FDA-approved drug",
    ),
    "sitagliptin": molecule(
        "Sitagliptin",
        "C1CN2C(=NN=C2C(F)(F)F)CN1C(=O)C[C@@H](CC3=CC(=C(C=C3F)F)F)N",
        cid=4369359,
        collection="FDA-approved drug",
    ),
    "sofosbuvir": molecule(
        "Sofosbuvir",
        "C[C@@H](C(=O)OC(C)C)N[P@](=O)(OC[C@@H]1[C@H]([C@@]([C@@H](O1)N2C=CC(=O)NC2=O)(C)F)O)OC3=CC=CC=C3",
        cid=45375808,
        collection="FDA-approved drug",
    ),
    "spironolactone": molecule(
        "Spironolactone",
        "CC(=O)S[C@@H]1CC2=CC(=O)CC[C@@]2([C@@H]3[C@@H]1[C@@H]4CC[C@]5([C@]4(CC3)C)CCC(=O)O5)C",
        cid=5833,
        collection="FDA-approved drug",
    ),
    "tamoxifen": molecule(
        "Tamoxifen",
        "CC/C(=C(\\C1=CC=CC=C1)/C2=CC=C(C=C2)OCCN(C)C)/C3=CC=CC=C3",
        cid=2733526,
        collection="FDA-approved drug",
    ),
    "tadalafil": molecule(
        "Tadalafil",
        "CN1CC(=O)N2[C@@H](C1=O)CC3=C([C@H]2C4=CC5=C(C=C4)OCO5)NC6=CC=CC=C36",
        cid=110635,
        collection="FDA-approved drug",
    ),
    "tenofovir": molecule(
        "Tenofovir",
        "C[C@H](CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O",
        cid=464205,
        collection="FDA-approved drug",
    ),
    "topiramate": molecule(
        "Topiramate",
        "CC1(O[C@@H]2CO[C@@]3([C@H]([C@@H]2O1)OC(O3)(C)C)COS(=O)(=O)N)C",
        cid=5284627,
        collection="FDA-approved drug",
    ),
    "tramadol": molecule(
        "Tramadol",
        "CN(C)C[C@H]1CCCC[C@@]1(C2=CC(=CC=C2)OC)O",
        cid=33741,
        collection="FDA-approved drug",
    ),
    "trazodone": molecule(
        "Trazodone",
        "C1CN(CCN1CCCN2C(=O)N3C=CC=CC3=N2)C4=CC(=CC=C4)Cl",
        cid=5533,
        collection="FDA-approved drug",
    ),
    "valacyclovir": molecule(
        "Valacyclovir",
        "CC(C)[C@@H](C(=O)OCCOCN1C=NC2=C1N=C(NC2=O)N)N",
        cid=135398742,
        collection="FDA-approved drug",
    ),
    "valproic_acid": molecule(
        "Valproic acid", "CCCC(CCC)C(=O)O", cid=3121, collection="FDA-approved drug"
    ),
    "venlafaxine": molecule(
        "Venlafaxine",
        "CN(C)CC(C1=CC=C(C=C1)OC)C2(CCCCC2)O",
        cid=5656,
        collection="FDA-approved drug",
    ),
    "warfarin": molecule(
        "Warfarin",
        "CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O",
        cid=54678486,
        collection="FDA-approved drug",
    ),
    "zolpidem": molecule(
        "Zolpidem",
        "CC1=CC=C(C=C1)C2=C(N3C=C(C=CC3=N2)C)CC(=O)N(C)C",
        cid=5732,
        collection="FDA-approved drug",
    ),
}

KNOWN_MOLECULE_KEYS = tuple(
    sorted(
        KNOWN_MOLECULES,
        key=lambda key: (KNOWN_MOLECULES[key].collection, KNOWN_MOLECULES[key].name),
    )
)
