import numpy as np
from rdkit import Chem

from amsr.conf import (
    GetConformer,
    _build_nonbonded_clash_data,
    _has_serious_nonbonded_clash,
)
from amsr.decode import ToMol


def test_vdw_thresholds_treat_one_four_pairs_less_strictly():
    butane = Chem.MolFromSmiles("CCCC")
    pairs, squared_thresholds = _build_nonbonded_clash_data(butane)
    assert pairs.tolist() == [[0, 3]]
    assert np.sqrt(squared_thresholds[0]) == 0.65 * (1.70 + 1.70)

    pentane = Chem.MolFromSmiles("CCCCC")
    pairs, squared_thresholds = _build_nonbonded_clash_data(pentane)
    threshold_by_pair = {
        tuple(pair): np.sqrt(threshold)
        for pair, threshold in zip(pairs, squared_thresholds, strict=True)
    }
    assert threshold_by_pair[0, 4] == 0.75 * (1.70 + 1.70)


def test_vdw_overlap_check_uses_precomputed_threshold():
    mol = Chem.MolFromSmiles("CCCC")
    clash_data = _build_nonbonded_clash_data(mol)
    coordinates = np.zeros((4, 3), dtype=np.float64)
    coordinates[3, 0] = 2.20
    assert _has_serious_nonbonded_clash(coordinates, clash_data)
    coordinates[3, 0] = 2.22
    assert not _has_serious_nonbonded_clash(coordinates, clash_data)


def test_conformer_selection_rejects_clashing_chembl3d_candidate():
    # CHEMBL237460_0 conformer 130, the worst ChEMBL3D geometry-audit case.
    amsr_string = r"cccn+cccc6cc6@__N__C<<(6).[CN]......\_N\>C`\_C^^C^\C/^C__5......[CN]"
    dihedrals = {}
    topology = ToMol(amsr_string, stringent=True, dihedral=dihedrals)

    decoded = GetConformer(topology, dihedral=dihedrals, max_confs=100)

    assert not _has_serious_nonbonded_clash(
        np.asarray(decoded.GetConformer().GetPositions()),
        _build_nonbonded_clash_data(decoded),
    )
