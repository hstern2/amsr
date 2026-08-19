import numpy as np
from rdkit import Chem

from amsr.conf import GetConformer
from amsr.decode import ToMol


def _serious_nonbonded_clashes(mol):
    coordinates = np.asarray(mol.GetConformer().GetPositions())
    euclidean = np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=-1)
    graph_distances = Chem.GetDistanceMatrix(mol)
    return int(np.count_nonzero(np.triu((graph_distances > 2) & (euclidean < 1.0), k=1)))


def test_conformer_selection_rejects_clashing_chembl3d_candidate():
    # CHEMBL237460_0 conformer 130, the worst ChEMBL3D geometry-audit case.
    amsr_string = r"cccn+cccc6cc6@__N__C<<(6).[CN]......\_N\>C`\_C^^C^\C/^C__5......[CN]"
    dihedrals = {}
    topology = ToMol(amsr_string, stringent=True, dihedral=dihedrals)

    decoded = GetConformer(topology, dihedral=dihedrals, max_confs=100)

    assert _serious_nonbonded_clashes(decoded) == 0
