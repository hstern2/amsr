use gamma::graph::{DefaultGraph, Graph};
use gamma::matching::{Pairing, maximum_matching};
use crate::molecule::Molecule;
use crate::atom::Atom;
use crate::bond::Bond;
use std::collections::HashSet;

/// Build a pi subgraph for aromatic atoms
/// This creates a graph where nodes are aromatic atoms that can participate in π bonds
/// and edges represent potential π bonds between them
pub fn build_pi_subgraph(mol: &Molecule, aromatic_atoms: &HashSet<usize>) -> DefaultGraph {
    let mut pi_graph = DefaultGraph::new();

    for &i in aromatic_atoms {
        pi_graph.add_node(i).expect("add node");
    }

    for &i in aromatic_atoms {
        for &neighbor_idx in &mol.get_neighbors(i) {
            if aromatic_atoms.contains(&neighbor_idx) && i < neighbor_idx {
                pi_graph.add_edge(i, neighbor_idx).expect("add edge");
            }
        }
    }

    pi_graph
}

/// Kekulize aromatic bonds in a molecule
/// This assigns single/double bonds to aromatic systems using maximum matching
pub fn kekulize_molecule(mol: &mut Molecule, aromatic_atoms: &HashSet<usize>) -> Result<(), Box<dyn std::error::Error>> {
    let pi_graph = build_pi_subgraph(mol, aromatic_atoms);

    if pi_graph.order() == 0 {
        // No aromatic atoms, nothing to do
        return Ok(());
    }

    // Find maximum matching in the pi subgraph
    let mut pairing = Pairing::new();
    maximum_matching(&pi_graph, &mut pairing);

    // Check if we have a perfect matching (all aromatic atoms matched)
    if pairing.order() != pi_graph.order() {
        return Err("Cannot kekulize aromatic system - no valid Kekulé structure found".into());
    }

    // Apply the matching: matched edges become double bonds, others remain single
    let matched_edges: std::collections::HashSet<_> = pairing.edges().collect();

    // Collect all aromatic atom pairs (i, neighbor_idx)
    let mut aromatic_pairs = Vec::new();
    for &i in aromatic_atoms {
        for &neighbor_idx in &mol.get_neighbors(i) {
            if aromatic_atoms.contains(&neighbor_idx) && i < neighbor_idx {
                aromatic_pairs.push((i, neighbor_idx));
            }
        }
    }

    // Now mutate bonds
    for (i, neighbor_idx) in aromatic_pairs {
        let edge = if i < neighbor_idx { (i, neighbor_idx) } else { (neighbor_idx, i) };
        let is_double = matched_edges.contains(&edge);
        if let Some(bond) = mol.get_bond_mut(i, neighbor_idx) {
            if is_double {
                bond.bond_type = crate::bond::BondType::Double;
            } else {
                bond.bond_type = crate::bond::BondType::Single;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles_decode::decode_smiles;

    #[test]
    fn test_benzene_kekulization() {
        let mut mol = decode_smiles("c1ccccc1").unwrap();
        let aromatic_atoms = HashSet::from([0, 1, 2, 3, 4, 5]);
        assert!(kekulize_molecule(&mut mol, &aromatic_atoms).is_ok());

        // After kekulization, should have alternating single/double bonds
        let bonds = &mol.bonds;
        assert_eq!(bonds.len(), 6);

        // Check that we have 3 double bonds and 3 single bonds
        let double_bonds = bonds.iter().filter(|(_, _, bond)| bond.bond_type == crate::bond::BondType::Double).count();
        let single_bonds = bonds.iter().filter(|(_, _, bond)| bond.bond_type == crate::bond::BondType::Single).count();

        assert_eq!(double_bonds, 3);
        assert_eq!(single_bonds, 3);
    }

    #[test]
    fn test_non_aromatic_molecule() {
        let mut mol = decode_smiles("CCO").unwrap();
        let aromatic_atoms = HashSet::from([0, 1, 2]);
        assert!(kekulize_molecule(&mut mol, &aromatic_atoms).is_ok());
        // Should not change anything for non-aromatic molecules
    }
}