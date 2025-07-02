use crate::errors::{AMSRResult, AMSRError};
use crate::atom::Atom;
use crate::bond::Bond;

use std::collections::{HashMap, HashSet, VecDeque};

pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<(usize, usize, Bond)>, // (atom1_idx, atom2_idx, bond)
    pub adjacency_list: HashMap<usize, Vec<usize>>, // atom_idx -> list of connected atom indices
}

impl Molecule {
    pub fn new() -> Self {
        Molecule {
            atoms: Vec::new(),
            bonds: Vec::new(),
            adjacency_list: HashMap::new(),
        }
    }
    

    
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(atom);
        self.adjacency_list.insert(idx, Vec::new());
        idx
    }
    
    pub fn add_bond(&mut self, atom1_idx: usize, atom2_idx: usize, bond: Bond) -> AMSRResult<()> {
        if atom1_idx >= self.atoms.len() || atom2_idx >= self.atoms.len() {
            return Err(AMSRError::GraphError("Invalid atom indices".to_string()));
        }
        
        // Add bond
        self.bonds.push((atom1_idx, atom2_idx, bond));
        
        // Update adjacency list
        self.adjacency_list.get_mut(&atom1_idx).unwrap().push(atom2_idx);
        self.adjacency_list.get_mut(&atom2_idx).unwrap().push(atom1_idx);
        
        // Update atom neighbor counts
        {
            let (left, right) = self.atoms.split_at_mut(atom1_idx + 1);
            let atom1 = &mut left[atom1_idx];
            let atom2 = &mut right[atom2_idx - atom1_idx - 1];
            atom1.add_bond_to(atom2);
        }
        
        Ok(())
    }
    
    pub fn get_neighbors(&self, atom_idx: usize) -> Vec<usize> {
        self.adjacency_list.get(&atom_idx).cloned().unwrap_or_default()
    }
    
    pub fn get_bond(&self, atom1_idx: usize, atom2_idx: usize) -> Option<&Bond> {
        for (a1, a2, bond) in &self.bonds {
            if (*a1 == atom1_idx && *a2 == atom2_idx) || (*a1 == atom2_idx && *a2 == atom1_idx) {
                return Some(bond);
            }
        }
        None
    }
    
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }
    
    pub fn num_bonds(&self) -> usize {
        self.bonds.len()
    }
    
    pub fn is_connected(&self) -> bool {
        if self.atoms.is_empty() {
            return true;
        }
        
        let mut visited = HashSet::new();
        let mut queue = VecDeque::new();
        queue.push_back(0);
        visited.insert(0);
        
        while let Some(current) = queue.pop_front() {
            for &neighbor in &self.get_neighbors(current) {
                if !visited.contains(&neighbor) {
                    visited.insert(neighbor);
                    queue.push_back(neighbor);
                }
            }
        }
        
        visited.len() == self.atoms.len()
    }
    
    pub fn find_rings(&self) -> Vec<Vec<usize>> {
        let mut rings = Vec::new();
        let mut visited = HashSet::new();
        
        for start in 0..self.atoms.len() {
            if !visited.contains(&start) {
                self.dfs_find_rings(start, start, &mut Vec::new(), &mut visited, &mut rings);
            }
        }
        
        rings
    }
    
    fn dfs_find_rings(&self, current: usize, start: usize, path: &mut Vec<usize>, 
                     visited: &mut HashSet<usize>, rings: &mut Vec<Vec<usize>>) {
        path.push(current);
        visited.insert(current);
        
        for &neighbor in &self.get_neighbors(current) {
            if neighbor == start && path.len() > 2 {
                // Found a ring
                rings.push(path.clone());
            } else if !visited.contains(&neighbor) {
                self.dfs_find_rings(neighbor, start, path, visited, rings);
            }
        }
        
        path.pop();
        visited.remove(&current);
    }
    

    

}

impl Default for Molecule {
    fn default() -> Self {
        Molecule::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_molecule() {
        let mut mol = Molecule::new();
        let c1 = mol.add_atom(Atom::new("C").unwrap());
        let c2 = mol.add_atom(Atom::new("C").unwrap());
        let o = mol.add_atom(Atom::new("O").unwrap());
        
        mol.add_bond(c1, c2, Bond::single()).unwrap();
        mol.add_bond(c2, o, Bond::single()).unwrap();
        
        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.num_bonds(), 2);
        assert!(mol.is_connected());
    }

    #[test]
    fn test_smiles_parsing() {
        let mol = crate::smiles_decode::decode_smiles("CCO").unwrap();
        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.num_bonds(), 2);
    }

    #[test]
    fn test_benzene() {
        let mol = crate::smiles_decode::decode_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.num_atoms(), 6);
        assert_eq!(mol.num_bonds(), 6);
        
        let rings = mol.find_rings();
        assert!(!rings.is_empty());
    }
} 
