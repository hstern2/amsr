use crate::errors::AMSRResult;
use crate::molecule::Molecule;
use std::collections::{HashMap, HashSet};

pub struct SMILESEncoder {
    canonical: bool,
    isomeric: bool,
    kekule: bool,
}

impl SMILESEncoder {
    pub fn new() -> Self {
        SMILESEncoder {
            canonical: false,
            isomeric: true,
            kekule: false,
        }
    }
    
    pub fn canonical(mut self, canonical: bool) -> Self {
        self.canonical = canonical;
        self
    }
    
    pub fn isomeric(mut self, isomeric: bool) -> Self {
        self.isomeric = isomeric;
        self
    }
    
    pub fn kekule(mut self, kekule: bool) -> Self {
        self.kekule = kekule;
        self
    }
    
    pub fn encode(&self, mol: &Molecule) -> AMSRResult<String> {
        if mol.atoms.is_empty() {
            return Ok(String::new());
        }
        
        let mut result = String::new();
        let mut visited = HashSet::new();
        let mut ring_numbers = HashMap::new();
        let mut next_ring_number = 1;
        
        // Find connected components
        let components = self.find_connected_components(mol);
        
        for (i, component) in components.iter().enumerate() {
            if i > 0 {
                result.push('.');
            }
            
            if let Some(&start_atom) = component.first() {
                self.encode_dfs(mol, start_atom, &mut visited, &mut ring_numbers, 
                              &mut next_ring_number, &mut result)?;
            }
        }
        
        Ok(result)
    }
    
    fn find_connected_components(&self, mol: &Molecule) -> Vec<Vec<usize>> {
        let mut components = Vec::new();
        let mut visited = HashSet::new();
        
        for start in 0..mol.atoms.len() {
            if !visited.contains(&start) {
                let mut component = Vec::new();
                let mut stack = Vec::new();
                stack.push(start);
                visited.insert(start);
                
                while let Some(current) = stack.pop() {
                    component.push(current);
                    for &neighbor in &mol.get_neighbors(current) {
                        if !visited.contains(&neighbor) {
                            visited.insert(neighbor);
                            stack.push(neighbor);
                        }
                    }
                }
                
                components.push(component);
            }
        }
        
        components
    }
    
    fn encode_dfs(&self, mol: &Molecule, current: usize, visited: &mut HashSet<usize>,
                  ring_numbers: &mut HashMap<(usize, usize), usize>, 
                  next_ring_number: &mut usize, result: &mut String) -> AMSRResult<()> {
        visited.insert(current);
        
        // Add atom
        result.push_str(&mol.atoms[current].to_smiles());
        
        // Process neighbors
        let neighbors = mol.get_neighbors(current);
        for &neighbor in &neighbors {
            if !visited.contains(&neighbor) {
                // New atom - add bond and recurse
                if let Some(bond) = mol.get_bond(current, neighbor) {
                    result.push_str(&bond.to_smiles());
                }
                self.encode_dfs(mol, neighbor, visited, ring_numbers, next_ring_number, result)?;
            } else {
                // Ring closure
                let bond_key = if current < neighbor { (current, neighbor) } else { (neighbor, current) };
                if !ring_numbers.contains_key(&bond_key) {
                    ring_numbers.insert(bond_key, *next_ring_number);
                    *next_ring_number += 1;
                }
                result.push_str(&ring_numbers[&bond_key].to_string());
            }
        }
        
        Ok(())
    }
}

impl Default for SMILESEncoder {
    fn default() -> Self {
        SMILESEncoder::new()
    }
}

pub fn encode_smiles(mol: &Molecule) -> AMSRResult<String> {
    SMILESEncoder::new().encode(mol)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::Bond;

    #[test]
    fn test_simple_encoding() {
        let mut mol = Molecule::new();
        let c1 = mol.add_atom(Atom::new("C").unwrap());
        let c2 = mol.add_atom(Atom::new("C").unwrap());
        let o = mol.add_atom(Atom::new("O").unwrap());
        
        mol.add_bond(c1, c2, Bond::single()).unwrap();
        mol.add_bond(c2, o, Bond::single()).unwrap();
        
        let smiles = encode_smiles(&mol).unwrap();
        assert_eq!(smiles, "CCO");
    }

    #[test]
    fn test_ring_encoding() {
        let mol = crate::smiles_decode::decode_smiles("C1CCCCC1").unwrap();
        let smiles = encode_smiles(&mol).unwrap();
        assert!(smiles.contains("1"));
    }
} 
