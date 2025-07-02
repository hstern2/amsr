use crate::errors::AMSRResult;
use crate::molecule::Molecule;
use crate::tokens::{DOT, MOLSEP};
use std::collections::{HashSet, VecDeque};

pub struct AMSREncoder {
    use_groups: bool,
    stringent: bool,
    randomize: bool,
    canonical: bool,
    use_stereo: bool,
}

impl AMSREncoder {
    pub fn new() -> Self {
        AMSREncoder {
            use_groups: true,
            stringent: true,
            randomize: false,
            canonical: false,
            use_stereo: true,
        }
    }
    
    pub fn use_groups(mut self, use_groups: bool) -> Self {
        self.use_groups = use_groups;
        self
    }
    
    pub fn stringent(mut self, stringent: bool) -> Self {
        self.stringent = stringent;
        self
    }
    
    pub fn randomize(mut self, randomize: bool) -> Self {
        self.randomize = randomize;
        self
    }
    
    pub fn canonical(mut self, canonical: bool) -> Self {
        self.canonical = canonical;
        self
    }
    
    pub fn use_stereo(mut self, use_stereo: bool) -> Self {
        self.use_stereo = use_stereo;
        self
    }
    
    pub fn encode(&self, mol: &Molecule) -> AMSRResult<String> {
        let tokens = self.encode_to_tokens(mol)?;
        Ok(tokens.into_iter().collect::<String>())
    }
    
    pub fn encode_to_tokens(&self, mol: &Molecule) -> AMSRResult<Vec<String>> {
        if mol.atoms.is_empty() {
            return Ok(Vec::new());
        }
        
        let mut tokens = Vec::new();
        let mut visited = HashSet::new();
        let mut seen_bonds = HashSet::new();
        let mut n_seen_atoms = 0;
        
        // Find connected components
        let components = self.find_connected_components(mol);
        
        for component in components {
            if !tokens.is_empty() {
                tokens.push(MOLSEP.to_string());
            }
            
            // Start with the first atom in the component
            if let Some(&start_atom) = component.first() {
                tokens.push(mol.atoms[start_atom].symbol.clone());
                self.encode_dfs(mol, start_atom, &mut visited, &mut seen_bonds, 
                              &mut n_seen_atoms, &mut tokens)?;
            }
        }
        
        // Remove trailing dots
        while tokens.last() == Some(&DOT.to_string()) {
            tokens.pop();
        }
        
        Ok(tokens)
    }
    
    fn find_connected_components(&self, mol: &Molecule) -> Vec<Vec<usize>> {
        let mut components = Vec::new();
        let mut visited = HashSet::new();
        
        for start in 0..mol.atoms.len() {
            if !visited.contains(&start) {
                let mut component = Vec::new();
                let mut queue = VecDeque::new();
                queue.push_back(start);
                visited.insert(start);
                
                while let Some(current) = queue.pop_front() {
                    component.push(current);
                    for &neighbor in &mol.get_neighbors(current) {
                        if !visited.contains(&neighbor) {
                            visited.insert(neighbor);
                            queue.push_back(neighbor);
                        }
                    }
                }
                
                components.push(component);
            }
        }
        
        components
    }
    
    fn encode_dfs(&self, mol: &Molecule, current: usize, visited: &mut HashSet<usize>,
                  seen_bonds: &mut HashSet<(usize, usize)>, n_seen_atoms: &mut usize,
                  tokens: &mut Vec<String>) -> AMSRResult<()> {
        visited.insert(current);
        *n_seen_atoms += 1;
        
        // Get neighbors sorted by search order
        let mut neighbors = mol.get_neighbors(current);
        neighbors.sort_by(|&a, &b| self.search_order(mol, current, a, b));
        
        for &neighbor in &neighbors {
            let bond_key = if current < neighbor { (current, neighbor) } else { (neighbor, current) };
            
            if seen_bonds.contains(&bond_key) {
                continue;
            }
            
            seen_bonds.insert(bond_key);
            
            if visited.contains(&neighbor) {
                // Ring closure
                if let Some(bond) = mol.get_bond(current, neighbor) {
                    if let Some(ref symbol) = bond.symbol {
                        tokens.push(symbol.clone());
                    }
                }
                
                // Find ring size
                let ring_size = self.find_ring_size(mol, current, neighbor, visited);
                if ring_size >= 3 {
                    tokens.push(ring_size.to_string());
                }
            } else {
                // New atom
                if let Some(bond) = mol.get_bond(current, neighbor) {
                    if let Some(ref symbol) = bond.symbol {
                        tokens.push(symbol.clone());
                    }
                }
                
                tokens.push(mol.atoms[neighbor].symbol.clone());
                self.encode_dfs(mol, neighbor, visited, seen_bonds, n_seen_atoms, tokens)?;
            }
        }
        
        // Add saturation dot if atom can still bond
        if mol.atoms[current].can_bond() {
            tokens.push(DOT.to_string());
        }
        
        Ok(())
    }
    
    fn search_order(&self, mol: &Molecule, current: usize, a: usize, b: usize) -> std::cmp::Ordering {
        // 1. Seen atoms before unseen atoms (i.e. rings)
        let a_seen = mol.atoms[a].n_neighbors > 0;
        let b_seen = mol.atoms[b].n_neighbors > 0;
        
        if a_seen != b_seen {
            return a_seen.cmp(&b_seen);
        }
        
        // 2. Aromatic bonds first (simplified)
        let a_aromatic = mol.atoms[a].is_aromatic();
        let b_aromatic = mol.atoms[b].is_aromatic();
        
        if a_aromatic != b_aromatic {
            return a_aromatic.cmp(&b_aromatic).reverse();
        }
        
        // 3. Small rings before larger (for seen) .. otherwise atom index (for unseen)
        if a_seen {
            // For seen atoms, prefer smaller ring sizes
            let a_ring_size = self.estimate_ring_size(mol, current, a);
            let b_ring_size = self.estimate_ring_size(mol, current, b);
            a_ring_size.cmp(&b_ring_size)
        } else {
            // For unseen atoms, use atom index
            a.cmp(&b)
        }
    }
    
    fn estimate_ring_size(&self, mol: &Molecule, current: usize, target: usize) -> usize {
        // Simplified ring size estimation
        // In a full implementation, this would use BFS to find the shortest path
        let mut queue = VecDeque::new();
        let mut visited = HashSet::new();
        queue.push_back((current, 0));
        visited.insert(current);
        
        while let Some((atom, distance)) = queue.pop_front() {
            if atom == target {
                return distance;
            }
            
            for &neighbor in &mol.get_neighbors(atom) {
                if !visited.contains(&neighbor) {
                    visited.insert(neighbor);
                    queue.push_back((neighbor, distance + 1));
                }
            }
        }
        
        usize::MAX
    }
    
    fn find_ring_size(&self, mol: &Molecule, start: usize, end: usize, visited: &HashSet<usize>) -> usize {
        // Simplified ring size calculation
        // In a full implementation, this would use BFS to find the shortest path
        let mut queue = VecDeque::new();
        let mut local_visited = HashSet::new();
        queue.push_back((start, 0));
        local_visited.insert(start);
        
        while let Some((atom, distance)) = queue.pop_front() {
            if atom == end {
                return distance + 1; // +1 for the bond back to start
            }
            
            for &neighbor in &mol.get_neighbors(atom) {
                if !local_visited.contains(&neighbor) && visited.contains(&neighbor) {
                    local_visited.insert(neighbor);
                    queue.push_back((neighbor, distance + 1));
                }
            }
        }
        
        0
    }
}

impl Default for AMSREncoder {
    fn default() -> Self {
        AMSREncoder::new()
    }
}

pub fn encode_molecule(mol: &Molecule) -> AMSRResult<String> {
    AMSREncoder::new().encode(mol)
}

pub fn encode_molecule_to_tokens(mol: &Molecule) -> AMSRResult<Vec<String>> {
    AMSREncoder::new().encode_to_tokens(mol)
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
        
        let amsr = encode_molecule(&mol).unwrap();
        assert!(!amsr.is_empty());
    }

    #[test]
    fn test_ring_encoding() {
        let mol = crate::smiles_decode::decode_smiles("c1ccccc1").unwrap();
        let amsr = encode_molecule(&mol).unwrap();
        assert!(!amsr.is_empty());
    }
} 
