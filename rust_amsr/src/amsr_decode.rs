use crate::errors::{AMSRResult, AMSRError};
use crate::molecule::Molecule;
use crate::atom::Atom;
use crate::bond::Bond;
use crate::tokens::{to_tokens, Token};
use std::collections::HashMap;

pub struct AMSRDecoder {
    stringent: bool,
}

impl AMSRDecoder {
    pub fn new() -> Self {
        AMSRDecoder {
            stringent: true,
        }
    }
    
    pub fn stringent(mut self, stringent: bool) -> Self {
        self.stringent = stringent;
        self
    }
    
    pub fn decode(&self, amsr: &str) -> AMSRResult<Molecule> {
        let mut mol = Molecule::new();
        let tokens = to_tokens(amsr);
        
        let mut atom_stack = Vec::new();
        let mut bond_stack = Vec::new();
        let mut ring_map = HashMap::new();
        let _next_ring_number = 1;
        
        for token in tokens {
            match token {
                Token::Atom(symbol) => {
                    let atom = Atom::new(&symbol)?;
                    let atom_idx = mol.add_atom(atom);
                    atom_stack.push(atom_idx);
                },
                Token::Bond(symbol) => {
                    let bond = Bond::from_symbol(&symbol)?;
                    bond_stack.push(bond);
                },
                Token::Ring(ring_str) => {
                    if let Ok(ring_num) = ring_str.parse::<usize>() {
                        if ring_num < 3 {
                            return Err(AMSRError::InvalidAMSR(format!("Invalid ring size: {}", ring_num)));
                        }
                        
                        // Handle ring closure
                        if let Some(atom_idx) = atom_stack.last().copied() {
                            if let Some(target_idx) = ring_map.get(&ring_num) {
                                // Create ring bond
                                let bond = bond_stack.pop().unwrap_or_else(Bond::single);
                                mol.add_bond(atom_idx, *target_idx, bond)?;
                            } else {
                                // First occurrence of this ring number
                                ring_map.insert(ring_num, atom_idx);
                            }
                        }
                    }
                },
                Token::Saturate => {
                    // Mark atom as saturated
                    if let Some(atom_idx) = atom_stack.last() {
                        mol.atoms[*atom_idx].is_saturated = true;
                    }
                },
                Token::MolSep => {
                    // Start new molecule component
                    atom_stack.clear();
                    bond_stack.clear();
                    ring_map.clear();
                },
                Token::Ampersand => {
                    // Handle ampersand (not implemented in this simplified version)
                },
            }
        }
        
        // Process any remaining bonds
        while let (Some(atom1), Some(atom2), Some(bond)) = (
            atom_stack.pop(),
            atom_stack.pop(),
            bond_stack.pop()
        ) {
            mol.add_bond(atom1, atom2, bond)?;
        }
        
        // Saturate remaining atoms
        for atom in &mut mol.atoms {
            if atom.can_bond() {
                atom.is_saturated = true;
            }
        }
        
        Ok(mol)
    }
    
    pub fn decode_with_dihedrals(&self, amsr: &str) -> AMSRResult<(Molecule, HashMap<(usize, usize, usize, usize), i32>)> {
        let mut mol = Molecule::new();
        let mut dihedrals = HashMap::new();
        let tokens = to_tokens(amsr);
        
        let mut atom_stack = Vec::new();
        let mut bond_stack = Vec::new();
        let mut ring_map = HashMap::new();
        let mut bond_dihedral_map = HashMap::new();
        
        for token in tokens {
            match token {
                Token::Atom(symbol) => {
                    let atom = Atom::new(&symbol)?;
                    let atom_idx = mol.add_atom(atom);
                    atom_stack.push(atom_idx);
                },
                Token::Bond(symbol) => {
                    let bond = Bond::from_symbol(&symbol)?;
                    bond_stack.push(bond);
                },
                Token::Ring(ring_str) => {
                    if let Ok(ring_num) = ring_str.parse::<usize>() {
                        if ring_num < 3 {
                            return Err(AMSRError::InvalidAMSR(format!("Invalid ring size: {}", ring_num)));
                        }
                        
                        // Handle ring closure
                        if let Some(atom_idx) = atom_stack.last().copied() {
                            if let Some(target_idx) = ring_map.get(&ring_num) {
                                // Create ring bond
                                let bond = bond_stack.pop().unwrap_or_else(Bond::single);
                                let dihedral_angle = bond.dihedral_angle;
                                let bond_idx = mol.bonds.len();
                                mol.add_bond(atom_idx, *target_idx, bond)?;
                                
                                // Store dihedral information
                                if let Some(angle) = dihedral_angle {
                                    bond_dihedral_map.insert(bond_idx, angle);
                                }
                            } else {
                                // First occurrence of this ring number
                                ring_map.insert(ring_num, atom_idx);
                            }
                        }
                    }
                },
                Token::Saturate => {
                    // Mark atom as saturated
                    if let Some(atom_idx) = atom_stack.last() {
                        mol.atoms[*atom_idx].is_saturated = true;
                    }
                },
                Token::MolSep => {
                    // Start new molecule component
                    atom_stack.clear();
                    bond_stack.clear();
                    ring_map.clear();
                },
                Token::Ampersand => {
                    // Handle ampersand (not implemented in this simplified version)
                },
            }
        }
        
        // Process any remaining bonds
        while let (Some(atom1), Some(atom2), Some(bond)) = (
            atom_stack.pop(),
            atom_stack.pop(),
            bond_stack.pop()
        ) {
            let dihedral_angle = bond.dihedral_angle;
            let bond_idx = mol.bonds.len();
            mol.add_bond(atom1, atom2, bond)?;
            
            // Store dihedral information
            if let Some(angle) = dihedral_angle {
                bond_dihedral_map.insert(bond_idx, angle);
            }
        }
        
        // Convert bond dihedral map to atom quartet dihedral map
        for (bond_idx, angle) in bond_dihedral_map {
            if let Some((atom1, atom2, _)) = mol.bonds.get(bond_idx) {
                // Find neighboring atoms for dihedral definition
                if let (Some(neighbor1), Some(neighbor2)) = (
                    mol.get_neighbors(*atom1).iter().find(|&&n| n != *atom2),
                    mol.get_neighbors(*atom2).iter().find(|&&n| n != *atom1)
                ) {
                    dihedrals.insert((*neighbor1, *atom1, *atom2, *neighbor2), angle);
                }
            }
        }
        
        // Saturate remaining atoms
        for atom in &mut mol.atoms {
            if atom.can_bond() {
                atom.is_saturated = true;
            }
        }
        
        Ok((mol, dihedrals))
    }
}

impl Default for AMSRDecoder {
    fn default() -> Self {
        AMSRDecoder::new()
    }
}

pub fn decode_amsr(amsr: &str) -> AMSRResult<Molecule> {
    AMSRDecoder::new().decode(amsr)
}

pub fn decode_amsr_with_dihedrals(amsr: &str) -> AMSRResult<(Molecule, HashMap<(usize, usize, usize, usize), i32>)> {
    AMSRDecoder::new().decode_with_dihedrals(amsr)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_decoding() {
        let amsr = "CCO";
        let mol = decode_amsr(amsr).unwrap();
        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.num_bonds(), 2);
    }

    #[test]
    fn test_ring_decoding() {
        let amsr = "c1ccccc1";
        let mol = decode_amsr(amsr).unwrap();
        assert_eq!(mol.num_atoms(), 6);
        assert_eq!(mol.num_bonds(), 6);
    }

    #[test]
    fn test_invalid_ring_size() {
        let amsr = "c2ccccc2"; // Invalid ring size 2
        let result = decode_amsr(amsr);
        assert!(result.is_err());
    }
} 