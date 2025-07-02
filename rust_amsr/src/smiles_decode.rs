use crate::errors::AMSRResult;
use crate::molecule::Molecule;
use crate::atom::Atom;
use crate::bond::{Bond, BondStereo};
use std::collections::HashMap;

/// Decode a SMILES string into a Molecule graph.
pub fn decode_smiles(smiles: &str) -> AMSRResult<Molecule> {
    let mut mol = Molecule::new();
    let mut atom_stack: Vec<usize> = Vec::new();
    let mut branch_points: Vec<usize> = Vec::new();
    let mut ring_closures: HashMap<u32, usize> = HashMap::new();
    let mut chars = smiles.chars().peekable();
    let mut prev_atom: Option<usize> = None;
    let mut pending_bond: Option<Bond> = None;

    while let Some(c) = chars.next() {
        match c {
            // Branch start
            '(' => {
                if let Some(idx) = prev_atom {
                    branch_points.push(idx);
                }
            }
            // Branch end
            ')' => {
                if let Some(idx) = branch_points.pop() {
                    prev_atom = Some(idx);
                }
            }
            // Bond symbols
            '-' => pending_bond = Some(Bond::single()),
            '=' => pending_bond = Some(Bond::double()),
            '#' => pending_bond = Some(Bond::triple()),
            ':' => pending_bond = Some(Bond::aromatic()),
            '/' => pending_bond = Some(Bond::single().with_stereo(BondStereo::E)),
            '\\' => pending_bond = Some(Bond::single().with_stereo(BondStereo::Z)),
            // Ring closure (single digit)
            d if d.is_ascii_digit() => {
                let ring_num = d.to_digit(10).unwrap();
                if let Some(&other_idx) = ring_closures.get(&ring_num) {
                    if let Some(cur_idx) = prev_atom {
                        let bond = pending_bond.take().unwrap_or_else(Bond::single);
                        mol.add_bond(cur_idx, other_idx, bond)?;
                    }
                    ring_closures.remove(&ring_num);
                } else if let Some(cur_idx) = prev_atom {
                    ring_closures.insert(ring_num, cur_idx);
                }
            }
            // Atom (including aromatic)
            c if c.is_ascii_alphabetic() => {
                let mut symbol = c.to_string();
                // Handle two-letter elements (e.g., Cl, Br)
                if let Some(&next) = chars.peek() {
                    if next.is_ascii_lowercase() {
                        symbol.push(chars.next().unwrap());
                    }
                }
                let atom = Atom::new(&symbol)?;
                let idx = mol.add_atom(atom);
                if let Some(prev) = prev_atom {
                    let bond = pending_bond.take().unwrap_or_else(Bond::single);
                    mol.add_bond(prev, idx, bond)?;
                }
                prev_atom = Some(idx);
                atom_stack.push(idx);
            }
            // Ignore unsupported characters for now
            _ => {}
        }
    }
    Ok(mol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear() {
        let mol = decode_smiles("CCO").unwrap();
        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.num_bonds(), 2);
    }

    #[test]
    fn test_ring() {
        let mol = decode_smiles("C1CCCCC1").unwrap();
        assert_eq!(mol.num_atoms(), 6);
        assert_eq!(mol.num_bonds(), 6);
    }

    #[test]
    fn test_branch() {
        let mol = decode_smiles("CC(C)C").unwrap();
        assert_eq!(mol.num_atoms(), 4);
        assert_eq!(mol.num_bonds(), 3);
    }
    
    #[test]
    fn test_aromatic() {
        let mol = decode_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.num_atoms(), 6);
        assert_eq!(mol.num_bonds(), 6);
    }
} 
