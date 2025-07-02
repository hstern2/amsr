use crate::molecule::Molecule;
use crate::atom::Atom;
use crate::bond::{Bond, BondStereo};
use crate::aromaticity::kekulize_molecule;
use std::collections::{HashMap, HashSet};

/// Represents a pending ring closure with optional bond type
#[derive(Debug, Clone)]
struct RingClosure {
    atom_idx: usize,
    bond: Option<Bond>,
}

/// Decode a SMILES string into a Molecule graph.
pub fn decode_smiles(smiles: &str) -> Result<Molecule, Box<dyn std::error::Error>> {
    let mut mol = Molecule::new();
    let mut atom_stack: Vec<usize> = Vec::new();
    let mut branch_points: Vec<usize> = Vec::new();
    let mut ring_closures: HashMap<u32, RingClosure> = HashMap::new();
    let mut chars = smiles.chars().peekable();
    let mut prev_atom: Option<usize> = None;
    let mut pending_bond: Option<Bond> = None;
    let mut aromatic_atoms = HashSet::new();

    while let Some(c) = chars.next() {
        match c {
            // Branch star
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
                handle_ring_closure(&mut mol, &mut ring_closures, &mut prev_atom,
                                  &mut pending_bond, ring_num)?;
            }
            // Multi-digit ring closure (e.g., %10, %11)
            '%' => {
                if let Some(ring_num) = parse_ring_number(&mut chars) {
                    handle_ring_closure(&mut mol, &mut ring_closures, &mut prev_atom,
                                      &mut pending_bond, ring_num)?;
                }
            }
            // Atom (including aromatic)
            c if c.is_ascii_alphabetic() => {
                let mut symbol = c.to_string();
                // Handle two-letter elements (e.g., Cl, Br)
                if let Some(&next) = chars.peek() {
                    if next.is_ascii_lowercase() && !c.is_lowercase() {
                        // This is a two-letter element (e.g., Cl, Br)
                        symbol.push(chars.next().unwrap());
                    }
                }
                let atom = Atom::new(&symbol)?;
                let idx = mol.add_atom(atom);
                if c.is_lowercase() {
                    aromatic_atoms.insert(idx);
                }
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

    // Kekulize aromatic systems if present
    kekulize_molecule(&mut mol, &aromatic_atoms)?;

    Ok(mol)
}

/// Handle ring closure logic for both single and multi-digit ring numbers
fn handle_ring_closure(
    mol: &mut Molecule,
    ring_closures: &mut HashMap<u32, RingClosure>,
    prev_atom: &mut Option<usize>,
    pending_bond: &mut Option<Bond>,
    ring_num: u32,
) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(ring_closure) = ring_closures.get(&ring_num) {
        // This is the second occurrence of this ring number - create the bond
        if let Some(cur_idx) = *prev_atom {
            let bond = pending_bond.take().unwrap_or_else(|| {
                // Use the bond from the first occurrence if available, otherwise single bond
                ring_closure.bond.clone().unwrap_or_else(Bond::single)
            });
            mol.add_bond(cur_idx, ring_closure.atom_idx, bond)?;
        }
        ring_closures.remove(&ring_num);
    } else if let Some(cur_idx) = *prev_atom {
        // This is the first occurrence of this ring number - store it
        let bond = pending_bond.take();
        ring_closures.insert(ring_num, RingClosure {
            atom_idx: cur_idx,
            bond,
        });
    }
    // If prev_atom is None, this is an invalid SMILES (ring closure before atom)
    // We ignore it as per SMILES specification
    Ok(())
}

/// Parse multi-digit ring numbers (e.g., %10, %11, etc.)
fn parse_ring_number(chars: &mut std::iter::Peekable<std::str::Chars<'_>>) -> Option<u32> {
    if chars.next()? == '%' {
        let mut num_str = String::new();
        while let Some(&c) = chars.peek() {
            if c.is_ascii_digit() {
                num_str.push(chars.next().unwrap());
            } else {
                break;
            }
        }
        num_str.parse::<u32>().ok()
    } else {
        None
    }
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

    #[test]
    fn test_multi_digit_ring() {
        let mol = decode_smiles("C1CCCCC%101").unwrap();
        assert_eq!(mol.num_atoms(), 6);
        assert_eq!(mol.num_bonds(), 6);
    }

    #[test]
    fn test_ring_with_explicit_bonds() {
        let mol = decode_smiles("C1=CC=CC=C1").unwrap();
        assert_eq!(mol.num_atoms(), 6);
        assert_eq!(mol.num_bonds(), 6);
    }

    #[test]
    fn test_complex_ring_system() {
        let mol = decode_smiles("C1CC2CCCC2CC1").unwrap();
        assert_eq!(mol.num_atoms(), 8);
        assert_eq!(mol.num_bonds(), 9);
    }

    #[test]
    fn test_ring_closure_before_atom() {
        // This tests the case where a ring closure appears before an atom
        // This is technically invalid SMILES but should be handled gracefully
        // The leading '1' should be ignored, so we get CCCCC1 = 5 atoms
        let mol = decode_smiles("1CCCCC1").unwrap();
        assert_eq!(mol.num_atoms(), 5);
        assert_eq!(mol.num_bonds(), 4); // 4 bonds between 5 atoms in a chain
    }
}
