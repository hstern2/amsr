use crate::molecule::Molecule;
use std::collections::{HashMap, HashSet};

pub struct SMILES_Encoder;

#[derive(Debug)]
enum SMILES_Symbol {
    Atom(String),
    Bond(String),
    BranchStart,
    BranchEnd,
    RingClosure(usize),
}

#[derive(Debug, Clone)]
struct RingBond {
    atom1: usize,
    atom2: usize,
    bond: Option<String>,
}

impl SMILES_Encoder {
    pub fn new() -> Self {
        SMILES_Encoder
    }

    pub fn encode(&self, mol: &Molecule) -> Result<String, Box<dyn std::error::Error>> {
        if mol.atoms.is_empty() {
            return Ok(String::new());
        }

        let mut result = String::new();
        let mut visited = HashSet::new();

        // Find connected components
        let components = self.find_connected_components(mol);

        for (i, component) in components.iter().enumerate() {
            if i > 0 {
                result.push('.');
            }

            if let Some(&start_atom) = component.first() {
                let mut symbols = Vec::new();
                let mut ring_bonds = Vec::new(); // Track ring closure bonds
                let mut atom_positions = HashMap::new(); // Track where each atom was first encountered

                self.encode_dfs(
                    mol,
                    start_atom,
                    None,
                    &mut visited,
                    &mut symbols,
                    &mut ring_bonds,
                    &mut atom_positions,
                )?;

                // Assign ring numbers to ring bonds
                let mut ring_numbers = HashMap::new();
                let mut next_ring_number = 1;

                // Deduplicate ring bonds
                let mut unique_ring_bonds = HashSet::new();
                for ring_bond in &ring_bonds {
                    let key = if ring_bond.atom1 < ring_bond.atom2 {
                        (ring_bond.atom1, ring_bond.atom2)
                    } else {
                        (ring_bond.atom2, ring_bond.atom1)
                    };
                    unique_ring_bonds.insert(key);
                }

                println!("    Detected {} unique ring bonds:", unique_ring_bonds.len());
                for &(atom1, atom2) in &unique_ring_bonds {
                    println!("      Bond: {} - {}", atom1, atom2);
                }

                for &(atom1, atom2) in &unique_ring_bonds {
                    if !ring_numbers.contains_key(&(atom1, atom2)) {
                        ring_numbers.insert((atom1, atom2), next_ring_number);
                        println!("      Assigned ring number {} to bond {} - {}", next_ring_number, atom1, atom2);
                        next_ring_number += 1;
                    }
                }

                // Insert ring closures at the correct positions
                // For each ring bond, insert ring closure after the first atom and at the end
                for &(atom1, atom2) in &unique_ring_bonds {
                    let ring_num = ring_numbers[&(atom1, atom2)];
                    let pos1 = atom_positions.get(&atom1).unwrap_or(&0);
                    let pos2 = atom_positions.get(&atom2).unwrap_or(&0);

                    // Insert ring closure after the first atom (smaller position)
                    let first_pos = pos1.min(pos2);
                    println!("    Inserting ring closure {} after atom at position {}", ring_num, first_pos);
                    symbols.insert(first_pos + 1, SMILES_Symbol::RingClosure(ring_num));

                    // Insert bond type and ring closure at the end
                    if let Some(bond) = mol.get_bond(atom1, atom2) {
                        symbols.push(SMILES_Symbol::Bond(bond.to_smiles()));
                    }
                    println!("    Inserting bond and ring closure {} at the end", ring_num);
                    symbols.push(SMILES_Symbol::RingClosure(ring_num));
                }

                // Build the final SMILES string from symbols
                let smiles_string = self.build_smiles_string(&symbols);
                println!("    Final SMILES symbols: {:?}", symbols);
                println!("    Final SMILES string: {}", smiles_string);
                result.push_str(&smiles_string);
            }
        }

        Ok(result)
    }

    fn build_smiles_string(&self, symbols: &[SMILES_Symbol]) -> String {
        let mut result = String::new();
        for symbol in symbols {
            match symbol {
                SMILES_Symbol::Atom(s) => result.push_str(s),
                SMILES_Symbol::Bond(s) => result.push_str(s),
                SMILES_Symbol::BranchStart => result.push('('),
                SMILES_Symbol::BranchEnd => result.push(')'),
                SMILES_Symbol::RingClosure(num) => result.push_str(&num.to_string()),
            }
        }
        result
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

    fn encode_dfs(
        &self,
        mol: &Molecule,
        current: usize,
        parent: Option<usize>,
        visited: &mut HashSet<usize>,
        symbols: &mut Vec<SMILES_Symbol>,
        ring_bonds: &mut Vec<RingBond>,
        atom_positions: &mut HashMap<usize, usize>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        visited.insert(current);

        // Record the position where this atom was first encountered
        atom_positions.insert(current, symbols.len());

        // Add atom
        symbols.push(SMILES_Symbol::Atom(mol.atoms[current].to_smiles()));

        // Collect neighbors and sort them for consistent output
        let neighbors = mol.get_neighbors(current);
        let mut neighbors: Vec<_> = neighbors.iter().collect();
        neighbors.sort();

        // Process neighbors
        for &&neighbor in &neighbors {
            if Some(neighbor) == parent {
                continue;
            }

            if !visited.contains(&neighbor) {
                // New atom - add bond and recurse
                if let Some(bond) = mol.get_bond(current, neighbor) {
                    symbols.push(SMILES_Symbol::Bond(bond.to_smiles()));
                }

                // Check if this neighbor has multiple unvisited neighbors (branching)
                let neighbor_neighbors = mol.get_neighbors(neighbor);
                let unvisited_neighbors: Vec<_> = neighbor_neighbors
                    .iter()
                    .filter(|&&n| n != current && !visited.contains(&n))
                    .collect();

                if unvisited_neighbors.len() > 1 {
                    symbols.push(SMILES_Symbol::BranchStart);
                }

                self.encode_dfs(
                    mol,
                    neighbor,
                    Some(current),
                    visited,
                    symbols,
                    ring_bonds,
                    atom_positions,
                )?;

                if unvisited_neighbors.len() > 1 {
                    symbols.push(SMILES_Symbol::BranchEnd);
                }
            } else {
                // Ring closure - track this bond for later processing
                let bond_symbol = mol.get_bond(current, neighbor)
                    .map(|b| b.to_smiles())
                    .unwrap_or_else(|| "".to_string());

                ring_bonds.push(RingBond {
                    atom1: current,
                    atom2: neighbor,
                    bond: Some(bond_symbol),
                });
            }
        }

        Ok(())
    }
}

impl Default for SMILES_Encoder {
    fn default() -> Self {
        SMILES_Encoder::new()
    }
}

pub fn encode_smiles(mol: &Molecule) -> Result<String, Box<dyn std::error::Error>> {
    SMILES_Encoder::new().encode(mol)
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

    #[test]
    fn test_branching_encoding() {
        let mol = crate::smiles_decode::decode_smiles("CC(C)C").unwrap();
        let smiles = encode_smiles(&mol).unwrap();
        assert!(smiles.contains("(") && smiles.contains(")"));
    }
}
