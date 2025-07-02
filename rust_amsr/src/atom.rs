use crate::valence::get_valence;
use crate::tokens::{PLUS, MINUS, RADICAL, EXTRA_PI, BANG, CW, CCW};
use regex::Regex;
use lazy_static::lazy_static;

lazy_static! {
    static ref ISOTOPE_REGEX: Regex = Regex::new(r"\[(\d+)").unwrap();
}

#[derive(Debug, Clone, PartialEq)]
pub enum ChiralType {
    Unspecified,
    CW,
    CCW,
}

#[derive(Debug, Clone)]
pub struct Atom {
    pub symbol: String,
    pub atom_symbol: String,
    pub isotope: Option<u32>,
    pub charge: i32,
    pub radical_electrons: i32,
    pub bangs: i32,
    pub chiral_type: ChiralType,
    pub max_pi_bonds: i32,
    pub n_pi_bonds: i32,
    pub max_neighbors: i32,
    pub n_neighbors: i32,
    pub is_saturated: bool,
}

impl Atom {
    pub fn new(symbol: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let mut atom = Atom {
            symbol: symbol.to_string(),
            atom_symbol: String::new(),
            isotope: None,
            charge: 0,
            radical_electrons: 0,
            bangs: 0,
            chiral_type: ChiralType::Unspecified,
            max_pi_bonds: 0,
            n_pi_bonds: 0,
            max_neighbors: 0,
            n_neighbors: 0,
            is_saturated: false,
        };

        atom.parse_symbol()?;
        atom.calculate_properties()?;

        Ok(atom)
    }

    fn parse_symbol(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        // Parse isotope (e.g., [13C] -> 13)
        if let Some(captures) = ISOTOPE_REGEX.captures(&self.symbol) {
            if let Ok(isotope) = captures[1].parse::<u32>() {
                self.isotope = Some(isotope);
            }
        }

        // Parse charge (+ and - symbols)
        self.charge = self.symbol.matches(PLUS).count() as i32 - self.symbol.matches(MINUS).count() as i32;

        // Parse radical electrons (* symbols)
        self.radical_electrons = self.symbol.matches(RADICAL).count() as i32;

        // Parse bangs (! symbols for hypervalent atoms)
        self.bangs = self.symbol.matches(BANG).count() as i32;
        println!("Bangs for symbol '{}': {}", self.symbol, self.bangs);

        // Parse chirality (' and ` symbols)
        if self.symbol.contains(CCW) {
            self.chiral_type = ChiralType::CCW;
        } else if self.symbol.contains(CW) {
            self.chiral_type = ChiralType::CW;
        }

        // Extract atom symbol (letters only, removing all special characters including numbers)
        self.atom_symbol = self.symbol.chars()
            .filter(|c| c.is_ascii_alphabetic())
            .collect::<String>();

        // Debug prin
        println!("Symbol: '{}', Atom symbol: '{}'", self.symbol, self.atom_symbol);

        // Handle aromatic atoms (lowercase first letter)
        if !self.atom_symbol.is_empty() && self.atom_symbol.chars().next().unwrap().is_lowercase() {
            // Convert first letter to uppercase and set aromatic pi bonds
            let first_char = self.atom_symbol.chars().next().unwrap().to_uppercase().next().unwrap();
            self.atom_symbol = first_char.to_string() + &self.atom_symbol[1..];
            self.max_pi_bonds = 1;
        }

        // Add extra pi bonds from : symbols
        self.max_pi_bonds += 2 * self.symbol.matches(EXTRA_PI).count() as i32;

        Ok(())
    }

    fn calculate_properties(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        // Get valence from lookup table
        println!("Valence lookup: symbol='{}', charge={}, bangs={}", self.atom_symbol, self.charge, self.bangs);
        let valence = get_valence(&self.atom_symbol, self.charge, self.bangs)
            .ok_or_else(|| format!("Unknown valence for {}{}{}",
                self.atom_symbol, self.charge, self.bangs))?;

        // Calculate max neighbors: valence - radical electrons - pi bonds
        self.max_neighbors = valence - self.radical_electrons - self.max_pi_bonds;

        if self.max_neighbors < 0 {
            return Err(format!("Invalid valence state for atom {}", self.symbol).into());
        }

        Ok(())
    }

    pub fn add_bond_to(&mut self, other: &mut Atom) {
        self.n_neighbors += 1;
        other.n_neighbors += 1;
    }

    pub fn can_bond(&self) -> bool {
        !self.is_saturated && self.n_neighbors < self.max_neighbors
    }

    pub fn can_bond_with(&self, other: &Atom, stringent: bool) -> bool {
        if !stringent {
            return true;
        }

        // Oxygen-oxygen bonds are unstable
        if self.is_oxygen() && other.is_oxygen() {
            return false;
        }

        // Halogen-heteroatom bonds are often unstable
        if self.is_halogen() && other.is_hetero() {
            return false;
        }

        true
    }

    pub fn n_available_pi_bonds(&self) -> i32 {
        self.max_pi_bonds - self.n_pi_bonds
    }

    pub fn is_carbon(&self) -> bool {
        self.atom_symbol == "C"
    }

    pub fn is_sulfur(&self) -> bool {
        self.atom_symbol == "S"
    }

    pub fn is_oxygen(&self) -> bool {
        self.atom_symbol == "O"
    }

    pub fn is_hetero(&self) -> bool {
        !self.is_carbon()
    }

    pub fn is_halogen(&self) -> bool {
        matches!(self.atom_symbol.as_str(), "F" | "Cl" | "Br" | "I" | "At" | "Ts")
    }

    pub fn is_hnh(&self) -> bool {
        self.is_hetero() && !self.is_halogen()
    }

    pub fn is_aromatic(&self) -> bool {
        // If symbol starts with '[', look for the first alphabetic character after any isotope
        if self.symbol.starts_with('[') {
            let mut chars = self.symbol[1..].chars();
            // Skip digits (isotope)
            while let Some(c) = chars.next() {
                if c.is_ascii_alphabetic() {
                    return c.is_lowercase();
                }
            }
            return false;
        }
        // Otherwise, check if the first alphabetic character is lowercase
        for c in self.symbol.chars() {
            if c.is_ascii_alphabetic() {
                return c.is_lowercase();
            }
        }
        false
    }

    pub fn symbol_with(&self, suffix: &str) -> String {
        if self.symbol.starts_with('[') {
            format!("[{}{}{}]", &self.symbol[1..self.symbol.len()-1], suffix, "]")
        } else {
            format!("{}{}", self.symbol, suffix)
        }
    }

    pub fn to_smiles(&self) -> String {
        let mut result = String::new();

        // Add isotope if presen
        if let Some(isotope) = self.isotope {
            result.push_str(&isotope.to_string());
        }

        // Add atom symbol
        if self.is_aromatic() {
            result.push_str(&self.atom_symbol.to_lowercase());
        } else {
            result.push_str(&self.atom_symbol);
        }

        // Add charge
        if self.charge > 0 {
            result.push_str(&PLUS.repeat(self.charge as usize));
        } else if self.charge < 0 {
            result.push_str(&MINUS.repeat((-self.charge) as usize));
        }

        // Add radical electrons
        result.push_str(&RADICAL.repeat(self.radical_electrons as usize));

        resul
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_carbon_atom() {
        let atom = Atom::new("C").unwrap();
        assert_eq!(atom.atom_symbol, "C");
        assert_eq!(atom.charge, 0);
        assert_eq!(atom.radical_electrons, 0);
        assert_eq!(atom.bangs, 0);
        assert_eq!(atom.max_neighbors, 4);
        assert_eq!(atom.max_pi_bonds, 0);
    }

    #[test]
    fn test_charged_atom() {
        let atom = Atom::new("O+").unwrap();
        assert_eq!(atom.atom_symbol, "O");
        assert_eq!(atom.charge, 1);
        assert_eq!(atom.max_neighbors, 3); // valence(3) - radicals(0) - pi_bonds(0) = 3

        let atom = Atom::new("N-").unwrap();
        assert_eq!(atom.atom_symbol, "N");
        assert_eq!(atom.charge, -1);
        assert_eq!(atom.max_neighbors, 2); // valence(2) - radicals(0) - pi_bonds(0) = 2
    }

    #[test]
    fn test_aromatic_atom() {
        let atom = Atom::new("c").unwrap();
        assert_eq!(atom.atom_symbol, "C");
        assert_eq!(atom.max_pi_bonds, 1);
        assert_eq!(atom.max_neighbors, 3); // valence(4) - radicals(0) - pi_bonds(1) = 3
        assert!(atom.is_aromatic());
    }

    #[test]
    fn test_radical_atom() {
        let atom = Atom::new("C*").unwrap();
        assert_eq!(atom.atom_symbol, "C");
        assert_eq!(atom.radical_electrons, 1);
        assert_eq!(atom.max_neighbors, 3); // valence(4) - radicals(1) - pi_bonds(0) = 3
    }

    #[test]
    fn test_hypervalent_atom() {
        let atom = Atom::new("S!").unwrap();
        assert_eq!(atom.atom_symbol, "S");
        assert_eq!(atom.bangs, 1);
        assert_eq!(atom.max_neighbors, 4); // valence(4) - radicals(0) - pi_bonds(0) = 4
    }

    #[test]
    fn test_isotope_atom() {
        let atom = Atom::new("[13C]").unwrap();
        assert_eq!(atom.atom_symbol, "C");
        assert_eq!(atom.isotope, Some(13));
        assert_eq!(atom.max_neighbors, 4);
    }

    #[test]
    fn test_chiral_atom() {
        let atom = Atom::new("C'").unwrap();
        assert_eq!(atom.atom_symbol, "C");
        assert_eq!(atom.chiral_type, ChiralType::CW);

        let atom = Atom::new("C`").unwrap();
        assert_eq!(atom.atom_symbol, "C");
        assert_eq!(atom.chiral_type, ChiralType::CCW);
    }

    #[test]
    fn test_extra_pi_bonds() {
        let atom = Atom::new("C::").unwrap();
        assert_eq!(atom.atom_symbol, "C");
        assert_eq!(atom.max_pi_bonds, 4); // 2 * 2 colons = 4
        assert_eq!(atom.max_neighbors, 0); // valence(4) - radicals(0) - pi_bonds(4) = 0
    }

    #[test]
    fn test_complex_atom() {
        let atom = Atom::new("[13s+*:!']").unwrap();
        assert_eq!(atom.atom_symbol, "S");
        assert_eq!(atom.isotope, Some(13));
        assert_eq!(atom.charge, 1);
        assert_eq!(atom.radical_electrons, 1);
        assert_eq!(atom.bangs, 1);
        assert_eq!(atom.max_pi_bonds, 3); // aromatic(1) + extra(2) = 3
        assert_eq!(atom.chiral_type, ChiralType::CW);
        assert!(atom.is_aromatic());
    }

    #[test]
    fn test_halogen_atoms() {
        assert!(Atom::new("F").unwrap().is_halogen());
        assert!(Atom::new("Cl").unwrap().is_halogen());
        assert!(Atom::new("Br").unwrap().is_halogen());
        assert!(Atom::new("I").unwrap().is_halogen());
        assert!(!Atom::new("C").unwrap().is_halogen());
    }

    #[test]
    fn test_heteroatom_detection() {
        assert!(!Atom::new("C").unwrap().is_hetero());
        assert!(Atom::new("N").unwrap().is_hetero());
        assert!(Atom::new("O").unwrap().is_hetero());
        assert!(Atom::new("S").unwrap().is_hetero());
    }
}
