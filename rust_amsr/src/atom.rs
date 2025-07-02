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
    // SMILES doesn't need neighbor counting - that's AMSR-specific
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
        };

        atom.parse_symbol()?;

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
        }

        // Add extra pi bonds from : symbols
        Ok(())
    }

    // SMILES doesn't need valence calculations - that's AMSR-specific

    // SMILES doesn't need these AMSR-specific methods
    // Neighbor counting and valence checking are handled differently in SMILES

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
        result.push_str(&self.atom_symbol);

        // Add charge
        if self.charge > 0 {
            result.push_str(&PLUS.repeat(self.charge as usize));
        } else if self.charge < 0 {
            result.push_str(&MINUS.repeat((-self.charge) as usize));
        }

        // Add radical electrons
        result.push_str(&RADICAL.repeat(self.radical_electrons as usize));

        result
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
    }

    #[test]
    fn test_charged_atom() {
        let atom = Atom::new("O+").unwrap();
        assert_eq!(atom.atom_symbol, "O");
        assert_eq!(atom.charge, 1);

        let atom = Atom::new("N-").unwrap();
        assert_eq!(atom.atom_symbol, "N");
        assert_eq!(atom.charge, -1);
    }

    #[test]
    fn test_aromatic_atom() {
        let atom = Atom::new("c").unwrap();
        assert_eq!(atom.atom_symbol, "C");
    }

    #[test]
    fn test_radical_atom() {
        let atom = Atom::new("C*").unwrap();
        assert_eq!(atom.atom_symbol, "C");
        assert_eq!(atom.radical_electrons, 1);
    }

    #[test]
    fn test_hypervalent_atom() {
        let atom = Atom::new("S!").unwrap();
        assert_eq!(atom.atom_symbol, "S");
        assert_eq!(atom.bangs, 1);
    }

    #[test]
    fn test_isotope_atom() {
        let atom = Atom::new("[13C]").unwrap();
        assert_eq!(atom.atom_symbol, "C");
        assert_eq!(atom.isotope, Some(13));
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
    }

    #[test]
    fn test_complex_atom() {
        let atom = Atom::new("[13s+*:!']").unwrap();
        assert_eq!(atom.atom_symbol, "S");
        assert_eq!(atom.isotope, Some(13));
        assert_eq!(atom.charge, 1);
        assert_eq!(atom.radical_electrons, 1);
        assert_eq!(atom.bangs, 1);
        assert_eq!(atom.chiral_type, ChiralType::CW);
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
