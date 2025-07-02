use regex::Regex;
use lazy_static::lazy_static;

pub const DOT: &str = ".";
pub const DIHEDRALS: [&str; 12] = ["^", "^\\", "<\\", "<", "</", "_/", "_", "\\_", "\\>", ">", "/>", "/^"];
pub const Z: &str = DIHEDRALS[0];
pub const E: &str = DIHEDRALS[6];
pub const CW: &str = "'";
pub const CCW: &str = "`";
pub const PLUS: &str = "+";
pub const MINUS: &str = "-";
pub const EXTRA_PI: &str = ":";
pub const BANG: &str = "!";
pub const RADICAL: &str = "*";
pub const L_BRACKET: &str = "[";
pub const R_BRACKET: &str = "]";
pub const L_PAREN: &str = "(";
pub const R_PAREN: &str = ")";
pub const SKIP: &str = "@";
pub const MOLSEP: &str = ";";
pub const AMPERSAND: &str = "&";

lazy_static! {
    static ref DIHEDRAL_FOR_BOND_SYMBOL: std::collections::HashMap<&'static str, i32> = {
        let mut map = std::collections::HashMap::new();
        for (i, &symbol) in DIHEDRALS.iter().enumerate() {
            let angle = if i > 6 { 30 * (i as i32 - 12) } else { 30 * i as i32 };
            map.insert(symbol, angle);
        }
        map.insert(E, -180);
        map
    };

    static ref BOND_SYMBOL_FOR_DIHEDRAL: std::collections::HashMap<i32, &'static str> = {
        let mut map = std::collections::HashMap::new();
        for (k, v) in DIHEDRAL_FOR_BOND_SYMBOL.iter() {
            map.insert(*v, *k);
        }
        map
    };

    static ref AMSR_REGEX: Regex = {
        let dihedrals_escaped = DIHEDRALS.iter()
            .map(|s| regex::escape(s))
            .collect::<Vec<_>>()
            .join("|");
        
        let bond_pattern = format!(r"(?P<bond>{})", dihedrals_escaped);
        let atom_chars = format!(r"[{}{}{}{}{}{}{}]", 
            regex::escape(PLUS), regex::escape(MINUS), regex::escape(RADICAL),
            regex::escape(EXTRA_PI), regex::escape(BANG), regex::escape(CW), regex::escape(CCW));
        
        let atom_pattern = format!(
            r"(?P<atom>{}[0-9]*[A-Za-z]+{}*{}|{}[0-9A-Za-z]+{}|[A-Za-z]{}*)",
            regex::escape(L_BRACKET), atom_chars, regex::escape(R_BRACKET),
            regex::escape(L_PAREN), regex::escape(R_PAREN), atom_chars
        );
        
        let ring_pattern = format!(
            r"(?P<ring>({}[0-9]+{}|[3-9]){}*)",
            regex::escape(L_BRACKET), regex::escape(R_BRACKET), regex::escape(SKIP)
        );
        
        let saturate_pattern = format!(r"(?P<saturate>{})", regex::escape(DOT));
        let molsep_pattern = format!(r"(?P<molsep>{})", regex::escape(MOLSEP));
        let ampersand_pattern = format!(r"(?P<ampersand>{})", regex::escape(AMPERSAND));
        
        let full_pattern = format!(
            r"({}?({}|({})))|{}|{}|{}",
            bond_pattern, atom_pattern, ring_pattern, saturate_pattern, molsep_pattern, ampersand_pattern
        );
        
        Regex::new(&full_pattern).unwrap()
    };
}

#[derive(Debug, Clone, PartialEq)]
pub enum Token {
    Bond(String),
    Atom(String),
    Ring(String),
    Saturate,
    MolSep,
    Ampersand,
}

pub fn to_tokens(s: &str) -> Vec<Token> {
    let mut tokens = Vec::new();
    
    for cap in AMSR_REGEX.captures_iter(s) {
        if let Some(bond) = cap.name("bond") {
            tokens.push(Token::Bond(bond.as_str().to_string()));
        }
        
        if let Some(atom) = cap.name("atom") {
            tokens.push(Token::Atom(atom.as_str().to_string()));
        }
        
        if let Some(ring) = cap.name("ring") {
            let ring_str = ring.as_str();
            let ring_number = ring_str.chars()
                .filter(|c| c.is_ascii_digit())
                .collect::<String>();
            tokens.push(Token::Ring(ring_number));
            
            // Add skip tokens
            let skip_count = ring_str.chars().filter(|c| *c == SKIP.chars().next().unwrap()).count();
            for _ in 0..skip_count {
                tokens.push(Token::Atom(SKIP.to_string()));
            }
        }
        
        if cap.name("saturate").is_some() {
            tokens.push(Token::Saturate);
        }
        
        if cap.name("molsep").is_some() {
            tokens.push(Token::MolSep);
        }
        
        if cap.name("ampersand").is_some() {
            tokens.push(Token::Ampersand);
        }
    }
    
    tokens
}

pub fn get_dihedral_angle(symbol: &str) -> Option<i32> {
    DIHEDRAL_FOR_BOND_SYMBOL.get(symbol).copied()
}

pub fn get_bond_symbol(angle: i32) -> Option<&'static str> {
    BOND_SYMBOL_FOR_DIHEDRAL.get(&angle).copied()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_token_parsing() {
        let amsr = "CCO";
        let tokens = to_tokens(amsr);
        assert_eq!(tokens, vec![
            Token::Atom("C".to_string()),
            Token::Atom("C".to_string()),
            Token::Atom("O".to_string()),
        ]);
    }

    #[test]
    fn test_dihedral_angles() {
        assert_eq!(get_dihedral_angle("^"), Some(0));
        assert_eq!(get_dihedral_angle("_"), Some(180));
        assert_eq!(get_dihedral_angle("E"), Some(-180));
    }
} 