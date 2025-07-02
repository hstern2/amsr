use crate::errors::{AMSRResult, AMSRError};
use crate::tokens::{E, Z, get_dihedral_angle, get_bond_symbol};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BondType {
    Single,
    Double,
    Triple,
    Aromatic,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BondStereo {
    None,
    E,
    Z,
    Cis,
    Trans,
}

#[derive(Debug, Clone)]
pub struct Bond {
    pub bond_type: BondType,
    pub stereo: BondStereo,
    pub symbol: Option<String>,
    pub is_rotatable: bool,
    pub dihedral_angle: Option<i32>,
}

impl Bond {
    pub fn new() -> Self {
        Bond {
            bond_type: BondType::Single,
            stereo: BondStereo::None,
            symbol: None,
            is_rotatable: false,
            dihedral_angle: None,
        }
    }
    
    pub fn from_symbol(symbol: &str) -> AMSRResult<Self> {
        let mut bond = Bond::new();
        
        // Check if it's a dihedral symbol
        if let Some(angle) = get_dihedral_angle(symbol) {
            bond.dihedral_angle = Some(angle);
            bond.symbol = Some(symbol.to_string());
            
            // Determine if it's E/Z stereo
            if symbol == E {
                bond.stereo = BondStereo::E;
            } else if symbol == Z {
                bond.stereo = BondStereo::Z;
            }
        } else {
            return Err(AMSRError::InvalidBond(format!("Unknown bond symbol: {}", symbol)));
        }
        
        Ok(bond)
    }
    
    pub fn single() -> Self {
        Bond {
            bond_type: BondType::Single,
            stereo: BondStereo::None,
            symbol: None,
            is_rotatable: true,
            dihedral_angle: None,
        }
    }
    
    pub fn double() -> Self {
        Bond {
            bond_type: BondType::Double,
            stereo: BondStereo::None,
            symbol: None,
            is_rotatable: false,
            dihedral_angle: None,
        }
    }
    
    pub fn triple() -> Self {
        Bond {
            bond_type: BondType::Triple,
            stereo: BondStereo::None,
            symbol: None,
            is_rotatable: false,
            dihedral_angle: None,
        }
    }
    
    pub fn aromatic() -> Self {
        Bond {
            bond_type: BondType::Aromatic,
            stereo: BondStereo::None,
            symbol: None,
            is_rotatable: false,
            dihedral_angle: None,
        }
    }
    
    pub fn with_stereo(mut self, stereo: BondStereo) -> Self {
        self.stereo = stereo;
        self
    }
    
    pub fn with_dihedral(mut self, angle: i32) -> Self {
        self.dihedral_angle = Some(angle);
        if let Some(symbol) = get_bond_symbol(angle) {
            self.symbol = Some(symbol.to_string());
        }
        self
    }
    
    pub fn to_smiles(&self) -> String {
        match self.bond_type {
            BondType::Single => {
                match self.stereo {
                    BondStereo::E => "/".to_string(),
                    BondStereo::Z => "\\".to_string(),
                    _ => "".to_string(),
                }
            },
            BondType::Double => "=".to_string(),
            BondType::Triple => "#".to_string(),
            BondType::Aromatic => ":".to_string(),
        }
    }
    
    pub fn is_rotatable(&self) -> bool {
        self.is_rotatable
    }
    
    pub fn set_rotatable(&mut self, rotatable: bool) {
        self.is_rotatable = rotatable;
    }
}

impl Default for Bond {
    fn default() -> Self {
        Bond::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_bond() {
        let bond = Bond::single();
        assert_eq!(bond.bond_type, BondType::Single);
        assert!(bond.is_rotatable());
    }

    #[test]
    fn test_double_bond() {
        let bond = Bond::double();
        assert_eq!(bond.bond_type, BondType::Double);
        assert!(!bond.is_rotatable());
    }

    #[test]
    fn test_dihedral_bond() {
        let bond = Bond::from_symbol("^").unwrap();
        assert_eq!(bond.dihedral_angle, Some(0));
        assert_eq!(bond.symbol, Some("^".to_string()));
    }

    #[test]
    fn test_e_z_stereo() {
        let bond_e = Bond::from_symbol("E").unwrap();
        assert_eq!(bond_e.stereo, BondStereo::E);
        
        let bond_z = Bond::from_symbol("Z").unwrap();
        assert_eq!(bond_z.stereo, BondStereo::Z);
    }
} 