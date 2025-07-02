use std::collections::HashMap;
use lazy_static::lazy_static;

lazy_static! {
    static ref VALENCE: HashMap<(String, i32, i32), i32> = {
        let mut map = HashMap::new();
        
        // Hydrogen
        map.insert(("H".to_string(), 0, 0), 1);
        map.insert(("H".to_string(), 1, 0), 0);
        map.insert(("H".to_string(), -1, 0), 0);
        
        // Helium
        map.insert(("He".to_string(), 0, 0), 0);
        
        // Lithium
        map.insert(("Li".to_string(), 1, 0), 0);
        map.insert(("Li".to_string(), 0, 0), 1);
        
        // Beryllium
        map.insert(("Be".to_string(), 2, 0), 0);
        
        // Boron
        map.insert(("B".to_string(), 0, 0), 3);
        map.insert(("B".to_string(), -1, 0), 4);
        
        // Carbon
        map.insert(("C".to_string(), 0, 0), 4);
        map.insert(("C".to_string(), -1, 0), 3);
        map.insert(("C".to_string(), 1, 0), 3);
        
        // Nitrogen
        map.insert(("N".to_string(), 0, 0), 3);
        map.insert(("N".to_string(), 1, 0), 4);
        map.insert(("N".to_string(), -1, 0), 2);
        
        // Oxygen
        map.insert(("O".to_string(), 0, 0), 2);
        map.insert(("O".to_string(), -1, 0), 1);
        map.insert(("O".to_string(), -2, 0), 0);
        map.insert(("O".to_string(), 1, 0), 3);
        
        // Fluorine
        map.insert(("F".to_string(), 0, 0), 1);
        map.insert(("F".to_string(), -1, 0), 0);
        map.insert(("F".to_string(), 1, 0), 2);
        
        // Sodium
        map.insert(("Na".to_string(), 0, 0), 1);
        map.insert(("Na".to_string(), 1, 0), 0);
        
        // Magnesium
        map.insert(("Mg".to_string(), 0, 0), 2);
        map.insert(("Mg".to_string(), 0, 1), 6);
        map.insert(("Mg".to_string(), 1, 0), 1);
        map.insert(("Mg".to_string(), 2, 0), 0);
        
        // Aluminum
        map.insert(("Al".to_string(), 0, 0), 3);
        map.insert(("Al".to_string(), 3, 0), 0);
        map.insert(("Al".to_string(), -3, 0), 6);
        
        // Xenon
        map.insert(("Xe".to_string(), 0, 0), 0);
        
        // Silicon
        map.insert(("Si".to_string(), 0, 0), 4);
        map.insert(("Si".to_string(), -1, 0), 5);
        map.insert(("Si".to_string(), 4, 0), 0);
        
        // Phosphorus
        map.insert(("P".to_string(), 0, 0), 3);
        map.insert(("P".to_string(), 0, 1), 5);
        map.insert(("P".to_string(), 0, 2), 7);
        map.insert(("P".to_string(), 1, 0), 4);
        map.insert(("P".to_string(), -1, 1), 6);
        
        // Sulfur
        map.insert(("S".to_string(), 0, 0), 2);
        map.insert(("S".to_string(), 0, 1), 4);
        map.insert(("S".to_string(), 0, 2), 6);
        map.insert(("S".to_string(), 1, 0), 3);
        map.insert(("S".to_string(), 1, 1), 5);
        map.insert(("S".to_string(), -1, 0), 1);
        map.insert(("S".to_string(), -1, 1), 3);
        map.insert(("S".to_string(), -1, 2), 5);
        map.insert(("S".to_string(), -2, 0), 0);
        
        // Chlorine
        map.insert(("Cl".to_string(), 0, 0), 1);
        map.insert(("Cl".to_string(), -1, 0), 0);
        map.insert(("Cl".to_string(), 1, 0), 2);
        map.insert(("Cl".to_string(), 2, 0), 3);
        map.insert(("Cl".to_string(), 3, 0), 4);
        
        // Potassium
        map.insert(("K".to_string(), 0, 0), 1);
        map.insert(("K".to_string(), 1, 0), 0);
        
        // Calcium
        map.insert(("Ca".to_string(), 0, 0), 2);
        map.insert(("Ca".to_string(), 2, 0), 0);
        
        // Iron
        map.insert(("Fe".to_string(), 3, 0), 0);
        
        // Zinc
        map.insert(("Zn".to_string(), 0, 0), 2);
        map.insert(("Zn".to_string(), 1, 0), 1);
        map.insert(("Zn".to_string(), 2, 0), 1);
        map.insert(("Zn".to_string(), -2, 0), 2);
        map.insert(("Zn".to_string(), -2, 1), 4);
        
        // Arsenic
        map.insert(("As".to_string(), 0, 0), 3);
        map.insert(("As".to_string(), 0, 1), 5);
        map.insert(("As".to_string(), 0, 2), 7);
        map.insert(("As".to_string(), 1, 0), 4);
        map.insert(("As".to_string(), -1, 1), 6);
        
        // Selenium
        map.insert(("Se".to_string(), 0, 0), 2);
        map.insert(("Se".to_string(), 0, 1), 4);
        map.insert(("Se".to_string(), 0, 2), 6);
        map.insert(("Se".to_string(), 1, 0), 3);
        map.insert(("Se".to_string(), 1, 1), 5);
        map.insert(("Se".to_string(), -1, 0), 1);
        map.insert(("Se".to_string(), -2, 0), 0);
        
        // Cesium
        map.insert(("Cs".to_string(), 1, 0), 0);
        map.insert(("Cs".to_string(), 0, 0), 1);
        
        // Barium
        map.insert(("Ba".to_string(), 0, 0), 2);
        map.insert(("Ba".to_string(), 2, 0), 0);
        
        // Bismuth
        map.insert(("Bi".to_string(), 0, 0), 3);
        map.insert(("Bi".to_string(), 3, 0), 0);
        
        // Bromine
        map.insert(("Br".to_string(), 0, 0), 1);
        map.insert(("Br".to_string(), -1, 0), 0);
        map.insert(("Br".to_string(), 2, 0), 3);
        
        // Zirconium
        map.insert(("Zr".to_string(), 4, 0), 0);
        
        // Krypton
        map.insert(("Kr".to_string(), 0, 0), 0);
        
        // Rubidium
        map.insert(("Rb".to_string(), 0, 0), 1);
        map.insert(("Rb".to_string(), 1, 0), 0);
        
        // Strontium
        map.insert(("Sr".to_string(), 0, 0), 2);
        map.insert(("Sr".to_string(), 2, 0), 0);
        
        // Silver
        map.insert(("Ag".to_string(), 0, 0), 2);
        map.insert(("Ag".to_string(), 1, 0), 0);
        map.insert(("Ag".to_string(), -4, 0), 3);
        
        // Tellurium
        map.insert(("Te".to_string(), 0, 0), 2);
        map.insert(("Te".to_string(), 1, 0), 3);
        map.insert(("Te".to_string(), 0, 1), 4);
        map.insert(("Te".to_string(), 0, 2), 6);
        map.insert(("Te".to_string(), -1, 1), 3);
        map.insert(("Te".to_string(), -1, 2), 5);
        
        // Iodine
        map.insert(("I".to_string(), 0, 0), 1);
        map.insert(("I".to_string(), 0, 1), 3);
        map.insert(("I".to_string(), 0, 2), 5);
        map.insert(("I".to_string(), -1, 0), 0);
        map.insert(("I".to_string(), 1, 0), 2);
        map.insert(("I".to_string(), 2, 1), 3);
        map.insert(("I".to_string(), 3, 0), 4);
        
        // Astatine
        map.insert(("At".to_string(), 0, 0), 1);
        
        // Gold
        map.insert(("Au".to_string(), 0, 0), 0);
        
        // Radium
        map.insert(("Ra".to_string(), 0, 0), 2);
        map.insert(("Ra".to_string(), 2, 0), 0);
        
        map
    };
    
    static ref BANGS: HashMap<(String, i32, i32), i32> = {
        let mut map = HashMap::new();
        
        // Generate BANGS table from VALENCE table
        for ((sym, chg, bangs), val) in VALENCE.iter() {
            if *bangs > 0 {
                map.insert((sym.clone(), *chg, *val), *bangs);
            }
        }
        
        map
    };
}

pub fn get_valence(atom_symbol: &str, charge: i32, bangs: i32) -> Option<i32> {
    VALENCE.get(&(atom_symbol.to_string(), charge, bangs)).copied()
}

pub fn get_bangs(atom_symbol: &str, charge: i32, total_valence: i32) -> i32 {
    BANGS.get(&(atom_symbol.to_string(), charge, total_valence)).copied().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_carbon_valence() {
        assert_eq!(get_valence("C", 0, 0), Some(4));
        assert_eq!(get_valence("C", 1, 0), Some(3));
        assert_eq!(get_valence("C", -1, 0), Some(3));
    }

    #[test]
    fn test_oxygen_valence() {
        assert_eq!(get_valence("O", 0, 0), Some(2));
        assert_eq!(get_valence("O", 1, 0), Some(3));
        assert_eq!(get_valence("O", -1, 0), Some(1));
        assert_eq!(get_valence("O", -2, 0), Some(0));
    }

    #[test]
    fn test_sulfur_bangs() {
        assert_eq!(get_bangs("S", 0, 4), 1);
        assert_eq!(get_bangs("S", 0, 6), 2);
        assert_eq!(get_bangs("S", 1, 5), 1);
    }

    #[test]
    fn test_phosphorus_bangs() {
        assert_eq!(get_bangs("P", 0, 5), 1);
        assert_eq!(get_bangs("P", 0, 7), 2);
        assert_eq!(get_bangs("P", -1, 6), 1);
    }
} 
