pub mod atom;
pub mod bond;
pub mod molecule;
pub mod tokens;
pub mod valence;
pub mod errors;
pub mod amsr_encode;
pub mod amsr_decode;
pub mod smiles_encode;
pub mod smiles_decode;

use errors::AMSRResult;
pub use amsr_encode::{encode_molecule, encode_molecule_to_tokens, AMSREncoder};
pub use amsr_decode::{decode_amsr, decode_amsr_with_dihedrals, AMSRDecoder};
pub use smiles_encode::{encode_smiles, SMILESEncoder};
pub use smiles_decode::decode_smiles;

/// Convert a SMILES string to AMSR
pub fn smiles_to_amsr(smiles: &str) -> AMSRResult<String> {
    let mol = decode_smiles(smiles)?;
    encode_molecule(&mol)
}

/// Convert AMSR to a SMILES string
pub fn amsr_to_smiles(amsr: &str) -> AMSRResult<String> {
    let mol = decode_amsr(amsr)?;
    encode_smiles(&mol)
}

/// Convert a molecule graph to AMSR
pub fn molecule_to_amsr(mol: &molecule::Molecule) -> AMSRResult<String> {
    encode_molecule(mol)
}

/// Convert AMSR to a molecule graph
pub fn amsr_to_molecule(amsr: &str) -> AMSRResult<molecule::Molecule> {
    decode_amsr(amsr)
}

/// Convert a SMILES string to a molecule graph
pub fn smiles_to_molecule(smiles: &str) -> AMSRResult<molecule::Molecule> {
    decode_smiles(smiles)
}

/// Convert a molecule graph to SMILES
pub fn molecule_to_smiles(mol: &molecule::Molecule) -> AMSRResult<String> {
    encode_smiles(mol)
}

pub use molecule::Molecule;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_molecule() {
        let smiles = "CCO";
        let amsr = smiles_to_amsr(smiles).unwrap();
        let back_to_smiles = amsr_to_smiles(&amsr).unwrap();
        assert_eq!(smiles, back_to_smiles);
    }

    #[test]
    fn test_benzene() {
        let smiles = "c1ccccc1";
        let amsr = smiles_to_amsr(smiles).unwrap();
        let back_to_smiles = amsr_to_smiles(&amsr).unwrap();
        assert_eq!(smiles, back_to_smiles);
    }
} 