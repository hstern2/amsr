use std::process::Command;
use std::fs;

/// Convert a SMILES string to InChI using Open Babel
pub fn smiles_to_inchi(smiles: &str) -> Result<String, Box<dyn std::error::Error>> {
    let output = Command::new("obabel")
        .args(&["-:", smiles, "-oinchi", "-xFixedH"])
        .output()?;

    if !output.status.success() {
        return Err(format!("Open Babel failed: {}", String::from_utf8_lossy(&output.stderr)).into());
    }

    let output_str = String::from_utf8(output.stdout)?;
    let lines: Vec<&str> = output_str.lines().collect();

    // Find the line that starts with "InChI="
    let inchi = lines.iter()
        .find(|line| line.starts_with("InChI="))
        .ok_or("No InChI found in output")?
        .trim()
        .to_string();

    Ok(inchi)
}

/// Read SMILES from a .smi file
pub fn read_smiles_file(file_path: &str) -> Result<Vec<(String, String)>, Box<dyn std::error::Error>> {
    let content = fs::read_to_string(file_path)?;
    let mut molecules = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let smiles = parts[0].to_string();
            let name = parts[1..].join(" ");
            molecules.push((smiles, name));
        }
    }

    Ok(molecules)
}

/// Test the round-trip: SMILES -> Molecule -> SMILES -> InChI
pub fn test_smiles_roundtrip(smiles: &str, name: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("Testing: {} ({})", smiles, name);

    // Step 1: Convert original SMILES to InChI
    let original_inchi = smiles_to_inchi(smiles)?;
    println!("  Original InChI: {}", original_inchi);

    // Step 2: Parse SMILES with our code
    let mol = crate::smiles_decode::decode_smiles(smiles)?;
    println!("  Parsed molecule: {} atoms, {} bonds", mol.num_atoms(), mol.num_bonds());

    // Step 3: Re-encode to SMILES with our code
    let re_encoded_smiles = crate::smiles_encode::encode_smiles(&mol)?;
    println!("  Re-encoded SMILES: {}", re_encoded_smiles);

    // Step 4: Convert re-encoded SMILES to InChI
    let re_encoded_inchi = smiles_to_inchi(&re_encoded_smiles)?;
    println!("  Re-encoded InChI: {}", re_encoded_inchi);

    // Step 5: Compare InChI strings
    if original_inchi == re_encoded_inchi {
        println!("  ✓ PASS: InChI strings match");
        Ok(())
    } else {
        println!("  ✗ FAIL: InChI strings don't match");
        println!("    Original:  {}", original_inchi);
        println!("    Re-encoded: {}", re_encoded_inchi);
        Err(format!("InChI mismatch for {}: original='{}', re-encoded='{}'",
                   name, original_inchi, re_encoded_inchi).into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_smiles_file_roundtrip() {
        let test_file = "test1.smi";

        if !Path::new(test_file).exists() {
            eprintln!("Test file {} not found, skipping test", test_file);
            return;
        }

        let molecules = read_smiles_file(test_file).expect("Failed to read test file");

        for (smiles, name) in molecules {
            if let Err(e) = test_smiles_roundtrip(&smiles, &name) {
                panic!("Round-trip test failed for {}: {}", name, e);
            }
        }
    }

    #[test]
    fn test_individual_molecules() {
        let test_cases = vec![
            ("CCO", "ethanol"),
            ("c1ccccc1", "benzene"),
            ("CC(C)C", "isobutane"),
        ];

        for (smiles, name) in test_cases {
            if let Err(e) = test_smiles_roundtrip(smiles, name) {
                panic!("Round-trip test failed for {}: {}", name, e);
            }
        }
    }
}
