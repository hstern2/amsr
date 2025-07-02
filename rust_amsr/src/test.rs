// SMILES round-trip test runner
// Usage: cargo run --bin test -- <smiles_file.smi>
// Each line in the file should be: SMILES name

use std::process::Command;
use std::fs;
use std::io::Write;
use std::env;

/// Check if Open Babel is available
fn check_obabel() -> Result<(), Box<dyn std::error::Error>> {
    let output = Command::new("obabel").arg("--version").output();
    match output {
        Ok(_) => Ok(()),
        Err(_) => Err("Open Babel (obabel) not found. Please install Open Babel and ensure it's in your PATH.".into())
    }
}

/// Convert a SMILES string to InChI using Open Babel
pub fn smiles_to_inchi(smiles: &str) -> Result<String, Box<dyn std::error::Error>> {
    let output = Command::new("obabel")
        .arg(format!("-:{}", smiles))
        .args(&["-oinchi", "-xFixedH"])
        .output()?;

    let output_str = String::from_utf8_lossy(&output.stdout);
    let stderr_str = String::from_utf8_lossy(&output.stderr);
    let lines: Vec<&str> = output_str.lines().collect();

    if let Some(inchi_line) = lines.iter().find(|line| line.starts_with("InChI=")) {
        Ok(inchi_line.trim().to_string())
    } else {
        eprintln!("Open Babel failed for SMILES '{}':\nSTDOUT:\n{}\nSTDERR:\n{}", smiles, output_str, stderr_str);
        Err(format!("No InChI found in Open Babel output for SMILES: {}", smiles).into())
    }
}

/// Read SMILES from a .smi file
/// Expected format: one line per molecule, "SMILES name" (whitespace-separated)
/// Lines starting with # are treated as comments and ignored
/// Empty lines are ignored
pub fn read_smiles_file(file_path: &str) -> Result<Vec<(String, String)>, Box<dyn std::error::Error>> {
    let content = fs::read_to_string(file_path)?;
    let mut molecules = Vec::new();

    for (line_num, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let smiles = parts[0].to_string();
            let name = parts[1..].join(" ");
            molecules.push((smiles, name));
        } else if parts.len() == 1 {
            molecules.push((parts[0].to_string(), String::from("(no_name)")));
        } else {
            eprintln!("Warning: Skipping malformed line {}: '{}'", line_num + 1, line);
        }
    }

    Ok(molecules)
}

/// Test the round-trip: SMILES -> Molecule -> SMILES -> InChI
pub fn test_smiles_roundtrip(smiles: &str, name: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("Testing: {} ({})", smiles, name);

    // Step 1: Convert original SMILES to InChI
    let original_inchi = match smiles_to_inchi(smiles) {
        Ok(inchi) => inchi,
        Err(e) => return Err(format!("Failed to convert original SMILES to InChI: {}", e).into())
    };
    println!("  Original InChI: {}", original_inchi);

    // Step 2: Parse SMILES with our code
    let mol = match crate::smiles_decode::decode_smiles(smiles) {
        Ok(mol) => mol,
        Err(e) => return Err(format!("Failed to decode SMILES: {}", e).into())
    };
    println!("  Parsed molecule: {} atoms, {} bonds", mol.num_atoms(), mol.num_bonds());

    // Debug: Print adjacency list to verify structure
    println!("  Adjacency list:");
    for i in 0..mol.num_atoms() {
        let neighbors = mol.get_neighbors(i);
        println!("    Atom {}: connected to {:?}", i, neighbors);
    }

    // Debug: Print bond types
    println!("  Bond types:");
    for i in 0..mol.num_atoms() {
        for &neighbor in &mol.get_neighbors(i) {
            if i < neighbor {
                if let Some(bond) = mol.get_bond(i, neighbor) {
                    println!("    Bond {} - {}: {:?}", i, neighbor, bond);
                }
            }
        }
    }

    // Step 3: Re-encode to SMILES with our code
    let re_encoded_smiles = match crate::smiles_encode::encode_smiles(&mol) {
        Ok(smiles) => smiles,
        Err(e) => return Err(format!("Failed to re-encode SMILES: {}", e).into())
    };
    println!("  Re-encoded SMILES: {}", re_encoded_smiles);

    // Step 4: Convert re-encoded SMILES to InChI
    let re_encoded_inchi = match smiles_to_inchi(&re_encoded_smiles) {
        Ok(inchi) => inchi,
        Err(e) => return Err(format!("Failed to convert re-encoded SMILES to InChI: {}", e).into())
    };
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

/// Run round-trip tests on a SMILES file
fn run_smiles_file_tests(file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let molecules = read_smiles_file(file_path)?;
    println!("Testing {} molecules from '{}'", molecules.len(), file_path);
    println!("=");

    let mut n_pass = 0;
    let mut n_fail = 0;
    let mut failures = Vec::new();

    for (smiles, name) in molecules {
        match test_smiles_roundtrip(&smiles, &name) {
            Ok(_) => n_pass += 1,
            Err(e) => {
                let failure_msg = format!("FAIL: {} ({}) - {}", smiles, name, e);
                eprintln!("  {}", failure_msg);
                failures.push(failure_msg);
                n_fail += 1;
            }
        }
        println!("---");
    }

    println!("=");
    println!("Summary: {} passed, {} failed", n_pass, n_fail);

    // Write failures to a summary file
    if !failures.is_empty() {
        let summary_file = format!("test_failures_{}.txt", file_path.replace(".smi", "").replace("/", "_"));
        if let Ok(mut file) = fs::File::create(&summary_file) {
            writeln!(file, "SMILES Round-trip Test Failures").unwrap();
            writeln!(file, "Generated from: {}", file_path).unwrap();
            writeln!(file, "=").unwrap();
            for failure in &failures {
                writeln!(file, "{}", failure).unwrap();
            }
            println!("Failure details written to: {}", summary_file);
        } else {
            eprintln!("Warning: Could not write failure summary to {}", summary_file);
        }
    }

    if n_fail > 0 {
        return Err(format!("{} tests failed in {}", n_fail, file_path).into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_basic_smiles() {
        // Check if Open Babel is available
        if let Err(e) = check_obabel() {
            eprintln!("Skipping test: {}", e);
            return;
        }

        let test_file = "test1.smi";
        if !Path::new(test_file).exists() {
            eprintln!("Test file {} not found, skipping test", test_file);
            return;
        }

        run_smiles_file_tests(test_file).expect("Basic SMILES tests failed");
    }

    #[test]
    fn test_complex_smiles() {
        // Check if Open Babel is available
        if let Err(e) = check_obabel() {
            eprintln!("Skipping test: {}", e);
            return;
        }

        let test_file = "test_complex.smi";
        if !Path::new(test_file).exists() {
            eprintln!("Test file {} not found, skipping test", test_file);
            return;
        }

        run_smiles_file_tests(test_file).expect("Complex SMILES tests failed");
    }

    #[test]
    fn test_invalid_smiles() {
        // Check if Open Babel is available
        if let Err(e) = check_obabel() {
            eprintln!("Skipping test: {}", e);
            return;
        }

        let test_file = "test_invalid.smi";
        if !Path::new(test_file).exists() {
            eprintln!("Test file {} not found, skipping test", test_file);
            return;
        }

        // For invalid SMILES, we expect some failures, so we don't panic on errors
        if let Err(e) = run_smiles_file_tests(test_file) {
            println!("Invalid SMILES test completed with expected errors: {}", e);
        }
    }
}
