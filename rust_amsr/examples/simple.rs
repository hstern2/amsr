use amsr::{Molecule, smiles_to_amsr, amsr_to_smiles, molecule_to_amsr, amsr_to_molecule};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("AMSR Rust Library Example");
    println!("========================\n");

    // Example 1: Simple molecule
    let smiles = "CCO";
    println!("Example 1: Ethanol");
    println!("SMILES: {}", smiles);

    let amsr = smiles_to_amsr(smiles)?;
    println!("AMSR: {}", amsr);

    let back_to_smiles = amsr_to_smiles(&amsr)?;
    println!("Back to SMILES: {}", back_to_smiles);
    println!();

    // Example 2: Benzene
    let smiles = "c1ccccc1";
    println!("Example 2: Benzene");
    println!("SMILES: {}", smiles);

    let amsr = smiles_to_amsr(smiles)?;
    println!("AMSR: {}", amsr);

    let back_to_smiles = amsr_to_smiles(&amsr)?;
    println!("Back to SMILES: {}", back_to_smiles);
    println!();

    // Example 3: Working with Molecule objects
    println!("Example 3: Working with Molecule objects");
    let mol = amsr_to_molecule(&amsr)?;
    println!("Number of atoms: {}", mol.num_atoms());
    println!("Number of bonds: {}", mol.num_bonds());
    println!("Is connected: {}", mol.is_connected());

    let rings = mol.find_rings();
    println!("Number of rings: {}", rings.len());
    println!();

    // Example 4: Creating molecule from scratch
    println!("Example 4: Creating molecule from scratch");
    let mut mol = Molecule::new();

    // Add atoms
    let c1 = mol.add_atom(amsr::atom::Atom::new("C")?);
    let c2 = mol.add_atom(amsr::atom::Atom::new("C")?);
    let o = mol.add_atom(amsr::atom::Atom::new("O")?);

    // Add bonds
    mol.add_bond(c1, c2, amsr::bond::Bond::single())?;
    mol.add_bond(c2, o, amsr::bond::Bond::single())?;

    let amsr = molecule_to_amsr(&mol)?;
    println!("AMSR from scratch: {}", amsr);

    let smiles = amsr::molecule_to_smiles(&mol)?;
    println!("SMILES from scratch: {}", smiles);

    Ok(())
}
