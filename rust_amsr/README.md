# AMSR Rust Library

A Rust implementation of the AMSR (Atom-Mapped SMILES Representation) molecular graph library. This library provides functionality for encoding and decoding molecular structures in both SMILES and AMSR formats, with a focus on molecular graph representation.

## Features

- **Molecular Graph Representation**: Simple and efficient molecular graph data structure
- **SMILES Support**: Parse and generate SMILES strings
- **AMSR Support**: Encode and decode AMSR representations
- **Graph Algorithms**: Ring detection, connectivity analysis, and more
- **Stereochemistry**: Basic support for E/Z stereochemistry
- **Valence Rules**: Comprehensive valence and bonding rules for common elements
- **Error Handling**: Robust error handling with custom error types

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
amsr = { path = "./rust_amsr" }
```

## Quick Start

```rust
use amsr::{smiles_to_amsr, amsr_to_smiles, Molecule};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Convert SMILES to AMSR
    let smiles = "CCO";
    let amsr = smiles_to_amsr(smiles)?;
    println!("AMSR: {}", amsr);
    
    // Convert AMSR back to SMILES
    let back_to_smiles = amsr_to_smiles(&amsr)?;
    println!("SMILES: {}", back_to_smiles);
    
    Ok(())
}
```

## Core Components

### Molecule

The main data structure representing a molecular graph:

```rust
use amsr::{Molecule, atom::Atom, bond::Bond};

let mut mol = Molecule::new();

// Add atoms
let c1 = mol.add_atom(Atom::new("C")?);
let c2 = mol.add_atom(Atom::new("C")?);
let o = mol.add_atom(Atom::new("O")?);

// Add bonds
mol.add_bond(c1, c2, Bond::single())?;
mol.add_bond(c2, o, Bond::single())?;

// Analyze the molecule
println!("Atoms: {}", mol.num_atoms());
println!("Bonds: {}", mol.num_bonds());
println!("Connected: {}", mol.is_connected());
```

### Atom

Represents an atom with its properties:

```rust
use amsr::atom::Atom;

let carbon = Atom::new("C")?;
let oxygen = Atom::new("O")?;
let charged_nitrogen = Atom::new("[NH4+]")?;

println!("Carbon valence: {}", carbon.max_neighbors);
println!("Oxygen is hetero: {}", oxygen.is_hetero());
```

### Bond

Represents chemical bonds with stereochemistry:

```rust
use amsr::bond::{Bond, BondType, BondStereo};

let single_bond = Bond::single();
let double_bond = Bond::double();
let e_bond = Bond::from_symbol("E")?; // E-stereochemistry

println!("Bond type: {:?}", single_bond.bond_type);
println!("Is rotatable: {}", single_bond.is_rotatable());
```

## API Reference

### Main Functions

- `smiles_to_amsr(smiles: &str) -> AMSRResult<String>`: Convert SMILES to AMSR
- `amsr_to_smiles(amsr: &str) -> AMSRResult<String>`: Convert AMSR to SMILES
- `molecule_to_amsr(mol: &Molecule) -> AMSRResult<String>`: Convert molecule to AMSR
- `amsr_to_molecule(amsr: &str) -> AMSRResult<Molecule>`: Convert AMSR to molecule
- `smiles_to_molecule(smiles: &str) -> AMSRResult<Molecule>`: Convert SMILES to molecule
- `molecule_to_smiles(mol: &Molecule) -> AMSRResult<String>`: Convert molecule to SMILES

### Molecule Methods

- `new()`: Create a new empty molecule
- `add_atom(atom: Atom) -> usize`: Add an atom and return its index
- `add_bond(atom1: usize, atom2: usize, bond: Bond) -> AMSRResult<()>`: Add a bond
- `get_neighbors(atom_idx: usize) -> Vec<usize>`: Get neighboring atom indices
- `get_bond(atom1: usize, atom2: usize) -> Option<&Bond>`: Get bond between atoms
- `num_atoms() -> usize`: Get number of atoms
- `num_bonds() -> usize`: Get number of bonds
- `is_connected() -> bool`: Check if molecule is connected
- `find_rings() -> Vec<Vec<usize>>`: Find all rings in the molecule
- `to_smiles() -> AMSRResult<String>`: Convert to SMILES
- `to_amsr() -> AMSRResult<String>`: Convert to AMSR

### Atom Methods

- `new(symbol: &str) -> AMSRResult<Atom>`: Create atom from symbol
- `add_bond_to(other: &mut Atom)`: Add bond to another atom
- `can_bond() -> bool`: Check if atom can form more bonds
- `can_bond_with(other: &Atom, stringent: bool) -> bool`: Check bonding compatibility
- `is_carbon() -> bool`: Check if atom is carbon
- `is_hetero() -> bool`: Check if atom is heteroatom
- `is_halogen() -> bool`: Check if atom is halogen
- `to_smiles() -> String`: Convert to SMILES representation

### Bond Methods

- `new()`: Create a new bond
- `from_symbol(symbol: &str) -> AMSRResult<Bond>`: Create bond from symbol
- `single()`, `double()`, `triple()`, `aromatic()`: Create specific bond types
- `with_stereo(stereo: BondStereo) -> Bond`: Add stereochemistry
- `with_dihedral(angle: i32) -> Bond`: Add dihedral angle
- `to_smiles() -> String`: Convert to SMILES representation

## Supported Elements

The library supports common elements with their valence rules:

- **Carbon (C)**: 4 valence electrons
- **Nitrogen (N)**: 3 valence electrons (can be charged)
- **Oxygen (O)**: 2 valence electrons (can be charged)
- **Sulfur (S)**: 2, 4, or 6 valence electrons
- **Phosphorus (P)**: 3 or 5 valence electrons
- **Halogens (F, Cl, Br, I)**: 1 valence electron
- **Hydrogen (H)**: 1 valence electron

## Stereochemistry

Basic support for E/Z stereochemistry is included:

```rust
use amsr::bond::{Bond, BondStereo};

let e_bond = Bond::single().with_stereo(BondStereo::E);
let z_bond = Bond::single().with_stereo(BondStereo::Z);
```

## Error Handling

The library uses custom error types for robust error handling:

```rust
use amsr::errors::{AMSRResult, AMSRError};

fn process_molecule(smiles: &str) -> AMSRResult<()> {
    match smiles_to_amsr(smiles) {
        Ok(amsr) => {
            println!("Successfully converted to AMSR: {}", amsr);
            Ok(())
        },
        Err(AMSRError::InvalidSmiles(msg)) => {
            println!("Invalid SMILES: {}", msg);
            Ok(())
        },
        Err(e) => Err(e),
    }
}
```

## Examples

See the `examples/` directory for more detailed examples:

```bash
cargo run --example simple
```

## Limitations

This is a simplified implementation with the following limitations:

1. **Limited SMILES Support**: Basic SMILES parsing (no complex features like isotopes, stereochemistry, etc.)
2. **Simplified AMSR**: Not all AMSR features are implemented
3. **No 3D Conformers**: No support for 3D molecular conformations
4. **Basic Stereochemistry**: Limited stereochemistry support
5. **No Groups**: No support for molecular group abbreviations

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

This implementation is inspired by the original Python AMSR library and aims to provide similar functionality in Rust. 