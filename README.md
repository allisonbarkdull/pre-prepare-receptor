# pre-prepare-receptor

A Python tool to prepare receptor-ligand systems for docking.

## Features

- Step 0: Analyze ligand neighborhood in a receptor
  - Detect nearby waters, cofactors, histidines, polar rotors, amide flips, alternative locations, and other protonatable residues.
  - Generate PyMOL commands for visualization.

- Step 1: Prepare system for simulation
  - Fix missing atoms/residues using PDBFixer
  - Add hydrogens and handle alternative locations
  - Include selected waters and cofactors
  - Run a short equilibration with OpenMM, equilibrating hydrogen positions

---

## Installation

### Set up environment


Install the required packages using Micromamba:
```bash
micromamba install -c conda-forge openmm openmmtools openff-toolkit openmmforcefields espaloma pdbfixer parmed mdanalysis ambertools rdkit pandas deeptime pyemma
```
```bash
micromamba install -c mdtools cvpack
```
```bash
micromamba install -c conda-forge prody
```

# Clone and install dependencies

### Autopath
```bash
git clone git@github.com:forlilab/autopath.git
cd autopath
```
### Option 1: allison-change branch
```bash
git checkout allison-change
pip install -e .
```
### Option 2: lipids branch
```bash
git checkout lipids
pip install -e .
```
### Scrubber
```bash
git clone git@github.com:forlilab/scrubber.git
cd scrubber
pip install -e .
```
### Install pre-prepare-receptor

```bash
git clone git@github.com:allisonbarkdull/pre-prepare-receptor.git
cd pre-prepare-receptor
pip install -e .
```

---
---

## Usage

### Step 0: Analyze ligand neighborhood in a receptor
```bash
pre_prepare_receptor --mode step0 --input_pdb receptor.pdb --ligand_sdf ligand.sdf
```
### Step 1: Prepare system for simulation
```bash
pre_prepare_receptor --mode step1 --input_pdb receptor.pdb --ligand_sdf ligand.sdf
```

---

## Command-Line Interface (CLI)

### Required Arguments
--mode [step0|step1]  
  Choose the preparation step:  
    step0 → Analyze ligand neighborhood (waters, cofactors, histidines, alternative locations, flippable residues). Generates PyMOL script.  
    step1 → Prepare system: fix missing atoms/residues, add hydrogens, include waters/cofactors, set protonation, equilibrate with restraints.

--input_pdb <PDB_FILE>  
  Path to the receptor PDB file.

### Ligand Input Options
--ligand_sdf <SDF_FILE>  
  Path to the ligand SDF file.

OR

--ligand_pdb <SDF_FILE>  
  Path to the ligand PDB file. Hydrogens will be added with scrubber in step 1.

OR

--ligand_resname <RESNAME>  
  Ligand residue name in the PDB (e.g., ATP). Hydrogens will be added with scrubber in step 1.

--ligand_chain <CHAIN_ID>  
  Ligand chain ID in the PDB (e.g., A).

--ligand_resnum <RESNUM>  
  Ligand residue number (overrides resname selection).

--ligand_smiles <SMILES>  
  Optional SMILES string for the ligand. If SMILES is not passed will fetch from RCSB and hydrogens will be added with scrubber. If the SMILES has explicit hydrogens, this will be the protonaton state simluated

### Step 0 Neighborhood Options
--box_center <X Y Z>  
  Box center (x y z) in Å for Step 0 neighborhood analysis. Pass this instead of a ligand.

--box_lengths <X Y Z>  
  Box side lengths (x y z) in Å for Step 0 neighborhood analysis. Pass this instead of a ligand.

--cutoff <FLOAT>  
  Neighborhood cutoff distance in Å (default: 5.0).

### Waters and Cofactors
--water_residues <RESNUM_OR_CHAIN:RESNUM>  
  Residues of waters to keep (e.g., 101 A:105 Z:301).

--cofactor <NAME:VAL[:SMILES]>  
  Cofactor spec: NAME:VAL[:SMILES], VAL = sdf file or chain. Repeatable. 
                          "Examples: --cofactor NAD:cofactor_nad.sdf "
                            --cofactor EOH:B:'[O-]C([H])([H])C([H])([H])[H]' 
                            If the SMILES has explicit hydrogens, this will be the protonation state simulated.

--keep_all
keeps all waters and cofactors in the input pdb

### Protein Setup
--soluble_protein  
  Set up system for a soluble protein (no membrane). boo, tomato, tomato.

--add-missing-residues  
  Add missing residues using PDBFixer.

--keep_terminals  
  Build N- and C-terminal residues (may increase box size).

### Advanced Residue Handling
--altloc <CHAIN:RESNUM:ALT>  
  Specify alternative locations for residues. Repeatable.

--set_template <RESID:CHAIN:VARIANT>  
  Set residue variants for protonation states.
  Format: RESID:CHAIN:VARIANT (repeatable).
  RESID must be an integer, CHAIN a single letter, VARIANT the desired residue name.
  Example: --set_template 101:B:GLH 45:A HID
