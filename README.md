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

micromamba install -c conda-forge openmm openmmtools openff-toolkit openmmforcefields espaloma pdbfixer parmed mdanalysis ambertools rdkit pandas deeptime pyemma

micromamba install -c mdtools cvpack

micromamba install -c conda-forge prody

# Clone and install dependencies

### Autopath
git clone git@github.com:forlilab/autopath.git
cd autopath

# Option 1
git checkout allison-change
pip install -e .

# Option 2
git checkout lipids
pip install -e .

### Scrubber
git clone git@github.com:forlilab/scrubber.git
cd scrubber
pip install -e .

---

## Usage

### Step 0: Analyze ligand neighborhood in a receptor
python pre_prepare_receptor.py --mode step0 --input_pdb receptor.pdb --ligand_sdf ligand.sdf

### Step 1: Prepare system for simulation
python pre_prepare_receptor.py --mode step1 --input_pdb receptor.pdb --ligand_sdf ligand.sdf --output_dir prepared_system

Replace receptor.pdb and ligand.sdf with your files. The --output_dir specifies where to save prepared systems.

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

--ligand_resname <RESNAME>  
  Ligand residue name in the PDB (e.g., ATP).

--ligand_chain <CHAIN_ID>  
  Ligand chain ID in the PDB (e.g., A).

--ligand_resnum <RESNUM>  
  Ligand residue number (overrides resname selection).

--ligand_smiles <SMILES>  
  Optional SMILES string for ligand. If not provided, it will be fetched from RCSB and hydrogens added with scrubber. Explicit hydrogens define the protonation state.

### Step 0 Neighborhood Options
--box_center <X Y Z>  
  Box center for Step 0 neighborhood analysis (pass instead of ligand).

--box_lengths <X Y Z>  
  Box side lengths for Step 0 (pass instead of ligand).

--cutoff <FLOAT>  
  Neighborhood cutoff distance in Å (default: 5.0).

### Waters and Cofactors
--water_residues <RESNUM_OR_CHAIN:RESNUM>  
  Residues of waters to keep (e.g., 101 A:105 Z:301).

--cofactor <NAME:VAL[:SMILES]>  
  Cofactor specification. Repeatable.

### Protein Setup
--soluble_protein  
  Set up system for a soluble protein (no membrane).

--add-missing-residues  
  Add missing residues using PDBFixer.

--keep_terminals  
  Build N- and C-terminal residues (may increase box size).

### Advanced Residue Handling
--altloc <CHAIN:RESNUM:ALT>  
  Specify alternative locations for residues. Repeatable.

--set_template <RESID:CHAIN:VARIANT>  
  Set residue variants for protonation states. Repeatable.
