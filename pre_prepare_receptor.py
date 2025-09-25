#!/usr/bin/env python3
"""
to set up environment: 
micromamba install -c conda-forge openmm openmmtools openff-toolkit openmmforcefields espaloma pdbfixer parmed mdanalysis ambertools rdkit pandas deeptime pyemma
micromamba install -c mdtools cvpack
micromamba install prody
git clone git@github.com:forlilab/autopath.git
cd autopath
git checkout allison-change
pip install -e .
or
git checkout lipids
pip install -e .

git clone git@github.com:forlilab/scrubber.git
cd scrubber
pip install -e .
"""

from __future__ import annotations  # must be first

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from typing import List, Dict, Tuple, Optional
from pathlib import Path
import numpy as np
import requests

# RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    Chem = None
    AllChem = None

# OpenMM / PDBFixer / OpenFF / autopath
from openmm.app import PDBFile, Modeller, Topology
from openmm.unit import angstroms
from pdbfixer import PDBFixer
from openff.toolkit.topology import Molecule as OFFMolecule
from openff.toolkit import Topology as offTopology
from openff.units.openmm import to_openmm as offquantity_to_openmm
from autopath import SystemPreparation, Equilibration
from openmm import XmlSerializer

# ProDy
try:
    import prody as pr
    from prody import LOGGER
    LOGGER.verbosity = 'none'  # silence parser chatter
except ImportError:
    sys.stderr.write("[ERROR] prody is required. Install with: pip install prody\n")
    raise

# Scrubber (user-provided)
try:
    from scrubber import Scrub
except ImportError:
    Scrub = None

# Logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# Constants
STANDARD_AA = {
    'ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','PHE',
    'PRO','SER','THR','TRP','TYR','VAL','SEC','PYL'
}
HIS_VARIANTS = {'HIS','HID','HIE','HIP','HSD','HSE','HSP'}
WATER_NAMES = {'HOH','WAT','DOD','H2O'}
ResidueKey = Tuple[str, int, str]  # (chain, resi, resn)


# ---------------------------
# Utilities
# ---------------------------
def banner(title: str):
    print("\n" + "="*70)
    print(f"{title}")
    print("="*70 + "\n")



def run_reduce(input_pdb: str, output_pdb: str):
    output_pdb_path = Path(output_pdb)

    cmd = ["reduce", "-BUILD", "-FLIP", input_pdb]

    # Redirect stdout to the output file instead of the terminal
    with open(output_pdb_path, "w") as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)

    # Check if reduce exited cleanly
    # if result.returncode != 0:
    #     raise RuntimeError(f"Reduce failed:\n{result.stderr}")

    return str(output_pdb_path)


def get_ligand_smiles(ligand_id: str) -> Optional[str]:
    """Fetch canonical SMILES for a chemcomp from RCSB."""
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id.upper()}"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
    except Exception as e:
        logging.warning(f"Failed to fetch {ligand_id} from RCSB: {e}")
        return None

    data = r.json()
    # Try multiple possible keys
    for key in ("rcsb_chem_comp_descriptor", "purchasable_compound", "chem_comp"):
        smiles = data.get(key, {}).get("smiles") or data.get(key, {}).get("description")
        if smiles:
            return smiles
    logging.warning(f"SMILES not available for ligand {ligand_id}.")
    return None


def extract_residue_as_pdb(pdb_file: str, resname: str, chain: str, output_pdb: str) -> str:
    """
    Extracts residues matching (resname + chain) and writes to output_pdb (ProDy writePDB).
    If multiple matching residue copies present, all copies are written.
    """
    st = pr.parsePDB(pdb_file)
    if st is None:
        raise ValueError(f"Could not parse {pdb_file}")
    sel_str = f"chain {chain} and resname {resname}"
    sel = st.select(sel_str)
    if sel is None:
        raise ValueError(f"No residue {resname} on chain {chain} found in {pdb_file}")
    pr.writePDB(output_pdb, sel)
    logging.info(f"Extracted {resname} (chain {chain}) to {output_pdb}")
    return output_pdb


def create_rdkit_mol_from_smiles_and_pdb(smiles: str, pdb_file: str):
    """Create RDKit molecule from SMILES and PDB coordinates."""
    if Chem is None:
        raise RuntimeError("RDKit is required for this operation.")

    explicit_Hs = "[H]" in smiles
    template = Chem.MolFromSmiles(smiles)

    # first try: sanitize template
    try:
        Chem.SanitizeMol(template)
    except Exception:
        return _meeko_fallback(pdb_file)

    # load PDB molecule
    pdb_mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if pdb_mol is None:
        raise ValueError(f"Could not load PDB fragment: {pdb_file}")
    try:
        ligand = AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)
    except Exception:
        return _meeko_fallback(pdb_file)

    # skip scrubber if explicit Hs
    if explicit_Hs:
        return Chem.AddHs(ligand, addCoords=True)

    if Scrub is None:
        raise RuntimeError("Scrubber (Scrub) is required for protonation step.")

    scrub = Scrub(ph_low=7.4, ph_high=7.4, skip_gen3d=True,
                  skip_tautomers=True, skip_ringfix=True)
    scrubbed = next(scrub(ligand), None)
    if scrubbed is None:
        raise ValueError("Scrubber failed to produce a molecule.")

    mol_h = Chem.AddHs(scrubbed, addCoords=True)
    return mol_h


def _meeko_fallback(pdb_file):
    """Fallback using Meeko for tricky ligands like heme."""
    print(
        "Using Meeko fallback. Double check protonation state!"
    )
    from meeko import Polymer, ResidueChemTemplates, MoleculePreparation
    with open(pdb_file) as f:
        pdb_str = f.read()
    templates = ResidueChemTemplates.create_from_defaults()
    mk_prep = MoleculePreparation()
    polymer = Polymer.from_pdb_string(pdb_str, templates, mk_prep)

    chain, resid = next(iter(polymer.monomers.keys())).split(":")
    monomer = polymer.monomers[f"{chain}:{resid}"]
    return monomer.rdkit_mol

def make_sdf_from_residue(
    pdb_file: str, resname: str, chain: str, output_sdf: str,
    smiles: Optional[str] = None, ligand_pdb: Optional[str] = None
) -> str:
    """Extract a residue and save as SDF. Can accept a ligand PDB directly."""
    tmp_dir = tempfile.mkdtemp(prefix="res_extract_")
    if ligand_pdb:
        if smiles:
            logging.info(f"Using provided SMILES for {resname}:{chain}")
        else:
            st = pr.parsePDB(ligand_pdb)
            resname = list({res.getResname() for res in st.iterResidues()})[0]
            smiles = get_ligand_smiles(resname)
    try:
        if not ligand_pdb:
            tmp_pdb = os.path.join(tmp_dir, f"{resname}_{chain}.pdb")
            extract_residue_as_pdb(pdb_file, resname, chain, tmp_pdb)

            # Use passed SMILES if provided; otherwise fetch from RCSB
            if smiles:
                logging.info(f"Using provided SMILES for {resname}:{chain}")
            else:
                smiles = get_ligand_smiles(resname)
                if not smiles:
                    raise ValueError(f"Cannot fetch SMILES for {resname}; aborting SDF creation.")
        print(f' smiles: {smiles}')
        # Create RDKit molecule
        if ligand_pdb:
            rdkit_mol = create_rdkit_mol_from_smiles_and_pdb(smiles, ligand_pdb)
            if rdkit_mol is None:
                raise ValueError(f"Could not load ligand PDB: {ligand_pdb}")
        else:
            rdkit_mol = create_rdkit_mol_from_smiles_and_pdb(smiles, tmp_pdb)

        writer = Chem.SDWriter(output_sdf)
        writer.write(rdkit_mol)
        writer.close()
        logging.info(f"Wrote SDF for {resname}:{chain} -> {output_sdf}")
    finally:
        if not ligand_pdb:
            shutil.rmtree(tmp_dir)
    return output_sdf


# ---------------------------
# STEP 0 FUNCTIONS (neighborhood)
# ---------------------------
def ligand_coords_from_file(ligand_path: str) -> np.ndarray:
    """Read ligand coordinates from PDB or SDF and return as Nx3 NumPy array."""
    ext = os.path.splitext(ligand_path)[1].lower()

    if ext == ".sdf":
        mol = Chem.SDMolSupplier(ligand_path, removeHs=False)[0]
        if mol is None:
            raise ValueError(f"Could not read ligand from {ligand_path}")
        conf = mol.GetConformer()
        coords = np.array(
            [[conf.GetAtomPosition(i).x,
              conf.GetAtomPosition(i).y,
              conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())],
            dtype=float
        )
        return select

    elif ext == ".pdb":
        st = pr.parsePDB(ligand_path)
        if st is None:
            raise ValueError(f"Cannot parse ligand PDB: {ligand_path}")
        return np.asarray(st.getCoords(), dtype=float)

    else:
        raise ValueError(f"Unsupported ligand file type: {ligand_path}")


def get_residue_key(atom) -> ResidueKey:
    resn = atom.getResname().strip()
    try:
        resi = int(atom.getResnum())
    except Exception:
        resi = 0
    chain = atom.getChid() or ''
    return (chain, resi, resn)


def unique_residues_from_atoms(atoms: pr.AtomGroup) -> List[ResidueKey]:
    """Return unique residues from ProDy AtomGroup."""
    if atoms is None:
        return []
    seen = set()
    out = []
    for a in atoms:
        key = (a.getChid() or '', int(a.getResnum()), a.getResname().strip())
        if key[:2] in seen:
            continue
        seen.add(key[:2])
        out.append(key)
    return sorted(out, key=lambda x: (x[0], x[1]))


def atoms_within_cutoff(receptor: pr.AtomGroup, ligand_coords: np.ndarray, cutoff: float) -> pr.AtomGroup:
    rec_coords = receptor.getCoords()
    if rec_coords is None:
        raise ValueError("Receptor has no coordinates.")
    idxs = []
    chunk = 50000
    for start in range(0, rec_coords.shape[0], chunk):
        rc = rec_coords[start:start+chunk]
        d2 = np.min(np.sum((rc[:, None, :] - ligand_coords[None, :, :])**2, axis=2), axis=1)
        sel = np.where(d2 <= cutoff**2)[0]
        idxs.extend((start + sel).tolist())
    if not idxs:
        return pr.AtomGroup('empty')
    return receptor[idxs]


def select_neighborhood(receptor: pr.AtomGroup, args) -> Tuple[Dict[str,List[ResidueKey]], str]:
    cutoff = args.cutoff if args.mode == 'step0' else None
    ligand_selector = 'ligand_site'

    if args.box_center and args.box_lengths:
        cx, cy, cz = args.box_center
        lx, ly, lz = args.box_lengths
        mins = np.array([cx - lx/2, cy - ly/2, cz - lz/2])
        maxs = np.array([cx + lx/2, cy + ly/2, cz + lz/2])

        rec_coords = receptor.getCoords()
        in_box = np.all((rec_coords >= mins) & (rec_coords <= maxs), axis=1)
        idxs = np.where(in_box)[0]
        neigh_atoms = receptor[idxs] if idxs.size > 0 else pr.AtomGroup('empty')
        ligand_selector = f"(x between {mins[0]:.1f} and {maxs[0]:.1f} and " \
                          f"y between {mins[1]:.1f} and {maxs[1]:.1f} and " \
                          f"z between {mins[2]:.1f} and {maxs[2]:.1f})"
    elif args.ligand_sdf or args.ligand_pdb:
        lig_path = args.ligand_sdf if args.ligand_sdf else args.ligand_pdb
        lig_coords = ligand_coords_from_file(lig_path)
        neigh_atoms = atoms_within_cutoff(receptor, lig_coords, cutoff if cutoff else np.inf)
        ligand_selector = f"all within {args.cutoff:.1f} of (ligand)"
    else:
        # must have ligand_resname and ligand_chain (and optionally ligand_resnum)
        if not (args.ligand_resname and args.ligand_chain):
            raise RuntimeError("Step 0 requires either --ligand_sdf or (--ligand_resname + --ligand_chain).")
        if args.ligand_resnum:
            sel_str = f"chain {args.ligand_chain} and resnum {args.ligand_resnum}"
        else:
            sel_str = f"chain {args.ligand_chain} and resname {args.ligand_resname}"
        lig_atoms = receptor.select(sel_str)
        if lig_atoms is None or lig_atoms.numAtoms() == 0:
            raise ValueError(f"Could not find ligand with selection: {sel_str}")
        ligand_selector = sel_str
        if cutoff:
            neigh_atoms = receptor.select(f"within {cutoff} of ({sel_str})")
        else:
            neigh_atoms = receptor

    if neigh_atoms is None:
        neigh_atoms = pr.AtomGroup('empty')

    cats = {}
    water_atoms = (neigh_atoms if cutoff else receptor).select('water')
    cats['waters'] = unique_residues_from_atoms(water_atoms) if water_atoms is not None else []
    het_atoms = (neigh_atoms if cutoff else receptor).select('hetero and not water')
    cof_keys = []
    if het_atoms is not None:
        for rk in unique_residues_from_atoms(het_atoms):
            if rk[2] not in STANDARD_AA and rk[2] not in WATER_NAMES and  rk[2] != "DUM" :
                cof_keys.append(rk)
    cats['cofactors'] = cof_keys
    his_atoms = (neigh_atoms if cutoff else receptor).select('resname ' + ' '.join(HIS_VARIANTS))
    cats['histidines'] = unique_residues_from_atoms(his_atoms) if his_atoms is not None else []
    ysth_atoms = (neigh_atoms if cutoff else receptor).select('resname TYR SER THR ' + ' '.join(HIS_VARIANTS))
    cats['polar_rotors'] = unique_residues_from_atoms(ysth_atoms) if ysth_atoms is not None else []
    nq_atoms = (neigh_atoms if cutoff else receptor).select('resname ASN GLN')
    cats['amide_flips'] = unique_residues_from_atoms(nq_atoms) if nq_atoms is not None else []
    altloc_atoms = [a for a in neigh_atoms if getattr(a, 'getAltLoc', lambda: '')() not in ('', None)]
    cats['altloc_residues'] = unique_residues_from_atoms(pr.AtomGroup(altloc_atoms)) if altloc_atoms else []
    other_protonatable_resnames = ['ASP', 'ASH', 'GLU','GLH' 'CYS', 'CYX', 'LYN', 'LYS']  # plus neutral variants
    prot_atoms = (neigh_atoms if cutoff else receptor).select('resname ' + ' '.join(other_protonatable_resnames))
    cats['other_protonatable'] = unique_residues_from_atoms(prot_atoms) if prot_atoms is not None else []

    return cats, ligand_selector


def format_res_list(reslist: List[ResidueKey], category: str = "") -> str:
    if not reslist:
        return "  (none)"
    lines = []
    for chain, resi, resn in reslist:
        chain_disp = chain if chain else "_"
        if category.lower() == "waters":
            lines.append(f"  Water {resn} {chain_disp}:{resi}")
        else:
            lines.append(f"  {resn} {chain_disp}:{resi}")
    return "\n".join(lines)


def pymol_commands(
    cats: Dict[str, List[ResidueKey]], 
    ligand_selector: str, 
    cutoff: float = 5.0
) -> str:
    """
    Generate PyMOL commands for defining a ligand site and nearby categories.
    
    Parameters
    ----------
    cats : dict
        Dictionary of residue categories.
    ligand_selector : str
        PyMOL selection string for the ligand (e.g., "resn LIG").
    cutoff : float, optional
        Distance cutoff in Å for defining the ligand site.
    """
    lines = []
    suffix = "_near"

    # Define ligand site as everything within cutoff Å of the ligand
    lines.append(f"# Define ligand site")
    lines.append(f"select ligand_site, (byres ({ligand_selector}) around {cutoff})")
    def res_expr(reslist):
        parts = []
        for chain, resi, resn in reslist:
            if chain:
                parts.append(f"(chain {chain} and resi {resi} and resn {resn})")
            else:
                parts.append(f"(resi {resi} and resn {resn})")
        return " or ".join(parts) if parts else "none"
    lines.append(f"create waters{suffix}, {res_expr(cats['waters'])}")
    lines.append(f"create cofactors{suffix}, {res_expr(cats['cofactors'])}")
    lines.append(f"create histidines{suffix}, {res_expr(cats['histidines'])}")
    lines.append(f"create polar_rotors{suffix}, {res_expr(cats['polar_rotors'])}")
    lines.append(f"create amide_flips{suffix}, {res_expr(cats['amide_flips'])}")
    lines.append(f"create altloc_residues{suffix}, {res_expr(cats['altloc_residues'])}")
    lines.append(f"create other_protonatable{suffix}, {res_expr(cats['other_protonatable'])}")
    return "\n".join(lines)


# ---------------------------
# STEP 1 functions (PDB fix, waters/cofactors insertion)
# ---------------------------

def filter_altlocs(input_pdb: str, altloc_dict: Dict[Tuple[str,int], str], output_pdb: str):
    """Keep only specified altLocs; remove others."""
    with open(input_pdb) as f:
        lines = f.readlines()

    new_lines = []
    for l in lines:
        if l.startswith("ATOM") or l.startswith("HETATM"):
            chain = l[21].strip()
            resnum = int(l[22:26].strip())
            altloc = l[16].strip()  # column 17
            key = (chain, resnum)
            if key in altloc_dict:
                if altloc == altloc_dict[key] or altloc == " ":
                    new_lines.append(l)
            else:
                new_lines.append(l)
        else:
            new_lines.append(l)

    with open(output_pdb, "w") as f:
        f.writelines(new_lines)

def fix_pdb(
    input_pdb: str,
    output_pdb: str,
    add_missing_residues: bool = False,
    skip_terminal_residues: bool = True
):
    """
    Fixes a PDB file using PDBFixer.

    Args:
        input_pdb (str): Input PDB filename
        output_pdb (str): Output PDB filename
        add_missing_residues (bool): If True, attempt to add missing residues
        skip_terminal_residues (bool): If True, don't add missing residues
                                       at the N- or C-terminus
    """
    logging.info(f"Fixing PDB: {input_pdb}")
    fixer = PDBFixer(filename=input_pdb)

    # --- handle missing residues ---
    fixer.findMissingResidues()  # always run, ensures attribute exists
    if not add_missing_residues:
        fixer.missingResidues = {}  # clear everything
    elif skip_terminal_residues:
        to_remove = []
        for key, residues in fixer.missingResidues.items():
            chain_index, insert_index = key
            chain = list(fixer.topology.chains())[chain_index]
            if insert_index == 0 or insert_index == len(list(chain.residues())):
                to_remove.append(key)
        for key in to_remove:
            del fixer.missingResidues[key]

    # --- rest of standard PDBFixer pipeline ---
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    with open(output_pdb, "w") as out_file:
        PDBFile.writeFile(fixer.topology, fixer.positions, out_file, keepIds=True)

    logging.info(f"Fixed PDB written to: {output_pdb}")


WATER_NAMES = {"HOH", "WAT", "H2O"}  # extend if needed

def add_selected_waters_and_varients(
    original_pdb_path: str, 
    fixed_pdb_path: str, 
    output_path: str, 
    water_residues: List[Tuple[Optional[str], int]],
    variants: Optional[Dict[Tuple[str,int], str]] = None
):
    """
    Add selected waters from the original PDB to the fixed PDB,
    preserving original chains.

    water_residues: list of (chain, resnum)
        chain=None means keep all chains with that resnum.
    """
    logging.info(f"Adding selected waters (residues): {water_residues}")

    original_pdb = PDBFile(original_pdb_path)
    fixed_pdb = PDBFile(fixed_pdb_path)
    modeller = Modeller(fixed_pdb.topology, fixed_pdb.positions)

    selected_atoms = []
    selected_positions = []
    for chain in original_pdb.topology.chains():
        for res in chain.residues():
            try:
                resid_int = int(res.id)
            except Exception:
                continue
            if res.name.upper() not in WATER_NAMES:
                continue
            # Check if this water is in the user list
            keep = False
            for w_chain, w_resnum in water_residues:
                if w_chain is None:
                    if resid_int == w_resnum:
                        keep = True
                        break
                else:
                    if resid_int == w_resnum and chain.id.upper() == w_chain.upper():
                        keep = True
                        break
            if not keep:
                continue
            # Add atoms of this water
            for atom in res.atoms():
                selected_atoms.append(atom)
                selected_positions.append(original_pdb.positions[atom.index])

    if selected_atoms:
        water_top = Topology()
        chain_map = {}  # original chain -> new chain object
        for atom in selected_atoms:
            chain_id = atom.residue.chain.id
            if chain_id not in chain_map:
                chain_map[chain_id] = water_top.addChain(id=chain_id)
            current_chain = chain_map[chain_id]
            # Add residue if not already added
            try:
                current_res = next(r for r in current_chain.residues() if r.id == atom.residue.id)
            except StopIteration:
                current_res = water_top.addResidue(atom.residue.name, current_chain, id=atom.residue.id)
            water_top.addAtom(atom.name, atom.element, current_res)

        modeller.add(water_top, selected_positions)

    # Pass the variants dict to addHydrogens
    modeller.addHydrogens(variants=[
        variants.get((res.chain.id, int(res.id))) if variants else None
        for res in modeller.topology.residues()
    ])

    with open(output_path, "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)

    logging.info(f"Finished adding waters. Output written to: {output_path}")

def write_equilibration_json(water_residues_to_keep=None, cofactor_resnames=None, output_file="equil.json"):
    """
    Creates or updates equilibration JSON file with protein, optional waters, and cofactors.

    water_residues_to_keep: list of (chain, resid) tuples, where chain can be None for any chain
    cofactor_resnames: list of cofactor resnames
    """
    # Base selection: protein without hydrogens
    selection = "protein and not name H*"

    # Add water residues if provided
    if water_residues_to_keep:
        water_terms = []
        for chain, resid in water_residues_to_keep:
            water_terms.append(f"(resname HOH and resid {resid} and not name H*)")
        selection += " or (" + " or ".join(water_terms) + ")"

    # Add cofactors if provided
    if cofactor_resnames:
        cofactor_terms = [f"(resname {n} and not name H*)" for n in cofactor_resnames]
        selection += " or (" + " or ".join(cofactor_terms) + ")"

    # Build JSON structure
    equil_data = {
        "components_lookup": {
            "everything": selection
        },
        "minimization": [
            {"name": "Minimization Stage1", "forces": [100.0]},
            {"name": "Minimization Stage2", "forces": [500.0]},
            {"name": "Minimization Stage2", "forces": [500.0]},
            {"name": "Minimization Stage2", "forces": [500.0]},
            {"name": "Minimization Stage2", "forces": [500.0]},
        ],
        "warmup": {
            "T_initial": 100.0,
            "T_final": 310.0,
            "T_step": 5.0,
            "npt_flag": False,
            "nsteps": 100_000,
            "stepsize": 0.002
        },
        "equilibration": [
            {"name": "Equilibration Stage1", "forces": [500.0], "npt_flag": False, "nsteps": 50_000, "stepsize": 0.002},
            {"name": "Equilibration Stage2", "forces": [500.0], "npt_flag": True, "nsteps": 50_000, "stepsize": 0.002}
        ]
    }

    # Write JSON to file
    with open(output_file, "w") as f:
        json.dump(equil_data, f, indent=2)

    logging.info(f"Equilibration JSON written to {output_file} with selection: {selection}")

def combine_protein_and_ligands(input_pdb, ligands_to_parametrize, output_pdb, allow_undefined_stereo=True):
    """
    Combine a fixed protein PDB with ligands/cofactors from SDFs into a single PDB.
    
    Args:
        input_pdb (str): path to fixed protein PDB
        ligands_to_parametrize (list of tuples): [(resname, sdf_path), ...]
        output_pdb (str): path to write combined PDB
        allow_undefined_stereo (bool): OpenFF Molecule option
    """
    # Load protein
    pdb = PDBFile(input_pdb)
    modeller = Modeller(pdb.topology, pdb.positions)

    # Load each ligand SDF, convert to OpenFF molecule, then OpenMM topology
    for resname, sdf_path in ligands_to_parametrize:
        rdkit_mol = Chem.SDMolSupplier(sdf_path, removeHs=False)[0]
        ligand = OFFMolecule.from_rdkit(rdkit_mol, allow_undefined_stereo=allow_undefined_stereo)
        ligand_off_topology = offTopology.from_molecules(molecules=[ligand])
        ligand_omm_topology = ligand_off_topology.to_openmm()
        for res in ligand_omm_topology.residues():
            res.name = resname
        ligand_positions = offquantity_to_openmm(ligand.conformers[0])

        # Add ligand to modeller
        modeller.add(ligand_omm_topology, ligand_positions)

    # Write combined PDB
    with open(output_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)
    print(f"→ Combined protein + ligands written to: {output_pdb}")


def clean_pdb(input_pdb: str, output_pdb: str, water_residues: List[Tuple[Optional[str], int]] = None):
    """
    Clean a PDB file by removing unwanted residues while optionally keeping specific waters.

    Parameters
    ----------
    input_pdb : str
        Path to the input PDB file.
    output_pdb : str
        Path to the output cleaned PDB file.
    water_residues : List[Tuple[Optional[str], int]], optional
        List of water residues to keep. Each entry is a tuple of (chain, resid).
        - If `chain` is None, the residue number will be matched across all chains.

    Behavior
    --------
    - Removes all waters (HOH), except those explicitly requested in `water_residues`.
    - Removes POP, Na, and Cl residues.
    - Prints which waters are kept and how many atoms are removed.
    - Writes the cleaned structure to `output_pdb`.

    Raises
    ------
    ValueError
        If the input PDB cannot be parsed.
    """
    st = pr.parsePDB(input_pdb)
    if st is None:
        raise ValueError(f"Cannot parse {input_pdb}")

    water_residues = water_residues or []
    to_remove_sel = "(resname HOH) or resname POP or resname NA or resname CL" # Base selection of unwanted atoms
    if water_residues:
        keep_terms = []
        for chain, resid in water_residues:
            if chain:
                keep_terms.append(f"(resname HOH and chain {chain} and resnum {resid})")
            else:
                keep_terms.append(f"(resname HOH and resnum {resid})")
        keep_clause = " or ".join(keep_terms)
        to_remove_sel = f"(({to_remove_sel}) and not ({keep_clause}))"

        print("[INFO] Keeping waters:", ", ".join(
            [f"{chain or '_'}:{resid}" for chain, resid in water_residues]
        ))

    remove_atoms = st.select(to_remove_sel)
    if remove_atoms is not None and remove_atoms.numAtoms() > 0:
        print(f"[WARNING] Removing {remove_atoms.numAtoms()} atoms: waters not selected, POP, Na, Cl.")
        keep_sel = f"not ({to_remove_sel})"
        st = st.select(keep_sel)

    pr.writePDB(output_pdb, st)
    print(f"→ Cleaned PDB written to: {output_pdb}")


# ---------------------------
# MAIN
# ---------------------------
def main():

    base_out_dir = "pre_prepare_receptor_outputs"
    md_out_dir = os.path.join(base_out_dir, "md")
    os.makedirs(base_out_dir, exist_ok=True)
    parser = argparse.ArgumentParser(
        description=(
            "pre_prepare_receptor.py: A tool to prepare receptor-ligand systems for docking.\n\n"
            "This script has two main modes:\n"
            "  step0 → Analyze the neighborhood of a ligand in the receptor.\n"
            "           Detect nearby waters, cofactors, histidines, protonatable residues, alternative locations, and flippable residues.\n"
            "           Generates a PyMOL selection script for visualization.\n"
            "  step1 → Prepare the system and minimize hydrogen positions.\n"
            "           Fix missing atoms/residues, add hydrogens, include selected waters/cofactors, apply residue protonation states, handle alternative locations.\n"
            "           Build a solvated/membrane system and equilibrate with restraints.\n\n"
            "Usage examples are included in the script header."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    # Required inputs
    parser.add_argument(
        "--mode", 
        choices=["step0", "step1"], 
        required=True,
        help=(
            "Choose the preparation step to run:\n"
            "  step0 → Analyze the neighborhood of a ligand in the receptor.\n"
            "           Detect nearby waters, cofactors, histidines, protonatable residues, alternative locations, flippable residues.\n"
            "           Generates a PyMOL selection script for visualization.\n"
            "  step1 → Prepare the system and minimize hydrogen positions with OPENMM.\n"
            "           Fix missing atoms/residues, add hydrogens, include selected waters/cofactors, residue protonation states, alternative locations.\n"
            "           Build solvated/membrane system and equilibrate with restraints."
        )
    )
    parser.add_argument("--input_pdb", type=str, required=True, help='PDB of the protein to prepare')


    # --- Ligand input options ---
    parser.add_argument("--ligand_sdf", type=str, default=None,
                        help="Ligand SDF path. Input a ligand sdf or a ligand pdb or ligand_resname + ligand_chain")
    parser.add_argument("--ligand_pdb", type=str, default=None,
                    help="Ligand PDB path. Input a ligand sdf or a ligand pdb or ligand_resname + ligand_chain")                    
    parser.add_argument("--ligand_resname", type=str, default=None,
                        help="Ligand residue name in PDB (e.g., ATP)")
    parser.add_argument("--ligand_chain", type=str, default=None, help="Ligand chain ID (e.g., A)")
    parser.add_argument("--ligand_resnum", type=int, default=None,
                        help="Ligand residue number (optional, overrides resname selection)")
    parser.add_argument("--ligand_smiles", type=str, default=None,
                        help="Optional SMILES string for the ligand. If SMILES is not passed will fetch from RCSB and add hydrogens with scrubber. If the SMILES has explicit hydrogens, this will be the protonaton state simluated")

    # --- Step 0 neighborhood options ---
    parser.add_argument("--box_center", type=float, nargs=3, default=None,
                        help="Box center (x y z) in Å for Step 0 neighborhood analysis. Pass this instead of a ligand.")
    parser.add_argument("--box_lengths", type=float, nargs=3, default=None,
                        help="Box side lengths (x y z) in Å for Step 0 neighborhood analysis. Pass this instead of a ligand.")
    parser.add_argument("--cutoff", type=float, default=5.0, help="Neighborhood cutoff (Å) for Step 0")

    # --- Waters and cofactors ---
    parser.add_argument("--water_residues", type=str, nargs='*', default=[],
                    help="Residue numbers or CHAIN:RESNUM of waters to keep (e.g., --water_residues 101 A:105 Z:301)")
    parser.add_argument("--cofactor", action='append', default=[],
                        help=(
                            "Cofactor spec: NAME:VAL[:SMILES], VAL = sdf file or chain. Repeatable. "
                            "Examples: --cofactor NAD:cofactor_nad.sdf "
                            "--cofactor EOH:B:'[O-]C([H])([H])C([H])([H])[H]' "
                            "If the SMILES has explicit hydrogens, this will be the protonaton state simluated."))
    parser.add_argument("--keep_all", action="store_true", help="If set in step1, keep ALL waters and ALL HETATM residues (except the ligand) as cofactors.")

    # --- Protein setup ---
    parser.add_argument("--soluble_protein", action='store_true',
                        help="Set up the system for a soluble protein; no membrane will be added. boo.")
    parser.add_argument("--add-missing-residues", action="store_true",
                        help="Add missing residues using PDBFixer. By default, missing residues are not added.") #to do, I think this is a lie :(. I think I am not adding missing LOOPS but still adding missing structured residues
    parser.add_argument("--keep_terminals", action='store_true',
                        help="Use pdbfixer to build N and C terminals of your protein in addtion to missing loops. This can cause your box size to be unecessarilty large")

    # --- Advanced residue handling ---
    parser.add_argument("--altloc", action='append', default=[],
                        help="Specify alternative locations for residues. Format: CHAIN:RESNUM:ALT (repeatable). "
                             "Example: --altloc A:101:A --altloc B:205:B")
    parser.add_argument("--set_template", type=str, nargs='*', default=[],
                        help=(
                            "Set residue variants for protonation states. "
                            "Format: RESID:CHAIN:VARIANT (repeatable). "
                            "RESID must be an integer, CHAIN a single letter, VARIANT the desired residue name. "
                            "Example: --set_template 101:B:GLH 45:A HID"
                        ))

    args = parser.parse_args()

    # Prefer reduced receptor_reduced.pdb if present (created by Step0)
    reduced_pdb_path = f"{base_out_dir}/0_1_receptor_reduced.pdb"

    if args.mode == "step0":
        banner("STEP 0: Neighborhood analysis and initial hydrogenation")
        print("Step 0 will perform the following:")
        print("  1) Run Reduce to add hydrogens and optimize orientations.")
        print("     → Output: receptor_reduced.pdb  (hydrogenated receptor for Step 1)")
        print("  2) Analyze neighborhood of ligand or binding site:")
        print("     → Identify nearby waters, cofactors, histidines, polar rotors, amide flips, altlocs, other protonatable residues.")
        print("  3) Generate PyMOL commands for visualization:")
        print("     → Output: 0_pymol_cmds.pml  (contains selections for waters, cofactors, histidines, polar residues, etc.)\n")

        # validate ligand options
        if args.ligand_sdf:
            print(f"[INFO] Using ligand from SDF: {args.ligand_sdf}")
        elif args.ligand_pdb:
            print(f"[INFO] Using ligand from SDF: {args.ligand_pdb}")
        elif args.ligand_resname and args.ligand_chain:
            print(f"[INFO] Using ligand by residue selection: {args.ligand_resname} chain {args.ligand_chain} "
                f"(resnum override: {args.ligand_resnum})")
        elif args.box_center and args.box_lengths:
            print(f"[INFO] Using box-based binding site: center={args.box_center}, lengths={args.box_lengths}")
        else:
            print("[ERROR] Step 0 requires one of: --ligand_sdf, ligand_pdb, (--ligand_resname + --ligand_chain), or "
                "(--box_center + --box_lengths).")
            sys.exit(1)

        # run Reduce to add hydrogens and optimize sidechain orientations (write to receptor_reduced.pdb)
        print("→ Running Reduce to add hydrogens and optimize orientations (this requires 'reduce' in PATH)...")
        run_reduce(args.input_pdb, reduced_pdb_path)
        print(f"→ Hydrogenated receptor saved as: {reduced_pdb_path}")
        print("Use this reduced file as input to Step 1 (the script will auto-detect it).")

        # Use the reduced structure for neighborhood analysis
        receptor = pr.parsePDB(reduced_pdb_path)
        cats, ligand_sel = select_neighborhood(receptor, args)


        ##to do add if so we only print if that thing is found? mayb
        print(f"\nNeighborhood cutoff: {args.cutoff:.1f} Å around the ligand (if in docking mode).\n")
        print("Found waters within cutoff — decide whether to KEEP or REMOVE depending on structural role:")
        print(format_res_list(cats['waters'], category="waters"))
        print("You can keep waters in step 1 with --water_residues")
        print("\nFound other cofactors within cutoff — decide case-by-case (metals, lipids, ions). You can keep cofactors in step 1 with --cofactors. Also think about the protonation state - step1 will protonate with scrubber but this could be incorrect:")
        print(format_res_list(cats['cofactors']))
        print("\nHistidines within cutoff — choose proper protonation state (HID/HIE/HIP).")
        print("Reduce/PDBFixer will attempt canonical choices in Step 1, but PLEASE double-check. You can choose protoantion state in step 1 with --set_template")
        print(format_res_list(cats['histidines']))
        print("\nTyr/Ser/Thr/His within cutoff — orient phenolic/serine/threonine/imidazole hydrogens for plausible H-bonds. Step 1 will minimize these positions with MD:")
        print(format_res_list(cats['polar_rotors']))
        print("\nAsn/Gln within cutoff — consider flipping amide to improve H-bond geometry:  Step 1 will take the reduced structure as input if it exsits by default, so update this if you want.")
        print(format_res_list(cats['amide_flips']))
        if cats['altloc_residues']:
            print("\n[WARNING] Residues with alternative locations within cutoff. You an choose which location you want in step 1 with --altloc:")
            print(format_res_list(cats['altloc_residues']))

        if cats['other_protonatable']:
            print("\n[NOTE] Other protonatable residues (Asp, Glu, Cys, Lys) within cutoff. You can choose protoantion state in step 1 with --set_template:")
            print(format_res_list(cats['other_protonatable']))

        cmds = pymol_commands(cats, ligand_sel, args.cutoff)
        with open(f"{base_out_dir}/0_pymol_cmds.pml", "w") as f:
            f.write(cmds)
        print("\nSaved PyMOL commands to 0_pymol_cmds.pml for visualization.")
        print("Step 0 complete.\n")

    elif args.mode == "step1":
        os.makedirs(md_out_dir, exist_ok=True)

        banner("STEP 1: System preparation & equilibration")
        print("This step will prepare the protein + ligand + cofactor system and run a short MD simulation with all heavy atoms restrained with the goal of optimizing hydrogen positons for docking.")
        print("\nFiles that may be generated in Step 1:")
        filtered_altloc_pdb = os.path.join(base_out_dir, "1_0_altloc_filtered.pdb")
        print(f"  1. AltLoc-filtered PDB: {filtered_altloc_pdb}")
        print("       - Only generated if --altloc is specified")
        if hasattr(args, "altloc") and args.altloc:
            print(f"       - AltLocs included: {args.altloc}")
        fixed_pdb = os.path.join(base_out_dir, "1_1_input_pdb_fixer.pdb")
        print(f"  2. Protein PDB after PDBFixer: {fixed_pdb}")
        print("       - This will have all waters and cofactors removed. All missing protein atoms added, hydrogens added at pH 7.4 (user-specified protonation states not yet applied)")
        if args.add_missing_residues or args.keep_terminals:
            print("       - Missing loops and optionally terminal residues added")
        with_water_pdb = os.path.join(base_out_dir, "1_2_water_and_variants.pdb")
        print(f"  3. PDB with waters and side-chain variants: {with_water_pdb}")
        print("       - Only generated if -water_residues and/or --set_template is specified")
        if args.water_residues:
            print(f"       - Includes selected waters: {args.water_residues}")
        if args.set_template:
            print(f"       - Includes selected sidechain templates: {args.set_template}")

        combined_pdb = os.path.join(base_out_dir, "1_3_md_input_no_solvent.pdb")
        print(f"  4. PDB with altlocs, user-specifed protonaton, ligand, cofactor, waters: {combined_pdb}")
        print("       - Generated if any ligands or cofactors are provided")
        if args.cofactor:
            print(f"       - Includes specified cofactors: {args.cofactor}")
        if args.ligand_sdf:
            print(f"       - Includes ligand: {args.ligand_sdf}")
        if args.ligand_pdb:
            print(f"       - Includes ligand: {args.ligand_pdb}")            
        if args.ligand_resname:
            print(f"       - Includes ligand resname: {args.ligand_resname}")

        minimized_pdb = os.path.join(base_out_dir, "1_4_minimized.pdb")
        print(f"  5. Minimized system PDB: {minimized_pdb}")
        print("       - Output of MD minimization with all MD solvents (lipid, waters you didn't specify) removed")

        # 6. Cleaned equilibrated PDB
        final_pdb = os.path.join(base_out_dir, "1_5_final.pdb")
        print(f"  6. Equilibrated system PDB: {final_pdb}")
        print("       - Output of MD equilibration with all MD solvents (lipid, waters you didn't specify) removed. You should check this, but hopefully this is an approprite input for docking!! <3")

        print(f"Also, the md files (input .pdb with solvent/lipid, .xml, .prmtop, .dcd, and output .pdbs are saved in {md_out_dir})")

        print("------------------------------------------------------------\n")
        # Choose input PDB for fixing — prefer reduced if present
        input_pdb_for_fix = reduced_pdb_path if os.path.exists(reduced_pdb_path) else args.input_pdb
        if os.path.exists(reduced_pdb_path):
            print(f"→ Found reduced receptor from Step 0: {reduced_pdb_path} (using this as input).")
        else:
            print(f"→ No reduced receptor found; using provided PDB: {args.input_pdb}.")

        if args.altloc:
            filtered_pdb = os.path.join(base_out_dir, "1_0_altloc_filtered.pdb")
            filter_altlocs(input_pdb_for_fix, altloc_dict, filtered_pdb)
            input_pdb_for_fix = filtered_pdb


        print("→ Running PDBFixer to repair coordinates, add missing atoms/hydrogens...")
        fix_pdb(
        input_pdb=input_pdb_for_fix,
        output_pdb=fixed_pdb,
        add_missing_residues=args.add_missing_residues,
        skip_terminal_residues=not args.keep_terminals
    )

        template_variants: Dict[Tuple[str,int], str] = {}  # key=(chain,resid), value=variant

        for entry in args.set_template:
            try:
                resid_str, chain, variant = entry.split(":")
                resid = int(resid_str)
                template_variants[(chain, resid)] = variant
            except Exception as e:
                logging.warning(f"Skipping invalid --set_template entry '{entry}': expected format RESID:CHAIN:VARIANT ({e})")
        
        water_residues = []  
        for item in args.water_residues:
            if ':' in item:
                parts = item.split(':')
                if len(parts) != 2:
                    raise ValueError(f"Invalid water residue spec '{item}'; use RESNUM or CHAIN:RESNUM")
                chain, resnum = parts
                water_residues.append((chain.upper(), int(resnum)))
            else:
                water_residues.append((None, int(item)))

        if args.keep_all:
            cofactor = []
            logging.info("`--keep_all` enabled: keeping ALL waters and ALL non-water HETATM residues as cofactors.")

            # Parse the input PDB
            st = pr.parsePDB(args.input_pdb)
            if st is None:
                raise ValueError(f"Could not parse {args.input_pdb}")

            # Keep all waters
            for res in st.iterResidues():
                if res.getResname().upper() in WATER_NAMES:
                    water_residues.append((res.getChid(), res.getResnum()))

            # Keep all other HETATM residues that are not waters or the ligand
            for res in st.iterResidues():
                resname = res.getResname().upper()
                if resname in WATER_NAMES:
                    continue
                if args.ligand_resname and resname == args.ligand_resname.upper():
                    continue
                if res.ishetero:
                    cofactor.append(resname)

            args.water_residues = water_residues
            args.cofactor = cofactor
        #adding back in solvent
        if args.water_residues or template_variants or water_residues:
            print(f"→ Adding back selected waters and/or variants")
            with_water_pdb = os.path.join(base_out_dir, "1_2_water_and_variants.pdb")
            add_selected_waters_and_varients(input_pdb_for_fix, fixed_pdb, with_water_pdb, water_residues, template_variants)
            fixed_pdb = with_water_pdb

        # --- ligands & cofactors handling ---
        ligands_to_parametrize, cofactor_resnames = [], []

        # Ligand handling
        if args.ligand_pdb:
            out_sdf = os.path.join(base_out_dir, f"{os.path.splitext(os.path.basename(args.ligand_pdb))[0]}.sdf")
            st = pr.parsePDB(args.ligand_pdb)
            ligand_resname = list({res.getResname() for res in st.iterResidues()})[0]
            make_sdf_from_residue(args.input_pdb, args.ligand_resname, args.ligand_chain,
                                out_sdf, ligand_pdb=args.ligand_pdb)
            ligands_to_parametrize.append((ligand_resname, out_sdf))
            print(f"→ Using ligand pdb: {args.ligand_pdb}")
        elif args.ligand_sdf:
            print(f"→ Using ligand SDF: {args.ligand_sdf}")
            ligands_to_parametrize.append(("LIG", args.ligand_sdf))
        elif args.ligand_resname and args.ligand_chain:
            out_sdf = os.path.join(base_out_dir, f"{args.ligand_resname}_{args.ligand_chain}.sdf")
            make_sdf_from_residue(args.input_pdb, args.ligand_resname, args.ligand_chain, out_sdf, smiles=args.ligand_smiles)
            ligands_to_parametrize.append((args.ligand_resname, out_sdf))
        else:
            print("→ No ligand specified; preparator will be run without explicit ligands.")

        # Cofactor handling
        if args.cofactor:
            for entry in args.cofactor:
                parts = entry.split(":")
                if len(parts) < 2:
                    logging.error("Cofactor must be NAME:VAL[:SMILES]; skipping entry.")
                    continue
                name = parts[0].strip()
                val = parts[1].strip()
                smiles = parts[2].strip() if len(parts) == 3 else None

                if os.path.isfile(val) and val.lower().endswith(".sdf"):
                    print(f"→ Cofactor {name} provided as SDF file: {val}")
                    ligands_to_parametrize.append((name, val))
                    cofactor_resnames.append(name)
                else:
                    # val is chain_id
                    out_sdf = os.path.join(base_out_dir, f"{name}_{val}.sdf")
                    print(f"→ Cofactor {name} provided as residue on chain {val}; generating SDF -> {out_sdf}")
                    make_sdf_from_residue(args.input_pdb, name, val, out_sdf, smiles=smiles)
                    ligands_to_parametrize.append((name, out_sdf))
                    cofactor_resnames.append(name)

        # Decide whether to pass a str or a list into run()
        if not ligands_to_parametrize:
            ligands_arg = None
            print("No ligands or cofactors will be parametrized (continuing with protein only).")
        elif len(ligands_to_parametrize) == 1:
            ligands_arg = ligands_to_parametrize[0][1]  # just the sdf path string
            print(f"\nSingle ligand/cofactor to parametrize: {ligands_to_parametrize[0]}")
        else:
            ligands_arg = ligands_to_parametrize
            print("\nLigands/cofactors to be parametrized and added to the system:")
            for lig_name, lig_path in ligands_to_parametrize:
                print(f"  {lig_name} -> {lig_path}")
            print("The preparator will attempt to parametrize these and add them to the modeller.\n")

        combined_pdb = os.path.join(base_out_dir, "1_3_md_input_no_solvent.pdb")
        if ligands_to_parametrize:
            combine_protein_and_ligands(fixed_pdb, ligands_to_parametrize, combined_pdb)
        # Run SystemPreparation
        print("→ Running SystemPreparation (this may take a while)...")
        preparator = SystemPreparation(
                forcefield=["amber14-all.xml", "amber14/tip3pfb.xml"],
                lig_ff="espaloma",
                allow_undefined_stereo=True,
                hydrogenMass=None,
                boxShape="dodecahedron",
                padding=0.5,
                ionicStrength=0.0,
                is_membrane=not args.soluble_protein,
                lipid_type='POPC' if not args.soluble_protein else None,
                out_dir=md_out_dir
            )
        system, topology = preparator.run(protein=fixed_pdb, ligands=ligands_to_parametrize if ligands_to_parametrize else None)
        logging.info("System preparation complete.")

        # Update equilibration restraints

        write_equilibration_json(
        water_residues_to_keep=water_residues, cofactor_resnames=cofactor_resnames, output_file=os.path.join(md_out_dir, "equilibration.json"))
        # Load system.xml and run Equilibration
        pdb = PDBFile(f'{md_out_dir}/system.pdb')
        topology = pdb.topology
        xml_file = f'{md_out_dir}/system.xml'
        with open(xml_file, "r") as f:
            system = XmlSerializer.deserialize(f.read())

        equilibration = Equilibration(
            topology=topology,
            system=system,
            out_dir=md_out_dir,
            restrained_minimization=True,
            protocol_fname=os.path.join(md_out_dir, "equilibration.json"),
            timestep=0.002,
            is_membrane=True,
            verbose=2
        )
        equilibration.run(pdb_file=f'{md_out_dir}/system.pdb', run_id="sys")
        logging.info("Step 1 equilibration complete.")

        clean_pdb(f"{md_out_dir}/sys_minim.pdb", f"{base_out_dir}/1_4_minimized.pdb", water_residues)
        clean_pdb(f"{md_out_dir}/sys_equilibrated.pdb", f"{base_out_dir}/1_5_final.pdb", water_residues)


if __name__ == "__main__":
    main()
