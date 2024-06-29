'''
Example script for analyzing rosettafold all atom outputs


'''

import prody as pr 
import os 
import json
import numpy as np
from collections import defaultdict
import pandas as pd
from tqdm import tqdm
import torch # need torch to open the .pt files...
from packing_density import calc_packing_density

dataset_atom_order = {
    'G': [[]],
    'X': [['N', 'CA', 'C', 'O']],
    'A': [['CB']],
    'S': [['CB', 'OG']],
    'C': [['CB', 'SG']],
    'T': [['CB', 'OG1', 'CG2']],
    'P': [['CB', 'CG', 'CD']],
    'V': [['CB', 'CG1', 'CG2'], ['CB', 'CG2', 'CG1']],
    'M': [['CB', 'CG', 'SD', 'CE']],
    'N': [['CB', 'CG', 'OD1', 'ND2'], ['CB', 'CG', 'ND2', 'OD1']],
    'I': [['CB', 'CG1', 'CG2', 'CD1']], 
    'L': [['CB', 'CG', 'CD1', 'CD2'], ['CB', 'CG', 'CD2', 'CD1']], # flip 
    'D': [['CB', 'CG', 'OD1', 'OD2'], ['CB', 'CG', 'OD2', 'OD1']],
    'E': [['CB', 'CG', 'CD', 'OE1', 'OE2'], ['CB', 'CG', 'CD', 'OE2', 'OE1']],
    'K': [['CB', 'CG', 'CD', 'CE', 'NZ']],
    'Q': [['CB', 'CG', 'CD', 'OE1', 'NE2'], ['CB', 'CG', 'CD', 'NE2', 'OE1']],
    'H': [['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'], ['CB', 'CG', 'CD2', 'ND1', 'NE2', 'CE1']],
    'F': [['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'], ['CB', 'CG', 'CD2', 'CD1', 'CE2', 'CE1', 'CZ']], # flip
    'R': [['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2']],
    'Y': [['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'], ['CB', 'CG', 'CD2', 'CD1', 'CE2', 'CE1', 'CZ', 'OH']], # flip
    'W': [['CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'NE1', 'CZ2', 'CZ3', 'CH2']] 
}

three2one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

one2three = {v: k for k, v in three2one.items()}

bb_atoms = ['N','CA','C','O','HA','H']


def compute_RMSD(native_protein, designed_protein):
    """
    Compute the RMSD between two proteins using the alpha carbons.
    Inputs:
        native_protein: prody protein object
        designed_protein: prody protein object
    Outputs:
        calpha_rmsd: float, RMSD between the two proteins
        designed_protein: prody protein object, re-aligned by the native protein
    """
    pr.calcTransformation(designed_protein.ca.getCoords(), native_protein.ca.getCoords()).apply(designed_protein)
    calpha_rmsd = pr.calcRMSD(designed_protein.ca, native_protein.ca)
    return calpha_rmsd, designed_protein

def compute_ligand_RMSD(rfaa_protein, designed_protein, ligand_chain : str ='L'):
    '''
    Compute the RMSD between the ligand atoms of a designed protein and the same atoms of the folded protein.
    Inputs:
        native_protein: prody protein object, folded protein
        designed_protein: prody protein object, designed protein
        ligand_chain: str, chain of the ligand in the designed protein
    Outputs:
        rmsd: float, RMSD between the ligand heavy atoms of the designed and folded proteins
    '''
    ligand_heavy_atoms = designed_protein.select(f'chain {ligand_chain} and not element H')
    designed_ligand_coords = ligand_heavy_atoms.getCoords()
    rfaa_ligand_atoms = rfaa_protein.select(f'resname LG1 and not element H')
    rfaa_ligand_coords = rfaa_ligand_atoms.getCoords()

    # assert that the ordering of the elements is the same, if not print something
    assert np.array_equal(ligand_heavy_atoms.getElements(), rfaa_ligand_atoms.getElements()), "The ordering of ligand atoms may not be same between the designed and folded proteins! This will result in an incorrect ligand RMSD"
    rmsd = pr.calcRMSD(designed_ligand_coords, rfaa_ligand_coords)

    return rmsd


    

def compute_binding_shell_RMSD(pdb, 
                               af2_pdb, 
                               ligand_chain : str ='L',
                               water_chain : str ='W',
                               calc_rmsd_before_alignment: bool = False):
    '''
    Compute the RMSD between the N CA C O atoms of the binding site atoms of a designed protein and the same atoms of the folded protein.
    The binding site is calculated as all residues with heavy atoms within 5A of the ligand.
    Also returns the protein re-aligned by the binding site residues.

    Inputs:
        pdb: prody protein object, designed protein
        af2_pdb: prody protein object, folded protein
        ligand_sele: str, prody selection string for the ligand
        rmsd_before_alignment: bool, whether to return the rmsd before alignment on binding shell residues; turn on if you have already globally aligned by CA and want the rmsd before alignment. 
    Outputs:
        rmsd: float, RMSD between the binding site atoms of the designed and folded proteins
        af2_pdb: prody protein object, re-aligned by the binding site residues
    '''
    ligand_heavy_atoms = pdb.select(f'chain {ligand_chain} and not element H') # typically 'chain L and not element H'
    bs_sele = pdb.select(f'not element H and not chain {water_chain} and not chain {ligand_chain} and within 5 of lig', lig=ligand_heavy_atoms) # don't select waters which are chain W
    # take a set of all the resnums of the binding set and make a prody selection string of these resnums
    bs_resnums = ' '.join([str(i) for i in set(bs_sele.getResnums())])
    # get coordinates of N CA C O atoms of binding site
    bs_coords = []
    bs_coords_af2 = []
    for name in ['N','CA','C','O']:
        bs_coords.append(pdb.select(f"not chain {water_chain} and (resnum {bs_resnums}) and name {name}").getCoords()) # shape N, 3
        bs_coords_af2.append(af2_pdb.select(f"not chain {water_chain} and (resnum {bs_resnums}) and name {name}").getCoords()) # shape N, 3
        # alternatively, could use selection "protein and (resnum {bs_resnums}) and name {name}"
    bs_coords = np.array(bs_coords).reshape(-1,3) # shape 4N, 3
    bs_coords_af2 = np.array(bs_coords_af2).reshape(-1,3) # shape 4N, 3
    # calculate rmsd before alignment
    rmsd_before_alignment = pr.calcRMSD(bs_coords, bs_coords_af2)

    # re-align protein by binding site residues
    pr.calcTransformation(bs_coords_af2, bs_coords).apply(af2_pdb)

    # calculate rmsd
    rmsd = pr.calcRMSD(bs_coords, bs_coords_af2)

    if calc_rmsd_before_alignment:
        return rmsd_before_alignment, af2_pdb

    else:
        return rmsd, af2_pdb


def compute_first_shell_RMSD(pdb, 
                             af2_pdb,
                             ligand_chain : str ='L',
                             water_chain : str ='W',
                             include_backbone_contacts: bool = False,
                             return_first_shell_df: bool = False):
    '''
    Compute the RMSD between the first shell of a designed protein and the same atoms of the folded protein.
    For backbone contacts, (if included) the RMSD of the N CA C O atoms is taken. 
    For sidechain contacts, the RMSD of the sidechain atoms are taken, 
    and the sidechain is flipped if necessary to obtain the lowest RMSD possible.

    Does not do any superposition. Typically this should be calculated after superposing the proteins by the binding site residues.

    Inputs:
        pdb: prody protein object, designed protein
        af2_pdb: prody protein object, folded protein
        ligand_chain: str, chain of the ligand
        water_chain: str, chain of the water
        include_backbone_contacts: bool, whether to include backbone contacts in the RMSD calculation
        return_first_shell_df: bool, whether to return the first shell dataframe

    Outputs:
        rmsd: float, RMSD between the first shell atoms of the designed and folded proteins
        first_shell: dataframe, first shell atoms of the designed protein
            cols: name, resnum, bb
    '''
    
    #af2_pdb = pr.parsePDB(str(path_to_af2_pdb_file)).select('not chain L') # get rid of waters and ligands 
    #pdb = pr.parsePDB(str(laser_model_ouput_path)).select('not chain L')

    # get first shells from probe 
    ligand_heavy_atoms = pdb.select(f'chain {ligand_chain} and not element H')
    first_shell_sel = pdb.select(f'not element H and not chain {water_chain} and not chain {ligand_chain} and within 5 of lig', lig=ligand_heavy_atoms) # don't select waters
    # turn into dataframe
    first_shell = pd.DataFrame(first_shell_sel.getResnums(), first_shell_sel.getNames()).reset_index()
    first_shell.rename(columns={'index':'name', 0:'resnum'}, inplace=True)
    first_shell['bb'] = first_shell['name'].isin(bb_atoms) # check if the atom is a backbone atom

    first_shell_resnums_bb_contact_bool = first_shell[['resnum', 'bb']].drop_duplicates().values

    # get first shell atoms from af2 pdb
    pdb_fs_coords = []
    af2_fs_coords = [] 
    for resnum, bb_contact_bool in first_shell_resnums_bb_contact_bool:
        if include_backbone_contacts:
            if bb_contact_bool:
                # only select the bb atoms 
                for name in ['N','CA','C','O']:
                    pdb_fs_coords.append(pdb.select(f"resnum {resnum} and name {name}").getCoords()[0])
                    af2_fs_coords.append(af2_pdb.select(f"resnum {resnum} and name {name}").getCoords()[0])
        if not bb_contact_bool:
            # select the side chain atoms
            resname = pdb.select(f"resnum {resnum}").getResnames()[0]
            resname_1 = three2one[resname]
            sc_order = dataset_atom_order[resname_1]

            if len(sc_order) > 1: # need to check both flipped sidechain conformations
    
                for atom1, atom2 in zip(sc_order[0], sc_order[0]): # nonflipped
                    #print(atom1, atom2)
                    temp_pdb_fs_coords = []
                    temp_af2_fs_coords = []
                    temp_pdb_fs_coords.append(pdb.select(f"resnum {resnum} and name {atom1}").getCoords()[0])
                    temp_af2_fs_coords.append(af2_pdb.select(f"resnum {resnum} and name {atom2}").getCoords()[0])
                    temp_pdb_fs_coords = np.array(temp_pdb_fs_coords)
                    temp_af2_fs_coords = np.array(temp_af2_fs_coords)
                    rmsd1 = pr.calcRMSD(temp_pdb_fs_coords, temp_af2_fs_coords)

                for atom1, atom2 in zip(sc_order[0], sc_order[1]): # flipped
                    #print(atom1, atom2)
                    temp_pdb_fs_coords = []
                    temp_af2_fs_coords = []
                    temp_pdb_fs_coords.append(pdb.select(f"resnum {resnum} and name {atom1}").getCoords()[0])
                    temp_af2_fs_coords.append(af2_pdb.select(f"resnum {resnum} and name {atom2}").getCoords()[0])
                    temp_pdb_fs_coords = np.array(temp_pdb_fs_coords)
                    temp_af2_fs_coords = np.array(temp_af2_fs_coords)
                    rmsd2 = pr.calcRMSD(temp_pdb_fs_coords, temp_af2_fs_coords)
                
                # calculate rmsd of both and choose the one with the lower rmsd
                sc_order = sc_order[np.argmin([rmsd1, rmsd2])]
            else: # the sidechain can be aligned without trying both flipped conformations
                sc_order = sc_order[0]

            for atom in sc_order: 
                pdb_fs_coords.append(pdb.select(f"resnum {resnum} and name {atom}").getCoords()[0])
                af2_fs_coords.append(af2_pdb.select(f"resnum {resnum} and name {atom}").getCoords()[0])
                
    pdb_fs_coords = np.array(pdb_fs_coords)
    af2_fs_coords = np.array(af2_fs_coords)

    # calculate rmsd of first shell atoms
    rmsd = pr.calcRMSD(pdb_fs_coords, af2_fs_coords)
     
    if return_first_shell_df:
        return rmsd, first_shell
    else:
        return rmsd

def process_rfaa(
        rfaa_dir : str,
        designed_dir : str,
        to_fold : list = None,
        compute_sc_RMSD : bool = True,
        compute_bs_RMSD : bool = True,
        eval_pack_density=False
    ):

    if to_fold is not None:
        print(f'Taking subset of {len(to_fold)} pdbs provided in to_fold list.')
        designed_proteins = sorted([x for x in os.listdir(designed_dir) if x.endswith('.pdb') and x[:-4] in to_fold])
    else:
        designed_proteins = sorted([x for x in os.listdir(designed_dir) if x.endswith('.pdb')])

    rfaa_folded = sorted([x for x in os.listdir(rfaa_dir) if x.endswith('.pdb')])
    rfaa_confidences = sorted([x for x in os.listdir(rfaa_dir) if x.endswith('.pt')])

    assert len(rfaa_folded) == len(rfaa_confidences) == len(designed_proteins)
    print(f'Processing colabfold outputs for {len(rfaa_folded)} folded proteins.')

    df = defaultdict(list)

    for folded, confidence_pt, design in tqdm(zip(rfaa_folded, rfaa_confidences, designed_proteins), total=len(rfaa_confidences)):

        rfaa_pdb_path = os.path.join(rfaa_dir, folded)
        input_pdb_path = os.path.join(designed_dir, design)

        # Read designed protein
        input_protein = pr.parsePDB(input_pdb_path)

        # Read rfaa batch metadata file
        info_path = os.path.join(rfaa_dir, confidence_pt)
        info = torch.load(info_path, map_location="cpu")

        # This assumes the ligand is always at the end of the file
        protein_len = len(input_protein.select('protein').ca.getResnums())
        ligand_plddt = info['plddts'][0][protein_len:].mean().numpy()
        protein_plddt = info['plddts'][0][:protein_len].mean().numpy()
        # do we want PAEs? they are LxL shape

        # Overwrite with aligned backbone.
        rfaa_protein_backbone = pr.parsePDB(rfaa_pdb_path)
        calpha_rmsd, rfaa_protein_backbone = compute_RMSD(input_protein, rfaa_protein_backbone)
        pr.writePDB(rfaa_pdb_path, rfaa_protein_backbone)

        # Calculate binding site RMSD and re-align by binding site residues
        if compute_bs_RMSD:
            bs_rmsd, rfaa_protein_backbone = compute_binding_shell_RMSD(input_protein, rfaa_protein_backbone)

        # Calculate sidechain RMSD after af2 protein backbone has been aligned by binding site residues 
        if compute_sc_RMSD:
            if compute_bs_RMSD == False:
                print("Warning: computing sidechain RMSD without superposing on binding site atoms.")
            sc_rmsd = compute_first_shell_RMSD(input_protein, rfaa_protein_backbone)

        # Get the ligand rmsd.
        ligand_rmsd = compute_ligand_RMSD(rfaa_protein_backbone, input_protein)

        #Calculate packing density
        if eval_pack_density:
            # For speed and to avoid calculating packing density for proteins that fold into spaghetti, only calculate pack density if the rmsd is <2A
            if calpha_rmsd < 2:
                pack_density = calc_packing_density(rfaa_pdb_path, n_workers=1)
            else:
                pack_density = np.nan

    df['pdb_name'].append(input_protein.getTitle())
    df['protein_plddt'].append(protein_plddt)
    df['ligand_plddt'].append(ligand_plddt)
    df['rfaa_rmsd'].append(calpha_rmsd)
    if compute_sc_RMSD:
        df['rfaa_sc_rmsd'].append(sc_rmsd)
    if compute_bs_RMSD:
        df['rfaa_bs_rmsd'].append(bs_rmsd)
    if eval_pack_density:
        df['pack_density'].append(pack_density)
    df['ligand_rmsd'].append(ligand_rmsd)
    df['input_pdb_path'].append(os.path.abspath(input_pdb_path))
    df['rfaa_pdb_path'].append(os.path.abspath(rfaa_pdb_path))
    df['seq'].append(input_protein.protein.ca.getSequence())

    df = pd.DataFrame(df)
    return df
