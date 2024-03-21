import prody as pr
import os
import sys
COMBS_DIR='/nfs/polizzi/shared/programs/design/Combs2/'
sys.path.append(COMBS_DIR)
import numpy as np

from combs2.design.convex_hull import partition_res_by_burial
from collections import defaultdict 


def convert_to_pdb_ala(pdb_gly):
    '''
    convert a prody atomgroup of glycine residues to a prody object of all alanines
    pdb_gly can be created from any prody pdb using p.select('name N CA C O and not chain W').toAtomGroup()
    '''
    
    N = pdb_gly.select('name N').getCoords()
    CA = pdb_gly.select('name CA').getCoords()
    C = pdb_gly.select('name C').getCoords()
    O = pdb_gly.select('name O').getCoords()
    b = CA - N
    c = C - CA
    a = np.cross(b, c, axis=-1)
    CB = -0.58273431*a + 0.56802827*b - 0.54067466*c + CA
    bb_coords = np.stack([N, CA, C, O, CB], axis=1)
    bb_coords = bb_coords.reshape(-1, 3)
    
    # convert to float32
    pdb_ala = pr.AtomGroup()
    pdb_ala.setCoords(bb_coords)
    pdb_ala.setResnums(np.repeat(pdb_gly.ca.getResnums(),5))
    pdb_ala.setChids(np.repeat(pdb_gly.ca.getChids(),5))
    pdb_ala.setNames(['N','CA','C','O','CB']*len(pdb_gly.ca))
    pdb_ala.setElements(['N','C','C','O','C']*len(pdb_gly.ca))
    pdb_ala.setResnames(['ALA']*len(pdb_gly.ca)*5)
    pdb_ala.setSegnames(['A']*len(pdb_gly.ca)*5)

    return pdb_ala

# dictionary describing which residues are preferred for each burial state
BURIAL_RES_DICT={
    'exposed':'ARNDQEGHKPST',
    'buried':'AFGILMPSTVWY', 
    'intermediate':'ANGHILMFPSTWYV'
}

# dictionary describing which residues are preferred for each charge state (used for topology charge dict)
POLAR='TSHQNA'
POS='RK'
NEG='DE'
NONCHARGED = 'ACFGHILMNPQSTVWY'

CHARGE_RES_DICT = {
    int(-1): list(NEG+POLAR), # remove pos charged - K and R
    int(0): list(NONCHARGED),
    int(1): list(POS+POLAR)
}


def write_resfile(p: pr.AtomGroup,
                      alpha: float = 9, # adjust the size of alpha hull for burial calculation
                      distance_threshold: float = -1.0, # distance from alpha hull to residue. if the distance is greater than the threshold, the residue is considered exposed
                      topo_charge_dict_path: str = None,
                      designable_residues_path: str = None,
                      vdm_residues_path : str = None,
                      resfile_output_path: str = None,
                      vdm_rosetta_rotamer: str = 'NATAA', 
                      ligand_chain : str = 'L',
                      ligand_resnum : str = 10 ): # either NATRO or NATAA
    
    '''

    Inputs:
        p : AtomGroup of the protein
        alpha : Size of alpha hull for burial calculation
        distance_threshold : Distance from alpha hull to residue. If the distance is greater than the threshold, the residue is considered exposed
        topo_charge_dict_path : Path to a file that contains the calculated charge preferences for each resindex.
        designable_residues_path : Path to a file that contains the designable residues. Each line should contain the chain and resnum of the designable residue e.g. A 10
        vdm_residues_path : Path to a file that contains the vdM residues. Each line should contain the chain and resnum of the vdM residue e.g. A 10
        
    Intermediates:
        resfile_dict : 
            Dictionary that holds the preferred residues for each resindex.
            The keys are resindices and the values are lists, where each item is a string describing the residue preferences.
            e.g. {1:['ARNDQEGHKPST', 'DE'], 2:['AFGILMPSTVWY'], }
            Afterwards, the intersection of the strings in the list will be taken to get the final string for the resfile.

    Outputs:
        resfile_output_path : Path to write the resfile to

    '''


    resfile_dict = defaultdict(list) 


    # use combs partition_by_burial fxn which uses alpha hull to partition residues into exposed, intermediate, and buried
    # we need an alanine pdb to calculate the alpha hull - convert_to_pdb_ala will convert any input PDB (glycine or full atom) into all alanines
    pdb_ala = convert_to_pdb_ala(p.select('name N CA C O and (not chain W and not chain L)'))
    resindices_exposed, resindices_intermediate, resindices_buried = partition_res_by_burial(pdb_ala, alpha=alpha, ahull_ca=None, ahull_cb=None, 
                                                                                            assign_intermediate_by_distance=False,
                                                                                            distance_threshold=distance_threshold,)
    
    # for each resindex, add preferred residues for each burial state
    for r in resindices_exposed:
        resfile_dict[r].append(BURIAL_RES_DICT['exposed'])
    for r in resindices_intermediate:
        resfile_dict[r].append(BURIAL_RES_DICT['intermediate'])
    for r in resindices_buried:
        resfile_dict[r].append(BURIAL_RES_DICT['buried'])


    # if topo_charge_dict_path is provided, read it
    if topo_charge_dict_path:
        with open(topo_charge_dict_path,'r') as f: 
            topo_charge_dict = {int(k):int(v) for k, v in [line.split('\t') for line in f.read().split('\n')]}
            # topo_charge_dict is a dict of resindex: (-1, 0, 1) for each resindex in the pdb_ala depending on the charge. not all the resindices are included
            # i.e. {8:-1, 9:0}
        for resindex, charge in topo_charge_dict.items():
            resfile_dict[resindex].append(CHARGE_RES_DICT[charge]) # add the preferred residues for the charge state

    # if designable_residues_path is provided, read it
    if designable_residues_path: 
        designable_residues = []
        with open(designable_residues_path, 'r') as f:
            for line in f:
                chain, resnum = line.split()
                designable_residues.append((chain, int(resnum)))
    else:
        # all the resnums are designable 
        designable_residues = [(res.getChid(), res.getResnum()) for res in p.iterResidues() if res.getChid() != ligand_chain]

    # if vdms_residues_path is providedo, read it
    if vdm_residues_path:
        vdm_residues = []
        with open(vdm_residues_path, 'r') as f:
            for line in f:
                chain, resnum = line.split()
                vdm_residues.append((chain, int(resnum))) 

    # write the resfile
    with open(resfile_output_path, 'w') as f:
        f.write(f'{vdm_rosetta_rotamer} \n') 
        f.write('start \n')
        f.write(f'{ligand_resnum} {ligand_chain} NATRO \n') # don't move the ligand 
        for res in p.iterResidues():
            if res.getChid() == ligand_chain:
                continue
            if (res.getChid(), res.getResnum()) in vdm_residues: # skip the vdM residues
                continue 
            if (res.getChid(), res.getResnum()) in designable_residues:
                # get the preferred residues for this resindex
                preferred_res = resfile_dict[res.getResindex()]
                # get the intersection of the two strings in preferred_res list
                preferred_res = ''.join(list(set.intersection(*map(set, preferred_res))))
                f.write(f'{res.getResnum()} {res.getChid()} PIKAA {preferred_res} \n')
            else:
                f.write(f'{res.getResnum()} {res.getChid()} NATAA \n') 
                # if not a vdM or designable residue, don't change from existing residue? TODO: allow for natural sequence input and PIKAA 
