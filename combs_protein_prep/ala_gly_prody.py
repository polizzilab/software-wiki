'''
Functions to use prody instead of rosetta flexbb to convert proteins to all glycine or all alanine. This can be more convenient than using fixbb especially when in a python script.

Usage:

import prody as pr
from ala_gly_prody import convert_to_pdb_gly, convert_to_pdb_ala

pdb = pr.parsePDB('ABCD.pdb')
pdb_gly = convert_to_pdb_gly(pdb)
pdb_ala = convert_to_pdb_ala(pdb_gly)
'''

import numpy as np
import prody as pr 

def convert_to_pdb_gly(pdb):
    '''
    Convert any protein to all glycine using prody, simply by taking backbone atoms (N CA C O) and setting resnames to GLY.
    Gets rid of non-protein entities like ligands and waters.
    Also sets Segment names to 'A' but keeps existing Chain IDs.
    '''
    pdb = pr.parsePDB(pdb)
    pdb = pdb.protein # get rid of ligands and waters
    pdb_gly = pdb.select('name N CA C O').toAtomGroup()
    pdb_gly.setResnames(['GLY']*len(pdb_gly))

    return pdb_gly


def convert_to_pdb_ala(pdb_gly):
    '''
    Convert a prody atomgroup of glycine residues to a prody object of all alanines.
    Also sets Segment names to 'A' but keeps existing Chain IDs.
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