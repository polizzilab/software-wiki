import prody as pr
import pandas as pd 
import numpy as np 
import itertools


def check_for_hb(a1, a2):
    # don't count hbonds within the same residue
    # sometimes these are for water, so need to allow for these 
    if (a1.getChid(),a1.getResnum())==(a2.getChid(),a2.getResnum()): 
        return False
    # don't count hbonds within backbone atoms
    # backbone atoms are named N, CA, C, O and have the same chid
    names = (a1.getName(), a2.getName())
    chids = (a1.getChid(), a2.getChid())
    if names in list(itertools.product(('O','N','S'), repeat=2)): 
        # if the chains do not contain W, return false 
        if not 'W' in chids:
            return False
    # finally, check if the bond is between O, N 
    l = (a1.getName()[0],a2.getName()[0])
    if l in list(itertools.product(('O','N','S'), repeat=2)):
        return True
    return False 


def make_constraints(
        pdb, # prody object
        prody_selection_str : str,
        max_bond : float = 3.5,
        min_bond : float = 2.5,
        x0 : float = 3,
        tol : float = 0.4,
        sd : float = 0.3,
        x0_h : float = 1.9,
        tol_h : float = 0.3,
        sd_h : float = 0.3,
        output_constraint_file : str = 'constraints.cst',
        return_constraint_list : bool = False, # return a list of constraints: list of tuples (atom1_name, atom1_resnum, atom1_chid, atom2_name, atom2_resnum, atom2_chid)
    ):

    '''
    Generate a Rosetta constraints file to keep hydrogen bonds fixed. Hydrogens must be added.
    Constraints will be written between all hydrogen bonds in the selection. 
    Two constraints are written per hydrogen bond: (Heavy Atom - Hydrogen) and (Heavy Atom - Heavy Atom).

    Automatically generates a FLAT HARMONIC constraint: https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/constraint-file#:~:text=FLAT_HARMONIC%20x0%20sd%20tol 
    From Rosetta Wiki: Zero in the range of x0 - tol to x0 + tol. Harmonic with width parameter sd outside that range. Basically, a HARMONIC potential (see above) split at x0 with a 2*tol length region of zero inserted.

    Constraints are NOT written for:
        - hydrogen bonds within the same residue
        - hydrogen bonds where *both* heavy atoms are backbone atoms (atom names N, CA, C, O) for instance i+4 helix h-bonding
        - hydrogen bonds with S 
        - hydrogen bonds that are missing an H (constraint will only be written if there is a hydrogen attached to one heavy atom within Hbonding distance)

    Inputs:
        pdb : Parsed prody object
        prody_selection_str : Prody selection string to select the atoms to make constraints for.
        max_bond : Maximum bond length for inferring hbonds. Default 3.5
        min_bond : Minimum bond length for inferring hbonds. Default 2.5

        x0 : Mean distance for heavy atom constraints. Default 3
        tol : Tolerance for heavy atom constraints. From x0 - tol to x0 + tol there  is no penalty. Default 0.4 
        sd : Standard deviation for heavy atom constraints. Outside of the tolerance, an energy penalty with this sd is applied. Default 0.3
        
        x0_h : Mean distance for heavy atom - hydrogen constraints. Default 1.9
        tol_h : Tolerance for heavy atom - hydrogen constraints. Default 0.3
        sd_h : Standard deviation for heavy atom - hydrogen constraints. Default 0.3

        output_constraint_file : Path to write the constraints file. Default 'constraints.cst'
        return_constraint_list : Return a list of constraints. Default False. Needed for pymol visualization.
    '''

    vdms_and_lig = pdb.select(prody_selection_str).toAtomGroup()
    vdms_and_lig.inferBonds(max_bond=max_bond, min_bond=min_bond) # infer bonds setBonds=True

    # make a list of constraints later for pymol visualization
    cst_list = []

    with open(output_constraint_file,'w') as f:
        for bond in vdms_and_lig.iterBonds():
            atom1, atom2 = bond.getAtoms()
            if check_for_hb(atom1, atom2):
                # also check for attached hydrogens; if there are any, add those constraints too 
                atom1_h = vdms_and_lig.select('element H and within 1.3 of sel',sel=atom1)
                if atom1_h:
                    dist = pr.calcDistance(atom1_h, atom2)
                    if np.min(dist) < 2.7: # 11/10/2021 changed from 2.5 to 2.7
                        if len(dist) > 1:
                            atom1_h_name = atom1_h.getNames()[np.argmin(dist)]
                            atom1_h = atom1_h.select(f'name {atom1_h_name}')
                        f.write(f'AtomPair {atom1_h.getNames()[0]} {atom1_h.getResnums()[0]}{atom1_h.getChids()[0]} {atom2.getName()} {atom2.getResnum()}{atom2.getChid()} FLAT_HARMONIC {x0_h} {sd_h} {tol_h} \n')
                        cst_list.append((atom1_h.getNames()[0], atom1_h.getResnums()[0], atom1_h.getChids()[0], atom2.getName(), atom2.getResnum(), atom2.getChid())) # add tuple of (atom1_h_name, atom1_h_resnum, atom1_h_chid, atom2_name, atom2_resnum, atom2_chid)
                    else:
                        atom1_h = None

                atom2_h = vdms_and_lig.select('element H and within 1.3 of sel',sel=atom2)
                if atom2_h:
                    dist = pr.calcDistance(atom2_h, atom1) 
                    if np.min(dist) < 2.7: # 11/10/2021 changed from 2.5 to 2.7
                        if len(dist) > 1:
                            atom2_h_name = atom2_h.getNames()[np.argmin(dist)]
                            atom2_h = atom2_h.select(f'name {atom2_h_name}')
                        f.write(f'AtomPair {atom2_h.getNames()[0]} {atom2_h.getResnums()[0]}{atom2_h.getChids()[0]} {atom1.getName()} {atom1.getResnum()}{atom1.getChid()} FLAT_HARMONIC {x0_h} {sd_h} {tol_h} \n')
                        cst_list.append((atom2_h.getNames()[0], atom2_h.getResnums()[0], atom2_h.getChids()[0], atom1.getName(), atom1.getResnum(), atom1.getChid())) # add tuple 
                    else:
                        atom2_h = None
                # if either atom1 or atom2 has a valid hydrogen attached to it, write the heavy atom constraint 
                if atom1_h or atom2_h:
                    f.write(f'AtomPair {atom1.getName()} {atom1.getResnum()}{atom1.getChid()} {atom2.getName()} {atom2.getResnum()}{atom2.getChid()} FLAT_HARMONIC {x0} {sd} {tol} \n')
                    cst_list.append((atom1.getName(), atom1.getResnum(), atom1.getChid(), atom2.getName(), atom2.getResnum(), atom2.getChid()))

    if return_constraint_list:
        return cst_list
    
def pymol_visualize_constraints(
        pdb_path: str,
        cst_list: list,
        output_pse_file: str = 'constraints.pse'
    ):

    from pymol import cmd

    cmd.delete("all")
    cmd.load(pdb_path)
    for c in cst_list:
        cmd.distance(f'{c[0]}_{c[1]}{c[2]}-{c[3]}_{c[4]}{c[5]}', f'resi {c[1]} and name {c[0]} and chain {c[2]}, resi {c[4]} and name {c[3]} and chain {c[5]}')
    # show sticks
    cmd.show_as('sticks')

    # show only polar hydrogens
    cmd.hide('everything', 'ele h')
    cmd.show('lines', 'ele h and neighbor (ele n+o)')

    # save the session
    cmd.save(output_pse_file)
                
def pymol_visualize_constraints_print(
        cst_list: list,
    ):

    for c in cst_list:
        print('distance', f'{c[0]}_{c[1]}{c[2]}-{c[3]}_{c[4]}{c[5]}, resi {c[1]} and name {c[0]} and chain {c[2]}, resi {c[4]} and name {c[3]} and chain {c[5]}')