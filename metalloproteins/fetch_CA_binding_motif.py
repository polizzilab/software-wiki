import prody as pr
import numpy as np

dataset = '/nfs/polizzi/jvaldiviezo/ca-fr/ca_pdb_2p5_0p3.txt'

# Open the text file in read mode
with open(dataset, 'r') as file:
    # Read all lines from the file
    lines = file.readlines()

# Iterate over the lines
for line in lines:
    # Remove leading/trailing whitespace and newline characters
    pdb_path = line.strip()
    # Load the PDB file
    structure = pr.parsePDB(pdb_path)
    # Get the PDB name
    pdb_name = structure.getTitle()
    
    # Iterate over all chains
    chains = np.unique(structure.getChids())
    for chain in chains:

        calcium = structure.select('chain {} and name CA and not protein'.format(chain))
        # if no calcium present, skip chain        
        if calcium is None:
            continue
        calcium.getResnums()
        center = calcium.getCoords()

        protein = structure.select('protein and chain {}'.format(chain))
        # If no protein present, skip chain
        if protein is None:
            continue
        # select residues that are ~3.5 A away from calcium
        BS = protein.select('same residue as within 3.5 of center', center = center)
        # If there are no residues around, skip
        if BS is None:
            continue
        BS_resnum = list(set(BS.getResnums()))

        min_resnum = min(protein.getResnums())
        max_resnum =  max(protein.getResnums())
        interval = 6 # get 13 residues
        for resnum in BS_resnum:
            start_resnum = max(resnum - interval, min_resnum)
            # if index is negative, skip
            if start_resnum < 0:
                continue
            end_resnum = min(resnum + interval, max_resnum)
            frag=structure.select(f'chain {chain} and resnum {start_resnum} to {end_resnum}')
            # include residues and calcium
            motif = structure.select(f'chain {chain} and resnum {start_resnum} to {end_resnum} or name CA and not protein within 3.5 of frag', frag = frag)
            pr.writePDB(f'/nfs/polizzi/jvaldiviezo/ca-fr/fragments/{pdb_name}_{chain}_{resnum}.pdb',motif)
