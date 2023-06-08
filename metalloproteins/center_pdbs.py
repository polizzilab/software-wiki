import prody as pr
import numpy as np

dataset = '/nfs/polizzi/jvaldiviezo/ca-fr/ls.txt'
#dataset = '/nfs/polizzi/jvaldiviezo/ca-frag/list.txt'

# Open the text file in read mode
with open(dataset, 'r') as file:
    # Read all lines from the file
    lines = file.readlines()

# Iterate over the lines
for line in lines:
    # Remove leading/trailing whitespace and newline characters
    pdb_path = line.strip()
    print('pdb_path:', pdb_path)
    # Load the PDB file
    p = pr.parsePDB(pdb_path)
    print('structure:', p)
    # Get the PDB name
    pdb_name = p.getTitle()

    print(pdb_name)


    cal = p.select('hetero')
    coords = cal.getCoords()[0,:] #1,3

    oriCoords = p.getCoords()
    cal_coords = np.zeros(oriCoords.shape)
    N,_ = oriCoords.shape
    for n in range(N):
        cal_coords[n,:] = coords


    newCoords = oriCoords - cal_coords
    p.setCoords(newCoords)

    pr.writePDB(f'/nfs/polizzi/jvaldiviezo/ca-fr/fragments_clustered_aligned/{pdb_name}_aligned.pdb',p)





