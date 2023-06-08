import numpy as np
from combs2.design.cluster import Cluster
import prody as pr
import os

#dataset = '/nfs/polizzi/jvaldiviezo/ca-fr/ca_frags_pdb_2p5_0p3.txt'
dataset = '/nfs/polizzi/jvaldiviezo/ca-fr/test.txt'

# Open the text file in read mode
with open(dataset, 'r') as file:
    # Read all lines from the file
    lines = file.readlines()

pdb_frags = []
# Iterate over the lines
for line in lines:
    # Remove leading/trailing whitespace and newline characters
    pdb_path = line.strip()
    # Load the PDB file
    structure = pr.parsePDB(pdb_path)
    # Get the PDB name
    pdb_name = structure.getTitle()
    frags = structure.select('all')
    hetero = structure.select('hetero')

    pdb_frags.append(frags)

# make coord array of fragment backbone atoms
atom_names = ['N', 'CA', 'C', 'O']
pdb_frags_coords = []
pdb_frags_indices = []
for i, frag in enumerate(pdb_frags):
    incomplete_fragment = False
    frag_coords = []
    for resnum in sorted(set(frag.getResnums())):
        res = frag.select(f'resnum {resnum}')
        for atom_name in atom_names:
            atom = res.select(f'name {atom_name}')
            if atom is None:
                frag_coords.append(hetero.getCoords()[0])
                incomplete_fragment = False
            if atom is not None:
                frag_coords.append(atom.getCoords()[0])
            else:
                incomplete_fragment = True
                break
        if incomplete_fragment:
            break
    if len(frag_coords)==53:
        pdb_frags_coords.append(frag_coords)
        pdb_frags_indices.append(i)
pdb_frags_coords = np.array(pdb_frags_coords, dtype=np.float32)
pdb_frags_indices = np.array(pdb_frags_indices)

    
print(pdb_frags_coords.shape)
print(pdb_frags_indices.shape)



# generate rmsd matrix and run the clustering algorithm
cluster = Cluster()
cluster.rmsd_cutoff = 0.5
cluster.pdb_coords = pdb_frags_coords
cluster.make_pairwise_rmsd_mat(superpose=True)
cluster.make_square()
cluster.make_adj_mat()
cluster.fast_cluster(min_cluster_size=2)

print('clusters:', cluster.mems)
print('cluster centroids:', cluster.cents)



outdir_clusters = '/nfs/polizzi/jvaldiviezo/ca-fr/fragments_clustered/'
os.makedirs(outdir_clusters, exist_ok=True)
for i, (centroid, members) in enumerate(zip(cluster.cents, cluster.mems)):
    print('centroid:', centroid)
    print('members:', members)
    print('')
    # can grab the centroid
    centroid_frag = pdb_frags[centroid].copy()
    # superimpose the centroid onto the centroid of largest cluster
    transformation = pr.calcTransformation(pdb_frags_coords[centroid], pdb_frags_coords[cluster.cents[0]])
    centroid_frag_coords = pr.measure.transform._applyTransformation(transformation, pdb_frags_coords[centroid])
    # can superpose the members onto the centroid and then print them
    for j, member in enumerate(members):
        member_frag = pdb_frags[member].copy()
        pr.calcTransformation(pdb_frags_coords[member], centroid_frag_coords).apply(member_frag)
        if member == centroid:
            pr.writePDB(f'{outdir_clusters}cluster_{i}_mem_{j}_centroid.pdb', member_frag)
        else:
            pr.writePDB(f'{outdir_clusters}cluster_{i}_mem_{j}.pdb', member_frag)
