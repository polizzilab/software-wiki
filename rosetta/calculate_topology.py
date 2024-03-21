
from utils.topology import Topology
import numpy as np
import prody as pr


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

def calculate_and_write_topology(p: pr.AtomGroup,
                                 calpha_distance: float = 12,
                                 alpha: float = 20,
                                 num_iterations: int = 1000000,
                                 topo_charge_dict_path: str = None,):
    """
    Calculate the topology and write the topology charge dictionary to a file.
    Parameters:
        p (prody.AtomGroup)
            The input pdb.
        calpha_distance (float)
            TODO
        alpha (float)
            TODO
        num_iterations (int)
            The number of iterations for the monte carlo simulation.
        topo_charge_dict_path (str)
            The path to write the topology charge dictionary to.
    """


    print('calculating topology...')
    pdb_ala = convert_to_pdb_ala(p)
    
    top = Topology()
    top.load_pdb(pdb_ala, selection=None)
    top.load_pdb_ala(pdb_ala, selection=None)
    top.set_topologies(outdir=None)
    top.set_surface_res(alpha=alpha) # 20
    top.set_contacts(calpha_distance=calpha_distance) # 12
    top.run_mc(num_iterations=num_iterations) # 1000000
    seqs = [top.seqs[i] for i in range(len(top.seqs)) if np.sum(list(top.seqs[i].values())) <= 0]
    en_gaps = [top.en_gaps[i] for i in range(len(top.seqs)) if np.sum(list(top.seqs[i].values())) <= 0]
    en_fs = [top.en_fs[i] for i in range(len(top.seqs)) if np.sum(list(top.seqs[i].values())) <= 0]
    seq_rep_lens = [top.seq_rep_lens[i] for i in range(len(top.seqs)) if np.sum(list(top.seqs[i].values())) <= 0 ]
    top.seq_rep_lens = seq_rep_lens
    top.en_gaps = en_gaps
    top.en_fs = en_fs
    top.seqs = seqs
    top.find_pareto_front()
    top.find_nearest_utopian_pt(weight_en_f=0.75, weight_seq_rep_len=0.5)
    top.map_seq_resnums_to_pdb(p)
    top.set_charge_groups(top.seqs[top.nearest_utopian_pt])

    seq_dict = top.seqs[top.nearest_utopian_pt]
    topo_charge_dict = {top.resnum_conv[k]:v for k, v in seq_dict.items()} # renumber the seq_dict to match the pdb_ala

    with open(topo_charge_dict_path, 'w') as f:
        f.write('\n'.join([f'{k}\t{v}' for k,v in topo_charge_dict.items()]))
