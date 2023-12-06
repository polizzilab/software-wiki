'''
Analyze RFDiffusion outputs

Input:
    .pdb of RFDiffusion generated scaffold
Optional inputs:
    .pdb for any motifs used as input to RFDiffusion
    contig string for the motif
    .log file of the RFDiffusion output
Output:
    dict with information such as radius of gyration, dssp string, etc.

'''
from pathlib import Path
import os
import numpy as np
import prody as pr
from scipy.spatial.distance import cdist
import re
import json


class ContigStringParseException(Exception):
    pass

def parse_contig_part(part: str):
    '''
    Example input:
        E8-10
    Example output:
        ('E', 8, 10)
    '''
    m = re.match(r'([A-Za-z]+)(\d+)-(\d+)', part)
    if not m:
        raise ContigStringParseException(part)
    return m.group(1), int(m.group(2)), int(m.group(3))


def parse_contig_string(contig_str: str):
    '''Returns the starting and ending resnums of fixed pieces in a contig string

    Example input:
        '19-19/E8-10/26-26/E30-33/37-37/E78-81/17-17/0'
    Example output:
        [
            {'contig': 'E8-10', 'start': 20, 'end': 22}, 
            {'contig': 'E30-33', 'start': 49, 'end': 52}, 
            {'contig': 'E78-81', 'start': 90, 'end': 93}
        ]
    '''
    curr = 1
    pieces = []
    for part in contig_str.split('/'):
        if part == '0':
            # this is the end
            break

        if part[0].isdigit():
            # generated part
            length_1, length_2 = part.split('-')
            if length_1 != length_2:
                raise ContigStringParseException(part)
            curr += int(length_1)
        else:
            # fixed motif part
            chain, m_start, m_end = parse_contig_part(part)
            length = m_end - m_start + 1
            pieces.append({
                'contig' : part,
                'start' : curr,
                'end' : curr + length - 1,
            })
            curr += length

    total_length = curr - 1
    return pieces, total_length
    

def extract_info_from_log(log_filename):
    '''Extracts these from the RFDiffusion log printed:
        - Seq
        - Motif RMSD
    '''
    all_info = defaultdict(dict)
    for line in open(log_filename):
        if 'Making design' in line:
            curr_design = Path(line.strip().split()[-1]).stem
        if 'Sequence init' in line:
            all_info[curr_design]['seq'] = line.strip().split()[-1]
        if 'Sampled motif RMSD' in line:
            all_info[curr_design]['motif_rmsd'] = float(line.strip().split()[-1])

    return all_info

def extract_info_from_trb(trb_filename, designed_chains=None):
    '''Returns dict of relevant information from the .trb file output by RFDiffusion'''
    trb = np.load(trb_filename, allow_pickle=True)

    # plddt comes as a (num_time_steps, seq_length) numpy array;
    # we want the last timestep here
    plddt = trb['plddt'][-1,:]
    mask_motif = trb['mask_1d']

    # sampled_mask contains the spacer lengths actually sampled
    # we want to convert this into actual residue numbers
    # example: ['19-19/E8-10/26-26/E30-33/37-37/E78-81/17-17/0', 'G0-5/0']
    sampled_mask = {'ABCDEFGH'[i] : s for i,s in enumerate(trb['sampled_mask']) }
    chains = dict()
    start = 1
    for chain,contig_str in sampled_mask.items():
        pieces, length = parse_contig_string(contig_str)
        if designed_chains is None or chain in designed_chains:
            chains[chain] = {
                'start' : start,
                'end' : start + length - 1,
                'length' : length,
                'pieces' : pieces,
            }
        start += length

    return {
        'overall_plddt' : float(plddt.mean()),
        'motif_plddt' : float(plddt[mask_motif].mean()),
        'sampled_mask' : sampled_mask,
        'chains' : chains,
    }


def parse_dssp(dssp_out_filepath):
    """parses a DSSP file and returns a dictionary where
        - keys are (chain, resnum)
        - values are (dssp code)
    """
    fh = open(dssp_out_filepath, 'r')
    dssp_ss = dict()

    def parse_float(float_str):
        # dssp seems to use 360.0 to represent a nonexistent angle
        if float_str == ' 360.0':
            return np.nan
        else:
            return float(float_str)

    for line in fh:
        if line.startswith('  #  RESIDUE'):
            break
    for line in fh:
        if line[13] == '!': #chain break
            continue
        resnum = int(line[5:10])
        chain = line[11]
        # resname = line[13].upper() # disulflide bond S's are in lowercase
        dssp = line[16].replace(' ', '-')

        dssp_ss[(chain, resnum)] = dssp

    fh.close()
    return dssp_ss


def calculate_dssp_str(scaffold_path, dssp_filename, chains=None):
    '''
    scaffold_path:      pdb to run dssp for
    dssp_filename:      where the dssp output goes
    chains:             which chains to return the dssp string for;
                        if not specified it defaults to all chains

    Returns dict with
        keys = chains
        values = dssp strings
    '''
    if not os.path.exists(dssp_filename):
        os.system(f'dssp {scaffold_path} 2> /dev/null > {dssp_filename}')
    dssp_dict = parse_dssp(dssp_filename)

    # get the dssp str for each chain separately
    if chains is None:
        # if the user does not specify a chain, use all chains
        chains = sorted(set(chain for chain,resnum in dssp_dict.keys()))

    dssp_strs = dict()
    for chain in chains:
        resnums = sorted(resnum for c, resnum in dssp_dict.keys()
                         if c == chain)
        dssp_strs[chain] = ''.join(dssp_dict[chain,resnum] for resnum in resnums)
    return dssp_strs


def calculate_radius_of_gyration(X):
    # I think this is defined as the r.m.s distance to the center of mass
    X = X.copy()
    X -= X.mean(axis=0)
    dev = np.linalg.norm(X, axis=1)
    return dev.mean()

def calculate_contact_order(ca_coords, skip_number=2, contact_distance=8):
    'fxn from Nick'
    dm = cdist(ca_coords, ca_coords)
    i_resinds, j_resinds = np.where(np.triu(dm <= contact_distance, k=1))
    seq_seps = [np.abs(i-j) for i, j in zip(i_resinds, j_resinds) if np.abs(i-j) > skip_number]
    co = sum(seq_seps) / (len(ca_coords) * len(seq_seps))
    return co

def longest_run(s, c):
    'Returns longest run of characters c in string s'
    n = 0
    N = 0
    for k in s:
        if k in c:
            n += 1
            N = max(n, N)
        else:
            n = 0
    return N


def analyze_rfdiffusion_output(
        scaffold_path : str,
        calculate_dssp : bool = True,
        output_json : str = None,
        designed_chains : list = None
        ):
    '''
    Parameters:
        scaffold_path           .pdb filepath of the scaffold created by RFDiffusion

    Optional parameters:
        calculate_dssp      whether to run dssp (default: True)
        output_json         if provided, the parsed rfdiffusion output is written here as a .json
        designed_chains     which chains were designed, e.g. ['A'] (default: all chains)

    Returns: 
        dict with the following fields
            overall_plddt
            motif_plddt
            sampled_mask
    '''
    trb_file = str(scaffold_path).replace('.pdb', '.trb')
    dssp_filename = str(scaffold_path).replace('.pdb', '.dssp')

    scaffold = pr.parsePDB(scaffold_path)
    if designed_chains is None:
        ca_coords = scaffold.select('name CA').getCoords()
    else:
        ca_coords = scaffold.select(f'chain {" ".join(designed_chains)} name CA').getCoords()

    # Parse trb file
    info_dict = extract_info_from_trb(trb_file, designed_chains=designed_chains)
    if calculate_dssp:
        info_dict['dssp_str'] = calculate_dssp_str(scaffold_path, dssp_filename, chains=designed_chains)
    info_dict['contact_order'] = calculate_contact_order(ca_coords)
    info_dict['radius_of_gyration'] = calculate_radius_of_gyration(ca_coords)

    if output_json is not None:
        Path(output_json).parent.mkdir(exist_ok=True, parents=True)
        with open(output_json, 'w') as f:
            json.dump(info_dict, f, indent=2)

    return info_dict

if __name__=='__main__':
    analyze_rfdiffusion_output(
        'rfdiffusion_output/pdb/3dnj_0.pdb',
        output_json='rfdiffusion_output/3dnj_0.json',
        designed_chains='A',
    )
