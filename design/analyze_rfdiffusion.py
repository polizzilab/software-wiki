'''
Analyze RFDiffusion outputs

Input:
    a single .pdb (and .trb) of RFDiffusion generated scaffold
    or entire directory of .pdb and .trb files
Optional inputs:
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
from collections import defaultdict


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
            {'contig': 'E8-10', 'resnum_start': 20, 'resnum_end': 22}, 
            {'contig': 'E30-33', 'resnum_start': 49, 'resnum_end': 52}, 
            {'contig': 'E78-81', 'resnum_start': 90, 'resnum_end': 93}
        ]
    '''
    curr = 1
    pieces = []
    for part in contig_str.split('/'):
        if part == '0':
            # this is the resnum_end
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
                'resnum_start' : curr,
                'resnum_end' : curr + length - 1,
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

def extract_info_from_trb(trb_filename):
    '''Returns dict of relevant information from the .trb file output by RFDiffusion'''
    trb = np.load(trb_filename, allow_pickle=True)

    # plddt comes as a (num_time_steps, seq_length) numpy array;
    # we want the last timestep here
    plddt = trb['plddt'][-1,:]
    mask_motif = trb['mask_1d']

    # sampled_mask contains the spacer lengths actually sampled
    # we want to convert this into actual residue numbers
    # example: ['19-19/E8-10/26-26/E30-33/37-37/E78-81/17-17/0', 'G0-5/0']
    sampled_mask = trb['sampled_mask']
    chains = dict()
    resnum_start = 1
    for i,contig_str in enumerate(sampled_mask):
        pieces, length = parse_contig_string(contig_str)
        chain = 'ABCDEFGH'[i]
        chains[chain] = {
            'contig' : contig_str,
            'resnum_start' : resnum_start,
            'resnum_end' : resnum_start + length - 1,
            'length' : length,
            'pieces' : pieces,
        }
        resnum_start += length

    return {
        'overall_plddt' : float(plddt.mean()),
        'motif_plddt' : float(plddt[mask_motif].mean()),
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


def calculate_dssp_str(scaffold_path, dssp_filename):
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
    try:
        co = sum(seq_seps) / (len(ca_coords) * len(seq_seps))
    except ZeroDivisionError:
        co = np.nan
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
        log_filename : str = None,
        output_json : str = None,
        designed_chains : list = None
        ):
    '''
    Parses the .pdb and .trb file outputted by RFDiffusion
        (the .trb file is assumed to have the same stem as the .pdb file)

    Parameters:
        scaffold_path           .pdb filepath of the scaffold created by RFDiffusion

    Optional parameters:
        calculate_dssp      whether to run dssp (default: True)
        log_filename        if provided, extra fields are also parsed from the .log file of RFDiffusion output
                                (make sure you provide the correct log file!)
        output_json         if provided, the parsed rfdiffusion output is written here as a .json
        designed_chains     which chains were designed, e.g. ['A'] (default: all chains)

    Returns: 
        dict with the following fields
            scaffold_path
            scaffold_stem
            overall_plddt
            motif_plddt
            seq                     (only if log_filename is provided)
            motif_rmsd              (only if log_filename is provided)
            chains
                contig
                resnum_start
                resnum_end
                length
                dssp_str            (only for designed chains)
                contact_order       (only for designed chains)
                radius_of_gyration  (only for designed chains)
    '''
    trb_file = str(scaffold_path).replace('.pdb', '.trb')
    dssp_filename = str(scaffold_path).replace('.pdb', '.dssp')

    info_dict = dict()
    stem = Path(scaffold_path).stem
    info_dict['scaffold_stem'] = stem
    info_dict['scaffold_path'] = str(Path(scaffold_path).resolve())

    # Parse trb file
    info_dict.update(extract_info_from_trb(trb_file))

    # Parse log file
    if log_filename is not None:
        info_dict.update(extract_info_from_log(log_filename).get(stem, dict()))

    scaffold = pr.parsePDB(scaffold_path)
    if designed_chains is None:
        designed_chains = list(sorted(set(scaffold.getChids())))

    if calculate_dssp:
        dssp_dict = calculate_dssp_str(scaffold_path, dssp_filename)
    else:
        dssp_dict = dict()

    for chain in designed_chains:
        scaffold_chain_ca = scaffold.select(f'chain {chain} name CA')
        if scaffold_chain_ca is None or len(scaffold_chain_ca) == 0:
            raise Exception(f'Chain {chain} not found in {scaffold_path}')
        ca_coords = scaffold_chain_ca.getCoords()

        chain_dict = info_dict['chains'][chain]
        chain_dict['dssp_str'] = dssp_dict.get(chain, '')
        chain_dict['contact_order'] = calculate_contact_order(ca_coords)
        chain_dict['radius_of_gyration'] = calculate_radius_of_gyration(ca_coords)

    if output_json is not None:
        Path(output_json).parent.mkdir(exist_ok=True, parents=True)
        with open(output_json, 'w') as f:
            json.dump(info_dict, f, indent=2)

    return info_dict

def analyze_rfdiffusion_outputs_dir(
        scaffolds_dir : str,
        calculate_dssp : bool = True,
        log_filename : str = None,
        output_json : str = None,
        designed_chains : list = None
        ):
    '''
    This is same as `analyze_rfdiffusion_output` except it takes an entire *directory* of
    .pdb and .trb files as an input. If you pass a .log file, please make sure it is the .log
    file outputted by the RFDiffusion script when it generated that entire directory.
    '''
    if not Path(scaffolds_dir).is_dir():
        raise Exception(f'scaffolds_dir should be a directory: {scaffolds_dir}')

    scaffold_paths = [str(p) for p in Path(scaffolds_dir).glob('*.pdb')]
    scaffold_paths.sort()
    if len(scaffold_paths) == 0:
        raise Exception(f'no pdb files found in directory: {scaffolds_dir}')
    
    # Look at each of the scaffold pdbs in the dir
    out = []
    for scaffold_path in scaffold_paths:
        out.append(analyze_rfdiffusion_output(scaffold_path, 
                                              calculate_dssp=calculate_dssp, 
                                              log_filename=log_filename, 
                                              designed_chains=designed_chains))

    # Add info from log file if present
    if log_filename is not None:
        log_dict = extract_info_from_log(log_filename)
        for d in out:
            d.update(log_dict.get(d['scaffold_stem'], dict()))

    if output_json is not None:
        Path(output_json).parent.mkdir(exist_ok=True, parents=True)
        with open(output_json, 'w') as f:
            json.dump(out, f, indent=2)

    return out

if __name__=='__main__':
    # Analyze a single .pdb and .trb file
    analyze_rfdiffusion_output(
        'rfdiffusion_output/pdb/3dnj_0.pdb',
        log_filename='rfdiffusion_output/rfdiff.log',
        output_json='rfdiffusion_output/3dnj_0.json',
        designed_chains=['A']
    )

    # Analyze an entire directory of .pdb and .trb files
    analyze_rfdiffusion_outputs_dir(
        'rfdiffusion_output/pdb/',
        log_filename='rfdiffusion_output/rfdiff.log',
        output_json='rfdiffusion_output/3dnj.json',
        designed_chains=['A']
    )
