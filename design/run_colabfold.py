'''
Example script for using colabfold

Input: a .fasta of sequences (or a directory, or a .a3m file, etc.)
Parameters: see the output of `colabfold_batch --help`
Output: a directory of pdbs of predicted structures

- jchang 2023-12-06
'''
from pathlib import Path
import subprocess

COLABFOLD_EXEC = '/nfs/sbgrid/programs/x86_64-linux/colabfold/1.5.2/bin/colabfold_batch'

def run_colabfold(input_file : str,
                  output_dir : str,
                  flags_list : list = [],
                  gpu_id : int = None,
                  colabfold_exec : str = COLABFOLD_EXEC,
                  verbose : int = 0):
    '''
    Parameters:
        input_file          .fasta of sequences (or directory, .a3m file, etc.)
        output_dir          directory for where to write the output
        flags_list          list of parameters to pass into colabfold, e.g.
                            [ '--msa-mode', 'single_sequence', '--overwrite-existing-results']

    Optional parameters:
        gpu_id              if provided, runs on the GPU with this id
        colabfold_exec      what executable to run
        verbose             if > 1 it shows the whole output of the colabfold script

    Returns: 
        None
    '''
    gpu_env_vars = {"CUDA_DEVICE_ORDER": "PCI_BUS_ID"}
    if gpu_id is not None:
        gpu_env_vars["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    cmd_parts = [ str(colabfold_exec), str(input_file), str(output_dir) ] + flags_list

    if verbose > 1:
        print(' '.join(cmd_parts))
        subprocess.run(cmd_parts, env=gpu_env_vars)
    else:
        subprocess.run(cmd_parts, env=gpu_env_vars, stdout=subprocess.DEVNULL)

if __name__=='__main__':
    run_colabfold('colabfold_input/2w9r_seqs.fasta',
                  'colabfold_output/pdb',
                  ['--msa-mode', 'single_sequence'],
                  gpu_id=0,
                  verbose=2)
