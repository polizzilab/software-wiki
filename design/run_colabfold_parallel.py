'''
Example script for running colabfold in parallel using python's multiprocessing.Pool on our workstation's GPUs

Input:
    - a directory with .fasta files
Parameters: 
    - number of GPUs to run on
    - and all the colabfold parameters --- see `colabfold_batch --help`
Output:
    - a directory of directories containing pdb files
        each subdirectory contains the results of one .fasta

- jchang 2023-12-06
'''
from pathlib import Path
import subprocess
import multiprocessing
from tqdm import tqdm

from run_colabfold import run_colabfold, COLABFOLD_EXEC

queue = multiprocessing.Queue()

def run_colabfold_parallel(input_dir_of_fastas : str,
                           output_dir : str,
                           flags_list : list = [],
                           num_gpus : int = 8,
                           colabfold_exec : str = COLABFOLD_EXEC,
                           verbose : int = 0):
    '''
    Parameters:
        input_dir_of_fastas 
        output_dir          
        flags_list          list of parameters to pass into colabfold, e.g.
                            [ '--msa-mode', 'single_sequence', '--overwrite-existing-results']

    Optional parameters:
        num_gpus            how many GPUs to run on
        colabfold_exec      what executable to run
        verbose             How much output to show
                            0   just show a tqdm progress bar
                            1   show which .fasta files are sent to which gpus
                            2   show outputs from all commands

    Returns: 
        None
    '''
    fastas = [p for p in Path(input_dir_of_fastas).glob('*') 
              if p.suffix in {'.fa', '.fasta', '.a3m'} ]
    fastas.sort()
    output_subdirs = [Path(output_dir)/f.stem for f in fastas]

    args_list = [
        (fasta, output_subdir, flags_list, colabfold_exec, verbose)
        for (fasta, output_subdir) in zip(fastas, output_subdirs)
    ]

    # Distribute jobs over all GPUs
    num_proc_per_gpu = 1
    for gid in range(num_gpus):
        for _ in range(num_proc_per_gpu):
            queue.put(gid)
    num_processes = num_gpus * num_proc_per_gpu


    with multiprocessing.Pool(processes=num_processes) as p:
        iterable = p.imap(_worker_run_colabfold, args_list)
        if verbose == 0:
            # show the progress bar
            iterable = tqdm(iterable, total=len(args_list))
        for _ in iterable:
            pass

def _worker_run_colabfold(tup):
    input_fasta, output_dir, flags_list, colabfold_exec, verbose = tup

    # Get the next available gpu from the queue
    gpu_id = queue.get()
    try:
        if verbose > 0:
            import time
            print(time.ctime(), multiprocessing.current_process().name, 'on GPU', gpu_id, ': folding', Path(input_fasta).name)

        run_colabfold(input_fasta, output_dir, flags_list, gpu_id=gpu_id, colabfold_exec=colabfold_exec, verbose=verbose)
    finally:
        # Return gpu to queue when finished
        queue.put(gpu_id)

if __name__=='__main__':
    run_colabfold_parallel('colabfold_input/clps_fastas/',
                           'colabfold_output/clps/pdb/',
                           ['--msa-mode', 'single_sequence'],
                           num_gpus = 8,
                           verbose=1)

