'''
Example script for running RFdiffusion to generate backbones

Input: a dictionary of flags; see full options at: https://github.com/RosettaCommons/RFdiffusion/tree/main
Output: a directory of pdbs

Notes:
🎵 probe does not know about segment ID's, so make sure your pdb file has unique chains

- jchang 2023-12-05
'''
from pathlib import Path
import os

PYTHON_EXEC = '/nfs/polizzi/jmou/miniconda3/envs/SE3nv/bin/python'
RFDIFFUSION_SCRIPT = '/nfs/polizzi/shared/programs/design/RFdiffusion/scripts/run_inference.py'

def run_rf_diffusion(flags_dict : dict,
                     logfile : str = None,
                     gpu_id : int = None,
                     python_exec : str = PYTHON_EXEC,
                     rfdiffusion_script : str = RFDIFFUSION_SCRIPT,
                     verbose : bool = False):
    '''
    Parameters:
        flags_dict          parameters to pass into rfdiffusion, e.g.
                            {
                                'inference.input_pdb' : 'rfdiffusion_input/3dnj_trunc.pdb',
                                'contigmap.contigs': '[A40-63/21-40/A90-119/0 C1-3]',
                                ...
                            }

    Optional parameters:
        logfile             if provided, the output from rfdiffusion is written to this file
        gpu_id              if provided, runs on the GPU with this id
        python_exec         what executable (conda env) of python to run rfdiffusion with
        rfdiffusion_script  what rfdiffusion script to run
        verbose             whether to show what commmand you're running

    Returns: 
        None
    '''
    flags_list = [f'{key}={value}' for key,value in flags_dict.items()]
    args = [ PYTHON_EXEC, '-u', RFDIFFUSION_SCRIPT] + flags_list
    cmd = ' '.join(repr(a) for a in args)
    if logfile is not None:
        Path(logfile).parent.mkdir(exist_ok=True, parents=True)
        if verbose:
            cmd += f' 2>&1 | tee -a {logfile}'
        else:
            cmd += f' 2>&1 >> {logfile}'
    if verbose:
        print(cmd)
    if gpu_id is not None:
        os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
    os.system(cmd)

if __name__=='__main__':
    input_pdb =  "rfdiffusion_input/3dnj_trunc.pdb"
    output_prefix =  "rfdiffusion_output/pdb/3dnj"
    logfile =  "rfdiffusion_output/rfdiff.log"
    inference_flags = {
        "inference.input_pdb": input_pdb,
        "inference.output_prefix": output_prefix,
        "contigmap.contigs": "[A40-63/21-40/A90-119/0 C1-3]",
        "inference.write_trajectory": False,
        "denoiser.noise_scale_ca": 1,
        "denoiser.noise_scale_frame": 1,
        "inference.num_designs": 3,
        "diffuser.T": 50
    }
    run_rf_diffusion(inference_flags, logfile=logfile, gpu_id=1, verbose=True)
