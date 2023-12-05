'''
Example script for running reduce to add hydrogens

Input: .pdb structure without hydrogens
Output: .pdb structure with hydrogens
Notes:

🎵 the default SBGrid installation is old; this uses Ben's newly installed version on the workstation
🎵 the "-BUILD" flag suffices for most purposes
🎵 reduce can crash for malformed pdbs, so watch out

- jchang 2023-12-05
'''
import os
from pathlib import Path

def run_reduce(input_pdb : str,
               output_pdb : str,
               overwrite : bool = False,
               verbose : bool = False,
               flags : str = '-BUILD',
               reduce_exec : str ='/nfs/polizzi/bfry/programs/reduce/reduce',
               hetdict_path : str ='/nfs/polizzi/bfry/programs/reduce/reduce_wwPDB_het_dict.txt'):
    if not Path(input_pdb).exists():
        raise FileNotFoundError(input_pdb)
    if Path(output_pdb).exists() and not overwrite:
        raise FileExistsError(output_pdb, '(pass overwrite=True)')
    Path(output_pdb).parent.mkdir(exist_ok=True, parents=True)

    command = f'{reduce_exec} {flags} -DB {hetdict_path} {input_pdb} > {output_pdb} 2> /dev/null'
    if verbose:
        print(command)
    ret = os.system(command)
    # Oddly, reduce can give nonzero exit code even when it succeeds
    # if ret != 0:
    #     raise Exception(f'Reduce failed:\n\t{command}')

if __name__=='__main__':
    run_reduce('reduce_input/5dtx.pdb', 'reduce_output/5dtx_H.pdb', verbose=True, overwrite=True)
    run_reduce('reduce_input/5dtx.pdb', 'reduce_output/5dtx_H_noflip.pdb', flags='-BUILD -NOFLIP', verbose=True, overwrite=True)
    run_reduce('reduce_input/3ry2.pdb', 'reduce_output/3ry2_H.pdb', verbose=True, overwrite=True)
