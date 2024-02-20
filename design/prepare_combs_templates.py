'''
Usage: prepare_templates.py <template directory> 

This script prepares all-alanine and all-glycine templates for use with combs in parallel.

Inputs:
    args.template_dir: a directory containing .pdb files. All pdb files in this folder will be processed
    args.submit: if True, submit the jobs to the cluster. If False, only print the commands it would have run. Can also use this to remove the termini only.
    args.verbose: if True, the script will print the commands it would have run
    args.threads: number of threads to use
    args.chunk_size: number of templates to process at once

The outputs are written to all_gly/ and all_ala/ subdirectories of the input directory.
They are named the same as the input pdb file.

'''
import prody as pr
import os
import subprocess
import glob
import multiprocessing

FIXBB_EXECUTABLE = '/nfs/sbgrid/programs/x86_64-linux/rosetta/3.13/main/source/bin/fixbb.default.linuxgccrelease'
ROSETTA_DATABASE = '/nfs/sbgrid/programs/x86_64-linux/rosetta/3.13/main/database'

def prepare_templates(chunk_ID, template_path_ls):
    #infilename = ' '.join(template_path_ls)
    #outpath = Path(infilename).parent
    for infilename in template_path_ls:
        for xxx in ['gly', 'ala']:
            resfile = f'./resfile_{xxx}.txt'
            command = [
                    FIXBB_EXECUTABLE,
                    '-s', infilename,
                    '-database', ROSETTA_DATABASE,
                    '-resfile', resfile,
                    '-out:no_nstruct_label',
                    '-out:pdb',
                    '-out:path:all', f'{template_dir}/all_{xxx}/',
                    '-out:overwrite',
                    '-nstruct', '1',
            ]

            print('\t', ' '.join(command))
            if args.verbose:
                print(f'Worker {chunk_ID}: creating all_{xxx} template')
                if args.submit:
                    subprocess.run(command, check=True)
            else: 
                if args.submit:
                    subprocess.run(command, check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

            # remove the termini if you'd like
            if args.remove_termini:
                infile_base = os.path.basename(infilename)
                p = pr.parsePDB(f'{template_dir}/all_{xxx}/{infile_base}')
                p = p.select('not name 1H 2H 3H OXT')
                pr.writePDB(f'{template_dir}/all_{xxx}/{infile_base}', p)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--template_dir',type=str)
    parser.add_argument('--submit', action='store_true', default=True)
    parser.add_argument('--verbose', action='store_true', default=False)
    parser.add_argument('--threads', type=int, default=5)
    parser.add_argument('--chunk_size', type=int, default=5)
    parser.add_argument('--remove_termini', action='store_true', default=False) # removes 1H 2H 3H so combs does not use them
    args = parser.parse_args() 

    template_dir = args.template_dir

    os.makedirs(os.path.join(template_dir, 'all_ala'), exist_ok=True)
    os.makedirs(os.path.join(template_dir, 'all_gly'), exist_ok=True)

    template_path_ls = glob.glob(os.path.join(template_dir, '*.pdb'))
    
    print(template_path_ls)

    chunks = [template_path_ls[i:i+args.chunk_size] 
                for i in range(0, len(template_path_ls), args.chunk_size)]

    print(f'Splitting into {len(chunks)} chunks of size', args.chunk_size)

    if args.submit:
        with multiprocessing.Pool(args.threads) as p:
            p.starmap(prepare_templates, enumerate(chunks))
    



