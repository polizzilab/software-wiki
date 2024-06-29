# align two ligands using LS align

'''
usage: python align_ligands.py lig1.pdb lig2.pdb lig2_reordered.pdb --output_dir output_dir

This script aligns two ligands using LS align. The output is a pdb file with the second ligand reordered to match the first ligand. The second ligand is also renamed.
The output is saved in the output_dir.
'''

import argparse
import os
import sys
import prody as pr

LSalign_path ='/nfs/polizzi/shared/programs/chemistry/LSalign/LSalign'

def align_ligands(lig1, lig2,lig2_reordered, output_dir):

    # convert the pdb files to mol2
    lig1_mol2 = lig1.replace(".pdb", ".mol2")
    lig2_mol2 = lig2.replace(".pdb", ".mol2")

    os.system(f'obabel {lig1} -O {lig1_mol2}')
    os.system(f'obabel {lig2} -O {lig2_mol2}')
    lig1_name = os.path.basename(lig1).split('.')[0]
    lig2_name = os.path.basename(lig2).split('.')[0]
    output_txt = os.path.join(output_dir, f'align_{lig1_name}_{lig2_name}.txt')
    output_pdb = os.path.join(output_dir, f'align_{lig1_name}_{lig2_name}.pdb')
    print(f'{LSalign_path} {lig1_mol2} {lig2_mol2} -rf 1 -md 1 -H 1 -a {output_txt} -o {output_pdb} -d0 2')
    os.system(f'{LSalign_path} {lig1_mol2} {lig2_mol2} -rf 1 -md 1 -H 1 -a {output_txt} -o {output_pdb} -d0 2')

    align = open(output_txt, 'r')
    atom_indices = []

    for line in align:
        if 'Templ Atom Index:' in line:
            atom_index = line.split()[3:]
            # convert the string to integer
            atom_index = [int(i) for i in atom_index]
            atom_indices.extend(atom_index)
    
    # these indices are not 0 indexed so subtract 1
    atom_indices = [i-1 for i in atom_indices]

    import itertools
    import operator

    lig1 = pr.parsePDB(lig1)
    lig2 = pr.parsePDB(lig2)

    reordered = []
    for i in atom_indices:
        atom = lig2.select(f'index {i}').toAtomGroup()
        reordered.append(atom)
    reordered = list(itertools.accumulate(reordered, operator.add))[-1]
    # rename all the atoms too
    lig2.setNames(lig1.getNames()[0:len(reordered)])

    # set title
    reordered.setTitle(lig2_name)

    pr.writePDB(lig2_reordered, reordered)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("lig1", help="first pdb file")
    parser.add_argument("lig2", help="first pdb file")
    parser.add_argument("lig2_reordered")
    parser.add_argument("--output_dir", default='./', help="output directory")
    args = parser.parse_args()
    align_ligands(args.lig1, args.lig2, args.lig2_reordered, args.output_dir)

if __name__ == "__main__":
    main()