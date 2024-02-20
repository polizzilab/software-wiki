import prody as pr 
import numpy as np 
import argparse

# script to make sure that ligand pdbs have the same segment ID, resnum, and chain ID

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="")
parser.add_argument("--output", type=str)
parser.add_argument("--center", action="store_true")
args = parser.parse_args() 

lig = pr.parsePDB(args.input)
lig.setSegnames('L')
lig.setChids('L')
lig.setResnums(10)

if args.center:
    center = pr.calcCenter(lig)
    # move the center of the mass of the ligand to the origin
    pr.moveAtoms(lig, to=np.zeros(3))
pr.writePDB(args.output, lig)
