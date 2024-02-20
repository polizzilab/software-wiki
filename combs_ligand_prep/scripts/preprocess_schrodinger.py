import argparse 
from schrodinger import structure
from schrodinger.structutils import build

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="")
parser.add_argument("--output", type=str)
parser.add_argument("--add_h", action="store_true")
args = parser.parse_args() 

st = structure.StructureReader.read(args.input)

if args.add_h:
    build.add_hydrogens(st) 

with structure.StructureWriter(args.output) as writer:
    writer.append(st)

# use $SCHRODINGER/run preprocess_schrodinger.py 