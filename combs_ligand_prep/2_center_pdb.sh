# Use prody to center the ligand at origin and rename entire ligand to chain L segment L resnum 10.

INPUT="LIG_schrodinger.pdb"
OUTPUT="LIG_centered.pdb"

/nfs/polizzi/jmou/miniconda3/envs/pytorch/bin/python scripts/prody_fix.py \
    --input $INPUT \
    --output $OUTPUT \
    --center
