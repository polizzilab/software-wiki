# Use schrodinger to add hyodrogens and correct aromaticity issues.

INPUT="LIG.pdb" # can be mol2, pdb, any format
OUTPUT="LIG_schrodinger.pdb"

$SCHRODINGER/run scripts/preprocess_schrodinger.py  
    --input $INPUT
    --add_h 
    --output $OUTPUT