# After centering, make a mol2 file using schrodinger.
INPUT="LIG_centered.pdb" 
OUTPUT="LIG.mol2"

# schrodinger structconvert can convert between many formats. can also convert to sdf 
$SCHRODINGER/utilities/structconvert $INPUT $OUTPUT 