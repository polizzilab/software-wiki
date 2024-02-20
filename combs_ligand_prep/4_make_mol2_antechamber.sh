# Use antechamber to correct partial charges before make a params file.

INPUT="LIG.mol2"
OUTPUT="LIG_antechamber.mol2"

antechamber -i $INPUT -fi mol2 -o $OUTPUT -fo mol2 -c gas -at sybyl -an n -pf y 