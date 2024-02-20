# Use Rosetta's existing function to make a params file.

INPUT="LIG_antechamber.mol2"
OUTPUT_FOLDER="./"
LIG_NAME="LIG"

python /nfs/sbgrid/programs/x86_64-linux/rosetta/3.13/main/source/scripts/python/public/molfile_to_params.py \
    $INPUT \
    --keep-names -n $LIG_NAME \
    --pdb $OUTPUT_FOLDER/$LIG_NAME \
    --no-pdb
