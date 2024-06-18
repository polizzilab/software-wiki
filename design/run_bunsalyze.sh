#!/usr/bin/env bash
# Example of how to run bunsalyze with its command line arguments
# See the flags of `bunsalyze --help` for more details

mkdir -p bunsalyze_output

# By default results are written to stdout; if you want to save it, redirect it to a file
./bunsalyze bunsalyze_input/3dnj_prot_lig.pdb bunsalyze_input/3dnj_dons_accs.txt \
    --ligand_chain C \
    > bunsalyze_output/3dnj_buns.json

# Intermediate results can be saved with these flags
./bunsalyze bunsalyze_input/3dnj_prot_lig.pdb bunsalyze_input/3dnj_dons_accs.txt \
    --ligand_chain C \
    --probe_output_dir bunsalyze_output/3dnj \
    --freesasa_output_dir bunsalyze_output/3dnj \
    > bunsalyze_output/3dnj_buns.json

# Here is another example for streptavidin-biotin
./bunsalyze bunsalyze_input/3ry2_A_H.pdb bunsalyze_input/BTN_dons_and_accs.txt \
    > bunsalyze_output/3ry2_buns.json

# You can use flags to use COMBSian clash detection and burial detection, instead of probe/freesasa
# This may take a while to get going, because COMBS takes a long time to import.
./bunsalyze bunsalyze_input/3ry2_A_H.pdb bunsalyze_input/BTN_dons_and_accs.txt \
    --use_alpha_hull --use_combsian_contact --vdW_tolerance 0.15 \
    --ligand_resname BTN --ligand_params bunsalyze_input/BTN.params \
    > bunsalyze_output/3ry2_buns_combsian.json

# The script also works for a directory containing many .pdb files
# The -j flag will parallelize
# The output of this is a *.jsonl* (json lines format)
./bunsalyze bunsalyze_input/laserfrogs_OBI_parabundles/ \
    bunsalyze_input/OBI_dons_accs.txt \
    --run_reduce \
    --reduce_output_dir bunsalyze_output/laserfrogs_pdb_H/ \
    --use_alpha_hull \
    --freesasa_output_dir bunsalyze_output/laserfrogs_burial/ \
    --use_combsian_contact \
    --probe_output_dir bunsalyze_output/laserfrogs_contacts/ \
    --ligand_resname OBI --ligand_params bunsalyze_input/OBI.params \
    -j 4 \
    > bunsalyze_output/laserfrogs.jsonl
