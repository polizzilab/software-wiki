#!/bin/bash

#SBATCH -p short # specify which partition. short gets priority and is up to 12 hours, medium up to 5 days
#SBATCH -t 0-5 # time in DD-HH
#SBATCH -c 1 #number of cores-
#SBATCH --mem=2G # memory. start with smaller amount and go up if failed
#SBATCH -o slurm.%x.%A.%a.out # standard out will be printed to this file. %A is the job ID number, %a is the task ID, and %x is the job name
#SBATCH -e slurm.%x.%A.%a.err
#SBATCH --mail-user=jodymou@mit.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-100 # number of lines in input_file_list.txt; note the array should start at 1, not 0

module load gcc/4.8.5
module load rosetta 

FILENAMES="input_file_list.txt" 
# should be a list of names, one per line, of the PDB files to be processed
# often each run uses a different pdb input, xml script, and resfile
# if we name the files as PDB_NAME.pdb, PDB_NAME_flexbb.xml, PDB_NAME_resfile.txt then we can get the files using a list of PDB_NAMES
# make sure #SBATCH --array above is --array=1-NUMBER_OF_LINES_IN_FILE

PDB_NAME=$( awk "NR==$SLURM_ARRAY_TASK_ID" $FILENAMES) # indexes the line of input_file_list.txt e.g. 1st line for 1st job, 2nd line for 2nd job, etc.
PDB_INPUT="${PDB_NAME}.pdb"
LIG_PARAMS="LIG.params"
XML_SCRIPT="${PDB_NAME}_flexbb.xml"
RESFILE="${PDB_NAME}_resfile.txt"

OUTDIR="./"

NSTRUCT=5 # usually increase nstruct to 5 or greater; they will run in series. this means the time needs to be increased

/n/app/rosetta/3.13/source/bin/rosetta_scripts.default.linuxgccrelease \
    -database /n/app/rosetta/3.13/database/ \
    -s $PDB_INPUT \
    -nstruct $NSTRUCT \
    -extra_res_fa $LIG_PARAMS \
    -parser:protocol $XML_SCRIPT \
    -packing:resfile $RESFILE \
    -packing:multi_cool_annealer 10 \
    -packing:linmem_ig 10 \
    -out:path:all $OUTDIR \
    -out:pdb \
    -overwrite \
    -ignore_waters false \
    -beta \
    -water_hybrid_sf \

# if using water vdMs:
# rename all waters to TP5
# use the last 3 options above

