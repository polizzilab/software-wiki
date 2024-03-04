PDB="ABCD.pdb"
LIG_PARAMS="LIG.params"
XML_SCRIPT="flexbb.xml"
RESFILE="resfile.txt"

OUTDIR="./"

NSTRUCT=1

/nfs/sbgrid/programs/x86_64-linux/rosetta/3.13/main/source/bin/rosetta_scripts.default.linuxgccrelease \
    -database /nfs/sbgrid/programs/x86_64-linux/rosetta/3.13/main/database/ \
    -s $PDB \
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