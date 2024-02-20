
PDB_INPUT="protein.pdb"
RESFILE="resfile_ala.txt" # or resfile_gly.txt
OUTDIR="./"

/nfs/sbgrid/programs/x86_64-linux/rosetta/3.13/main/source/bin/rosetta_scripts.default.linuxgccrelease \
    -database /nfs/sbgrid/programs/x86_64-linux/rosetta/3.13/main/database/ \
    -s $PDB_INPUT \
    -resfile $RESFILE \
    -nstruct 1 \
    -out:path:all $OUTDIR \
    -out:no_nstruct_label \
    -out:pdb \
    -overwrite
