###
# Run graph matcher to get a minimal scaffold for the ligand.txt file. 
# This uses the pytorch conda environment.
# If running into conda issues, try export PATH="/nfs/polizzi/jmou/miniconda3/bin:$PATH",
# and then deactivate and reactivate the environment.
# you will have to manually fix the ligand.txt after.
# use this for help! https://github.com/npolizzi/Combs2/wiki/Reference:-ligand.txt-specification
###

/nfs/polizzi/jmou/miniconda3/envs/pytorch/bin/python \
    scripts/graph_matcher.py \
    --pdb_path LIG_schrodinger.pdb \
    --ligand_txt_out LIG_raw.txt 

