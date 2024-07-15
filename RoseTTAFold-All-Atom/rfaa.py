import os
import yaml
import subprocess

# remember to activate environment - run.sh takes care of that for you if you wish
# conda activate /nfs/polizzi/shared/miniconda3/envs/SE3nv

# Adjust based on your ligand
LIG = "SRO"
yaml_file = f"protein_sm_singleseq_{LIG}.yaml" 
sdf_file = f"{LIG}_ideal.sdf" 

# current directory
pwd = os.getcwd()
print(pwd)

# Path to change to
path = "/nfs/polizzi/shared/programs/structpred/RoseTTAFold-All-Atom_JV" # I made a copy for myself because I started having issues with my writing permit.
# Change the current working directory
os.chdir(path)

# Constants
FASTA_DIR = f"{pwd}/pdbs_fasta" # this is the directory where your fasta files are

YAML_PATH = f"{path}/rf2aa/config/inference/{yaml_file}"
print(YAML_PATH)

# Read the original YAML template
with open(YAML_PATH, 'r') as file:
    original_yaml = file.read()

# Iterate over each fasta file in the directory
for fasta_file in os.listdir(FASTA_DIR):
    if fasta_file.endswith('.fasta'):
        base_name = fasta_file[:-6]  # Removes the '.fasta' extension

        # Load the yaml content
        yaml_content = yaml.safe_load(original_yaml)

        # Update the YAML content
        yaml_content['job_name'] = f"{base_name}_rfaa"
        yaml_content['loader_params']['MAXCYCLE'] = 10 # adjust based on your ligand, so far 10 has worked fine for me
        yaml_content['output_path'] = f"{pwd}/rfaa_output"
        yaml_content['protein_inputs']['A']['fasta_file'] = f"{FASTA_DIR}/{fasta_file}"
        yaml_content['sm_inputs']['B']['input'] = f"{pwd}/{sdf_file}"

        # Write the updated YAML content back to the file
        with open(YAML_PATH, 'w') as file:
            yaml.dump(yaml_content, file)

        # Execute the script
        config_name = yaml_file.replace('.yaml', '')
        command = f"python -m rf2aa.run_inference --config-name {config_name}"
        subprocess.run(command, shell=True)

        # Restore the original YAML content after each run
        with open(YAML_PATH, 'w') as file:
            file.write(original_yaml)
