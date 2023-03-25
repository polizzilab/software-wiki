mpnn_path="/nfs/polizzi/shared/programs/design/ProteinMPNN/"

folder_with_pdbs="parsed_struct/"
path_for_parsed_chains=$folder_with_pdbs"/parsed_pdbs.jsonl"
path_for_assigned_chains=$folder_with_pdbs"/assigned_pdbs.jsonl"
path_for_fixed_positions=$folder_with_pdbs"/fixed_pdbs.jsonl"
chains_to_design="Y"
#design_only_positions="1 2 3 4 5 6 7"
#The first amino acid in the chain corresponds to 1 and not PDB residues index for now.

python ${mpnn_path}/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains


python ${mpnn_path}/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --num_seq_per_target 256 \
        --sampling_temp 0.1 \
        --seed 37 \
        --batch_size 128 \
        --save_score 1 \
        --save_probs 1 \
        --model_name v_48_010 \
        --backbone_noise 0.00 \

### brief description of relevant options 
### see protein_mpnn_run.py for full argparse options and descriptions

# out_folder - where to save outputs. will make subfolders probs/ seqs/ scores/ subfolders
# num_seq_per_target - how many sequences to output. note that this might be higher than the number of UNIQUE sequences 
# sampling_temp - 0.1 default. lower T, less diversity; higher T, more diversity
# batch_size - depends on CUDA memory available, bigger proteins need smaller batch sizes. affects speed 
# save_score, save_probs - set 1 for true 0 for false 
# model_name - differs in backbone noise added during training 
        #V_48_002 -> 0.02 Å
        #V_48_010 -> 0.1 Å	
        #V_48_020 -> 0.2 Å <-typically use this one
        #V_48_030 -> 0.3 Å
# backbone noise - gaussian noise to add to the inputs prior to design 
