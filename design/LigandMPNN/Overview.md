**How to run LigandMPNN in the Polizzi lab:**

**Example 1:** Redesign ClpS protein without any fixed residues

-With this example, I use LigandMPNN with a temperature of 0.2 to redesign the entire input structure and provide one output structure
-To run this python command, I first activate Jody's pytorch conda environment on the workstation. **conda activate pytorch**

`python /nfs/polizzi/shared/programs/design/LigandMPNN/run.py --model_type "ligand_mpnn" --checkpoint_ligand_mpnn /nfs/polizzi/shared/programs/design/LigandMPNN/model_params/ligandmpnn_v_32_020_25.pt --pdb_path /nfs/polizzi/adharani/npl-002-yeast-display-library/240226_LigandMPNN_specificity/3DNJ_test_ligandMPNN_specificity/3DNJ_FAA.pdb --out_folder /nfs/polizzi/adharani/npl-002-yeast-display-library/240226_LigandMPNN_specificity/3DNJ_test_ligandMPNN_specificity/test_LigandMPNN --seed 111 --batch_size 1 --number_of_batches 1 --temperature 0.2 --save_stats 1`

**Output from this command:**
`Skipping these pdb files since the outputs already exist: []
Running these pdb files: ['/nfs/polizzi/adharani/npl-002-yeast-display-library/240226_LigandMPNN_specificity/3DNJ_test_ligandMPNN_specificity/3DNJ_FAA.pdb']
Designing protein from this path: /nfs/polizzi/adharani/npl-002-yeast-display-library/240226_LigandMPNN_specificity/3DNJ_test_ligandMPNN_specificity/3DNJ_FAA.pdb
These residues will be redesigned:  ['A39', 'A40', 'A41', 'A42', 'A43', 'A44', 'A45', 'A46', 'A47', 'A48', 'A49', 'A50', 'A51', 'A52', 'A53', 'A54', 'A55', 'A56', 'A57', 'A58', 'A59', 'A60', 'A61', 'A62', 'A63', 'A64', 'A65', 'A66', 'A67', 'A68', 'A69', 'A70', 'A71', 'A72', 'A73', 'A74', 'A75', 'A76', 'A77', 'A78', 'A79', 'A80', 'A81', 'A82', 'A83', 'A84', 'A85', 'A86', 'A87', 'A88', 'A89', 'A90', 'A91', 'A92', 'A93', 'A94', 'A95', 'A96', 'A97', 'A98', 'A99', 'A100', 'A101', 'A102', 'A103', 'A104', 'A105', 'A106', 'A107', 'A108', 'A109', 'A110', 'A111', 'A112', 'A113', 'A114', 'A115', 'A116', 'A117', 'A118', 'A119']
These residues will be fixed:  []
The number of ligand atoms parsed is equal to: 21
Type: N, Coords [17.225 -6.148 12.885], Mask 1
Type: C, Coords [17.813 -6.646 11.629], Mask 1
Type: C, Coords [18.899 -5.722 11.135], Mask 1
Type: O, Coords [19.14  -4.662 11.702], Mask 1
Type: C, Coords [16.723 -6.801 10.585], Mask 1
Type: C, Coords [15.61  -7.703 11.029], Mask 1
Type: C, Coords [15.853 -9.045 11.277], Mask 1
Type: C, Coords [14.326 -7.204 11.228], Mask 1
Type: C, Coords [14.831 -9.89  11.695], Mask 1
Type: C, Coords [13.304 -8.035 11.634], Mask 1
Type: C, Coords [13.557 -9.391 11.868], Mask 1
Type: N, Coords [19.551 -6.148 10.072], Mask 1
Type: C, Coords [20.699 -5.429  9.527], Mask 1
Type: C, Coords [20.348 -4.762  8.219], Mask 1
Type: O, Coords [19.806 -5.403  7.313], Mask 1
Type: C, Coords [21.857 -6.432  9.393], Mask 1
Type: N, Coords [20.612 -3.446  8.119], Mask 1
Type: C, Coords [20.477 -2.767  6.834], Mask 1
Type: C, Coords [21.683 -3.018  5.961], Mask 1
Type: O, Coords [22.788 -3.148  6.495], Mask 1
Type: C, Coords [20.247 -1.273  7.118], Mask 1`

****Example 2: **** Redesign ClpS protein with three fixed residues

-With this example, I use LigandMPNN with a temperature of 0.2 to redesign the entire input structure with three residues fixed (Chain A and residues 47 A49 A79) and provide five output structures
-To run this python command, I first activate Jody's pytorch conda environment on the workstation. **conda activate pytorch**

`python /nfs/polizzi/shared/programs/design/LigandMPNN/run.py  --model_type "ligand_mpnn" --checkpoint_ligand_mpnn /nfs/polizzi/shared/programs/design/LigandMPNN/model_params/ligandmpnn_v_32_020_25.pt --pdb_path /nfs/polizzi/adharani/npl-002-yeast-display-library/240226_LigandMPNN_specificity/3DNJ_test_ligandMPNN_specificity/3DNJ_FAA.pdb --out_folder /nfs/polizzi/adharani/npl-002-yeast-display-library/240226_LigandMPNN_specificity/3DNJ_test_ligandMPNN_specificity/test_LigandMPNN_5_seqs --seed 111 --batch_size 5 --number_of_batches 1 --temperature 0.2 --save_stats 1 --fixed_residues "A47 A49 A79"`

**Output from this command:**
`Skipping these pdb files since the outputs already exist: []
Running these pdb files: ['/nfs/polizzi/adharani/npl-002-yeast-display-library/240226_LigandMPNN_specificity/3DNJ_test_ligandMPNN_specificity/3DNJ_FAA.pdb']
Designing protein from this path: /nfs/polizzi/adharani/npl-002-yeast-display-library/240226_LigandMPNN_specificity/3DNJ_test_ligandMPNN_specificity/3DNJ_FAA.pdb
These residues will be redesigned:  ['A39', 'A40', 'A41', 'A42', 'A43', 'A44', 'A45', 'A46', 'A48', 'A50', 'A51', 'A52', 'A53', 'A54', 'A55', 'A56', 'A57', 'A58', 'A59', 'A60', 'A61', 'A62', 'A63', 'A64', 'A65', 'A66', 'A67', 'A68', 'A69', 'A70', 'A71', 'A72', 'A73', 'A74', 'A75', 'A76', 'A77', 'A78', 'A80', 'A81', 'A82', 'A83', 'A84', 'A85', 'A86', 'A87', 'A88', 'A89', 'A90', 'A91', 'A92', 'A93', 'A94', 'A95', 'A96', 'A97', 'A98', 'A99', 'A100', 'A101', 'A102', 'A103', 'A104', 'A105', 'A106', 'A107', 'A108', 'A109', 'A110', 'A111', 'A112', 'A113', 'A114', 'A115', 'A116', 'A117', 'A118', 'A119']
These residues will be fixed:  ['A47', 'A49', 'A79']
The number of ligand atoms parsed is equal to: 21
Type: N, Coords [17.225 -6.148 12.885], Mask 1
Type: C, Coords [17.813 -6.646 11.629], Mask 1
Type: C, Coords [18.899 -5.722 11.135], Mask 1
Type: O, Coords [19.14  -4.662 11.702], Mask 1
Type: C, Coords [16.723 -6.801 10.585], Mask 1
Type: C, Coords [15.61  -7.703 11.029], Mask 1
Type: C, Coords [15.853 -9.045 11.277], Mask 1
Type: C, Coords [14.326 -7.204 11.228], Mask 1
Type: C, Coords [14.831 -9.89  11.695], Mask 1
Type: C, Coords [13.304 -8.035 11.634], Mask 1
Type: C, Coords [13.557 -9.391 11.868], Mask 1
Type: N, Coords [19.551 -6.148 10.072], Mask 1
Type: C, Coords [20.699 -5.429  9.527], Mask 1
Type: C, Coords [20.348 -4.762  8.219], Mask 1
Type: O, Coords [19.806 -5.403  7.313], Mask 1
Type: C, Coords [21.857 -6.432  9.393], Mask 1
Type: N, Coords [20.612 -3.446  8.119], Mask 1
Type: C, Coords [20.477 -2.767  6.834], Mask 1
Type: C, Coords [21.683 -3.018  5.961], Mask 1
Type: O, Coords [22.788 -3.148  6.495], Mask 1
Type: C, Coords [20.247 -1.273  7.118], Mask 1`
