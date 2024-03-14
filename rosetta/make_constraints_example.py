import prody as pr
from constraints import make_constraints, pymol_visualize_constraints_print

pdb_path = 'example_inputs/0_pose8_en_0p45_no_CG_top1_of_5.pdb'

p = pr.parsePDB(pdb_path)

sele_str = 'not resname GLY' # ideally select ligand and vdMs by resnum -- for now, just select all non-glycine backbones. this will miss glycine-backbone interactions

cst_list = make_constraints(
    p, 
    sele_str, 
    output_constraint_file='example_inputs/0_pose8_en_0p45_no_CG_top1_of_5.cst', 
    return_constraint_list=True
    )

# function that writes a pymol session file to visualize constraints on input PDB
# uses pymol cmd 

# pymol_visualize_constraints(
#     pdb_path=pdb_path,
#     cst_list = cst_list,
#     output_pml_file='example_inputs/0_pose8_en_0p45_no_CG_top1_of_5.pml'
#     )

# simpler vizualization that prints pymol commands 
# load input pdb and copy-paste the printed commands into pymol console
pymol_visualize_constraints_print(cst_list=cst_list)