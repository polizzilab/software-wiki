import prody as pr
from calculate_resfile import write_resfile
from calculate_topology import calculate_and_write_topology

# first calculate the topology
# this only needs to be done once for each unique backbone
p_helices = pr.parsePDB('example_inputs/helicesAll27_BvB_template.pdb')
topo_charge_dict_path = 'example_inputs/0_pose8_en_0p45_no_CG_top1_of_5_topo_charge_dict.txt'
#calculate_and_write_topology(p_helices, calpha_distance=12, alpha=20, num_iterations=1000000, topo_charge_dict_path=topo_charge_dict_path)

# then write the resfile
p = pr.parsePDB('example_inputs/0_pose8_en_0p45_no_CG_top1_of_5.pdb')

# output path for resfile
resfile_output_path = 'example_inputs/0_pose8_en_0p45_no_CG_top1_of_5_resfile.txt' 

# text file with one residue per line, "chain resnum" e.g. "A 10" describing vdM residues
vdm_residues_path = 'example_inputs/0_pose8_en_0p45_no_CG_top1_of_5_vdms_residues.txt'

write_resfile(p, 
              alpha=9, 
              distance_threshold=-1.0, 
              topo_charge_dict_path=topo_charge_dict_path, 
              designable_residues_path=None, # no input here, so all residues other than vdMs will be designable
              vdm_residues_path=vdm_residues_path, 
              resfile_output_path=resfile_output_path, 
              vdm_rosetta_rotamer='NATAA')