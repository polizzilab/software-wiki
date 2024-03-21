import prody as pr 

constraint_file_path = 'example_inputs/0_pose8_en_0p45_no_CG_top1_of_5.cst'
xml_output_path = 'example_inputs/0_pose8_en_0p45_no_CG_top1_of_5.xml'
xml_template = 'flexbb_water.xml'

p = pr.parsePDB('example_inputs/0_pose8_en_0p45_no_CG_top1_of_5.pdb')

resnums = p.ca.getResnums()
BB_START = min(resnums)
BB_END = max(resnums)

replace_dict = {
    'BB_START': str(BB_START) , 
    'BB_END': str(BB_END), 
    'LIG_START': str(BB_END+1),
    'CONSTRAINT_FILE': constraint_file_path}

with open(xml_template, 'r') as infile:
    with open(xml_output_path, 'w') as outfile:
        for line in infile:
            for o, n in replace_dict.items():
                line = line.replace(o, n)
            outfile.write(line)

