show cartoon
color gray80, elem C

# show sticks for sidechains that were designed
show sticks, segi X+Y or name CA

# hide the chemical group fragments of vdMs
hide everything, segi Y

# color ligand purple and sidechains green
color purple, elem C and segi L 
color green, elem C and segi X 

# show the chemical group fragments as cyan spheres, set smaller size
show spheres, segi Y 
set sphere_scale, 0.25
color cyan, elem C and segi Y 

# reveal and color any designed glycines 
show spheres, segi X and resn GLY 
color orange, segi X and resn GLY 

# label with the vdM scores in the B-factors
label name CB and segi X, "%s \n %.2f" % (resn, b)
label name HA2 and segi X, "%s \n %.2f" % (resn, b)
set label_size, -0.5
set label_position, (1.5,1.5,1.5)

