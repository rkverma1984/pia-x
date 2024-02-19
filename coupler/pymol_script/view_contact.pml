python

import glob
from pymol import cmd
import re
global name_tag

parameters = glob.glob('cont_freq_resids_list/diff_*Ga_*.py')

cmd.do('run contact_distance_general.pml')

color_template = [('orange','greencyan'), ('yellow', 'cyan'), ('wheat','green')]
color_pairs = []
for i in range(len(parameters)):
    color_pairs.append(color_template[i%3])
#####################################################
#---------- Loading contact residues list ----------#
#####################################################

for i in range(len(parameters)):
    positive_color = color_pairs[i][0]
    negative_color = color_pairs[i][1]
    p = parameters[i]
    cmd.do('run '+p)
    #cmd.do('run contact_distance_general.pml')
    
    cmd.do('run draw_arrow_line.py')
    #####################################################
    #---------- Assign values to the function ----------#
    #####################################################
    
    name = re.split(r'[/|\\]',p)[-1].split('.')[0][13:-7]
    name_tag = name[:-2]
    print (name_tag)

    mol_receptor = '/restrain_g//R'
    mol_other = '/restrain_g//A'
    
    cmd.extend('show_residues', show_residues)
    show_residues(mol_receptor, mol_other, receptor_list, galpha_list)
    
    cmd.extend('contact_arrow', contact_arrow)
    contact_arrow(mol_receptor,mol_other,receptor_list, galpha_list, radius_list, name_tag, positive_color, negative_color)


    cmd.group(name, name_tag+'R*')              # grab all distances to a group

python end


hide everything, elem H

set_view (\
     0.006823687,   -0.192558035,   -0.981261492,\
    -0.428242505,    0.886184752,   -0.176879108,\
     0.903637528,    0.421425343,   -0.076414376,\
    -0.000147857,    0.000267588, -150.596496582,\
   120.703155518,  116.674957275,  115.742874146,\
  -569.829589844,  871.022521973,   20.000000000 )
