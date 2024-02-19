#####################################################
#------- Loading template structure ----------------#
#####################################################

load 6ddf_xtal.pdb, template_6ddf_main
color yellow, template_6ddf_main
show cartoon,template_6ddf_main
disable template_6ddf_main

load Gwhole_restrain.pdb, restrain_g
align restrain_g and chain R and name CA, template_6ddf_main and chain R and name CA

color white,restrain_g
color grey60, chain R


#####################################################
#------- Loading Draw arrow lines Function ---------#
#####################################################

#run draw_arrow_line.py

######################################
#----------- cylinder MOR -----------#
######################################

set cartoon_cylindrical_helices,1
select tm1, chain R and resid 73-81
select tm2, chain R and resid 116-125
select tm3, chain R and resid 147-153
select tm4, chain R and resid 192-197
select tm5, chain R and resid 239-246
select tm6, chain R and resid 292-299
select tm7, chain R and resid 322-329
alter i. 77+121+150+195+239+292+325, ss='L'
######################################

deselect
color atomic, not elem C
hide everything, elem H

delete tm*
disable template_6ddf_main


