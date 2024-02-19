from pymol.cgo import *
from chempy import cpv

def show_residues(mol_receptor, mol_other, receptor_list, galpha_list):
    cmd.show('lines', mol_receptor+' and resid '+'+'.join(list(set(receptor_list))))
    cmd.show('lines', mol_other+' and resid '+'+'.join(list(set(galpha_list))))
    return

def contact_arrow(mol_receptor, mol_other, receptor_list, galpha_list, radius_list, name_tag, positive_color, negative_color):
    for i in range(len(receptor_list)):
        selection1 = mol_receptor +' & n. ca & i. '+receptor_list[i]
        selection2 = mol_other +' & n. ca & i. '+galpha_list[i]
        try:
            sel1 = cmd.get_model(selection1,state=1)
            coord1 = sel1.atom[0].coord
            sel2 = cmd.get_model(selection2,state=1)
            coord2 = sel2.atom[0].coord

            radius = radius_list[i]*0.68
            if radius > 0:
                #gap = 0.2
                gap = 0
                #hlength = radius*3.0
                hlength = 0
                #hradius = hlength*0.8
                hradius = 0
                normal = cpv.normalize(cpv.sub(coord1, coord2))
                
                diff = cpv.scale(normal, gap)
                coord1 = cpv.sub(coord1, diff)
                coord2 = cpv.add(coord2, diff)
                coord3 = cpv.add(cpv.scale(normal, hlength), coord2)

                color = list(cmd.get_color_tuple(positive_color))
                color = list(cmd.get_color_tuple(positive_color))

                # arrow direction is from coord1 to coord2
                cylinder_bing = [cgo.CYLINDER] + coord1 + coord3 + [radius] + color + color +\
                [cgo.CONE] + coord3 + coord2 + [hradius, 0.0] + color + color +\
                [1.0, 0.0]
                cmd.load_cgo(cylinder_bing, name_tag+'R'+str(receptor_list[i])+'-A'+str(galpha_list[i]))
            else:
                oppo_coord1 = coord2
                oppo_coord2 = coord1
                oppo_radius = 0 - radius
                #gap = oppo_radius
                gap = 0
                #hlength = oppo_radius*3.0
                hlength = 0
                #hradius = hlength*0.8
                hradius = 0
                normal = cpv.normalize(cpv.sub(oppo_coord1, oppo_coord2))
                
                diff = cpv.scale(normal, gap)
                oppo_coord1 = cpv.sub(oppo_coord1, diff)
                oppo_coord2 = cpv.add(oppo_coord2, diff)
                oppo_coord3 = cpv.add(cpv.scale(normal, hlength), oppo_coord2)

                color = list(cmd.get_color_tuple(negative_color))
                color = list(cmd.get_color_tuple(negative_color))

                cylinder_bing = [cgo.CYLINDER] + oppo_coord1 + oppo_coord3 + [oppo_radius] + color + color +\
                [cgo.CONE] + oppo_coord3 + oppo_coord2 + [hradius, 0.0] + color + color +\
                [1.0, 0.0]
                cmd.load_cgo(cylinder_bing, name_tag+'R'+str(receptor_list[i])+'-A'+str(galpha_list[i]))
                
        except:
            print (selection1,selection2, 'does not work')
            
            
    cmd.center('all')
    return

