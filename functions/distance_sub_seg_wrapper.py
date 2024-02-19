"""
Description:
    wrapper script to calculate the distance between the COM of CA from helical segments.
    
    Dependencies:
        distance_angle_calc.py
        center_of_mass_calc.py
        distance_angle_calc.py


Author
    Ravi kumar verma
"""

import numpy as np
from collections import OrderedDict
import center_of_mass_calc
import distance_angle_calc
import progress_bar


def sub_seg_distance(ltraj, p_name, st, traj, parms, matrix, coord_dict, dict_frames, dict_average,
                     individual_matrix, total_frames, separate_output_option, mass):
    """
    DESCRIPTION

        Calculate to pairwise distances between sub segments.

    ARGUMENTS
        ltraj = loaded trajectory
        p_name = protein name
        traj   = name of the trajectory
        dict_frames = dictionary with list of frame_ids (specified by the cluster to be used)
        parms.helix_index[p_name] = Dictionary containing start and end positions of the helical segments.
        parms.helix_names[p_name]  = Dictionary containing name of the helical segments.
        meta_dict   = Meta Dictionary used to store all the analysis.

    Author:
        Ravi Kumar Verma

    """
    
    for i in range(0, len(parms.helix_index[p_name])):
        tm = parms.helix_index[p_name][i]
        tm_n = parms.helix_names[p_name][i]
        
        if parms.atomselection_com_TM == 'CA':
            indices = ltraj.topology.select('protein and (resSeq ' + str(tm[0]) + ' to ' + str(tm[1]) + ' and name CA)')
        if parms.atomselection_com_TM == 'all':
            indices = ltraj.topology.select('protein and (resSeq ' + str(tm[0]) + ' to ' + str(tm[1]) + ')')
        if parms.atomselection_com_TM == 'all_heavy':
            indices = ltraj.topology.select('protein and (resSeq ' + str(tm[0]) + ' to ' + str(tm[1]) + ' and not element H)')
        
        # print indices
        # Calculate the center of the mass, change mass==None to mass=='atomic' to take care of atomic masses.
        coord_dict[p_name][st][traj][matrix][tm_n].update(center_of_mass_calc.calculate_com(ltraj, dict_frames[p_name][st][traj], indices, mass))
        # print (p_name)
        # print (tm, tm_n, center_of_mass_calc.calculate_com(ltraj, dict_frames[p_name][st][traj], indices, mass))
        # print (dict_frames[p_name][st][traj])
    
    for fr in dict_frames[p_name][st][traj]:
        temp = OrderedDict()
        arr_done = []
        sorted_keys = []
        
        for i in range(0, len(parms.helix_names[p_name])):
            
            t1 = parms.helix_names[p_name][i]
            for j in range(0, len(parms.helix_names[p_name])):
                t2 = parms.helix_names[p_name][j]
                sorted_keys.append(t1 + '_' + t2)
                
                dict_distance = OrderedDict()
                if t1 != t2:
                    if t1 + '_' + t2 not in arr_done and t2 + '_' + t1 not in arr_done:
                        # Calculate the pairwise distances between sub segments in a single frame and store in the dictionary.
                        dict_distance.update(distance_angle_calc.get_distance(coord_dict[p_name][st][traj][matrix][t1].get(fr),
                                                                              coord_dict[p_name][st][traj][matrix][t2].get(fr), fr))
                        arr_done.append(t1 + '_' + t2)
                        temp.update({t1 + '_' + t2: dict_distance})
                        temp.update({t2 + '_' + t1: dict_distance})
                elif t1 == t2:
                    dict_distance.update({fr: 0.0})
                    arr_done.append(t1 + '_' + t2)
                    temp.update({t1 + '_' + t2: dict_distance})
                    temp.update({t2 + '_' + t1: dict_distance})
        #                meta_dict[p_name][st][traj][matrix][t1 + '_' + t2].update(dict_distance)
        #                meta_dict[p_name][st][traj][matrix][t2 + '_' + t1].update(dict_distance)
        
        # Calculate the pairwise distances between sub segments in a single frame and store in the dictionary.
        
        concatenated_matrix = []
        for key in sorted_keys:
            # append all the pairwise distances for a single frame in a single array
            data = temp[key].get(fr)
            #            data = meta_dict[p_name][st][traj][matrix][key].get(fr)
            concatenated_matrix.append(data)
        
        # create a len(parms.helix_names[p_name]) * len(parms.helix_names[p_name]) matrix
        arr_matrix = [concatenated_matrix[i:i + len(parms.helix_names[p_name])] for i in
                      np.arange(0, len(concatenated_matrix), len(parms.helix_names[p_name]))]
        
        dict_average[st][matrix].append(np.array(arr_matrix))
        
        individual_matrix[p_name][st][matrix][traj][fr] = np.array(arr_matrix)
        
        if separate_output_option == 'yes':
            # write matrix for each of the frame for a given traj.
            np.savetxt(("%s/%s.%s_%s_%s_%s.csv" % (parms.outdir, p_name, st, matrix, traj, str(fr))), arr_matrix, delimiter=',', fmt="%s")
        elif separate_output_option == 'no':
            pass
        
        progress_bar.cli_progress(dict_average, total_frames, parms.metrics, bar_length=20)
    
    return dict_average
