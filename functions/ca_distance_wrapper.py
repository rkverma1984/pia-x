"""
Description:
    wrapper script to calculate the pairwise distances between the CA of binding site residues.
    
    Dependencies:
        distance_angle_calc
        center_of_mass_calc

Author
    Ravi kumar verma
"""


from collections import OrderedDict
import numpy as np
import center_of_mass_calc
import distance_angle_calc
import progress_bar


def binding_res_ca_distance(ltraj, p_name, st, traj, matrix, parms, coord_dict, dict_frames, dict_average, individual_matrix, total_frames,
                            separate_output_option, mass):
    """
    DESCRIPTION

        Calculate to pairwise distances.

    ARGUMENTS
        
        :param ltraj:loaded trajectory
        :param st: state id working on
        :param traj: name of the traj working on.
        :param matrix: matrix calculation to be performed.
        :param parms: parameters file
        :param coord_dict: dictionary to store coordinates or center of mass data.
        :param dict_frames: dictionary with list of frame_ids (specified by the cluster to be used)
        :param dict_average: dictionary to save data for individual state for each matrix
        :param individual_matrix: dictionary to store individual matrices
        :param total_frames: total no of frames on which analysis is to be performed
        :param separate_output_option: write individual matrices for each frame if 'yes'.
        :param mass:

        :type p_name: object
        binding_index[p_name] = Dictionary containing residue numbering for the binding site residues.
        binding_BW[p_name]  = Dictionary containing Ballesteros-Wienstein numbering for the binding site residues.
        
    Author:
        Ravi Kumar Verma
    """

    for i in range(0, len(parms.binding_index[p_name])):
        b_res = parms.binding_index[p_name][i]
        b_n = str(parms.binding_BW[p_name][i])
        indices = ltraj.topology.select('protein and (resSeq ' + str(b_res) + ' and name CA)')
        
        # Calculate the center of the mass, change mass==None to mass=='atomic' to take care of atomic masses.
        coord_dict[p_name][st][traj][matrix][b_n].update(
            center_of_mass_calc.calculate_com(ltraj, dict_frames[p_name][st][traj], indices, mass))
    
    for fr in dict_frames[p_name][st][traj]:
        temp = OrderedDict()
        arr_done = []
        sorted_keys = []
        for i in range(0, len(parms.binding_BW[p_name])):
            b1 = str(parms.binding_BW[p_name][i])
            for j in range(0, len(parms.binding_BW[p_name])):
                b2 = str(parms.binding_BW[p_name][j])
                sorted_keys.append(b1 + '_' + b2)
                
                dict_distance = OrderedDict()
                if b1 != b2:
                    if b1 + '_' + b2 not in arr_done and b2 + '_' + b1 not in arr_done:
                        # calculate the pairwise distances between CA atoms in a single frame and store in the dictionary.
                        dict_distance.update(
                            distance_angle_calc.get_distance(coord_dict[p_name][st][traj][matrix][b1].get(fr),
                                                             coord_dict[p_name][st][traj][matrix][b2].get(fr), fr))
                        arr_done.append(b1 + '_' + b2)
                        temp.update({b1 + '_' + b2: dict_distance})
                        temp.update({b2 + '_' + b1: dict_distance})
                
                elif b1 == b2:
                    dict_distance.update({fr: 0.0})
                    arr_done.append(b1 + '_' + b2)
                    temp.update({b1 + '_' + b2: dict_distance})
                    temp.update({b2 + '_' + b1: dict_distance})
#                meta_dict[p_name][st][traj][matrix][b1 + '_' + b2].update(dict_distance)
#                meta_dict[p_name][st][traj][matrix][b2 + '_' + b1].update(dict_distance)
        
        concatenated_matrix = []
        for key in sorted_keys:
            # append all the pairwise distances for a single frame in a single array
            data = temp[key].get(fr)
#            data = meta_dict[p_name][st][traj][matrix][key].get(fr)
            concatenated_matrix.append(data)
        
        # create a len(parms.binding_BW[p_name]) * len(parms.binding_BW[p_name]) matrix
        arr_matrix = [concatenated_matrix[i:i + len(parms.binding_BW[p_name])] for i in np.arange(0, len(concatenated_matrix), len(parms.binding_BW[p_name]))]
        dict_average[st][matrix].append(np.array(arr_matrix))
        individual_matrix[p_name][st][matrix][traj][fr] = np.array(arr_matrix)
        
        if separate_output_option == 'yes':
            # write matrix for each of the frame for a given traj.
            # print ("%s/%s.%s_%s_%s_%s.csv" % (parms.outdir, p_name, st, matrix, traj, str(fr)))
            np.savetxt(("%s/%s.%s_%s_%s_%s.csv" % (parms.outdir, p_name, st, matrix, traj, str(fr))), arr_matrix, delimiter=',', fmt="%s")
        elif separate_output_option == 'no':
            pass

        progress_bar.cli_progress(dict_average, total_frames, parms.metrics, bar_length=20)
            
    return dict_average
