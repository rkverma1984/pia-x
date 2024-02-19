"""
Description:
    wrapper script to calculate the angle between the helical segments.
    
    Dependencies:
        distance_angle_calc.py
        helix_angle_calc.py


Author
    Ravi kumar verma
"""

from collections import OrderedDict

import numpy as np

import distance_angle_calc
import helix_angle_calc
import progress_bar


def sub_seg_angle(ltraj, p_name, st, matrix, traj, parms, coord_dict, dict_frames, dict_average, individual_matrix,
                  total_frames, separate_output_option):
    """
    DESCRIPTION

        Calculate to helical vector.

    ARGUMENTS
        :type individual_matrix: object
        :param coord_dict: dictionary to store coordinates or COM data.
        :param ltraj: = loaded trajectory
        :param dict_average: dictionary to save data for individual state for each matrix
        :param parms: parameters file
        :param individual_matrix: dictionary to store individual matrices
        :param p_name: protein name
        :param st: state
        :param traj: name of the trajectory
        :param matrix: matrix type
        :param dict_frames: dictionary with list of frame-ids (specified by the cluster to be used)
        :param total_frames: total no of frames on which analysis is to be done
        :param separate_output_option: write individual matrices for each frame if 'yes'.
        binding_index[p_name] = Dictionary containing residue numbering for the binding site residues.
        binding_BW[p_name]  = Dictionary containing Ballesteros-Wienstein numbering for the binding site residues

    Author:
        Ravi Kumar Verma

    """
    
    # print ltraj, p_name, st, matrix, traj, parms, coord_dict
    for i in range(0, len(parms.helix_index[p_name])):
        tm = parms.helix_index[p_name][i]
        tm_n = parms.helix_names[p_name][i]
        # print tm, tm_n
        # Calculate the helical vector for each of the helical sub-segment.
        coord_dict[p_name][st][traj][matrix][tm_n].update(helix_angle_calc.angle_between_helices(ltraj, int(tm[0]), int(tm[1]), dict_frames[p_name][
            st][traj], parms.sigma_cutoff))

    for fr in dict_frames[p_name][st][traj]:
        temp = OrderedDict()
        arr_done = []
        sorted_keys = []
        for i in range(0, len(parms.helix_names[p_name])):
            t1 = parms.helix_names[p_name][i]
            for j in range(0, len(parms.helix_names[p_name])):
                t2 = parms.helix_names[p_name][j]
                sorted_keys.append(t1 + '_' + t2)

                assert isinstance(t1, object)
                assert isinstance(t2, object)

                dict_angle = OrderedDict()
                if t1 != t2:
                    if t1 + '_' + t2 not in arr_done and t2 + '_' + t1 not in arr_done:
                        # calculate the pairwise distances between sub segments in a single frame and store in the dictionary.
                        dict_angle.update({fr: distance_angle_calc.get_angle(
                            coord_dict[p_name][st][traj][matrix][t1].get(fr),
                            coord_dict[p_name][st][traj][matrix][t2].get(fr))})
                        arr_done.append(t1 + '_' + t2)
                        temp.update({t1 + '_' + t2: dict_angle})
                        temp.update({t2 + '_' + t1: dict_angle})

                elif t1 == t2:
                    dict_angle.update({fr: 0.0})
                    arr_done.append(t1 + '_' + t2)
                    temp.update({t1 + '_' + t2: dict_angle})
                    temp.update({t2 + '_' + t1: dict_angle})

        concatenated_matrix = []
        for key in sorted_keys:
            # append all the pairwise angles for a single frame in a single array
            data = distance_angle_calc.convert_angle(temp[key].get(fr))
            # data = distance_angle_calc.convert_angle(meta_dict[p_name][st][traj][matrix][key].get(fr))
            concatenated_matrix.append(data)
        # create a len(parms.helix_names[p_name]) * len(parms.helix_names[p_name]) matrix
        arr_matrix = [concatenated_matrix[i:i + len(parms.helix_names[p_name])] for i in
                      np.arange(0, len(concatenated_matrix), len(parms.helix_names[p_name]))]
        dict_average[st][matrix].append(np.array(arr_matrix))

        individual_matrix[p_name][st][matrix][traj][fr] = np.array(arr_matrix)

        if separate_output_option == 'yes':
            # write matrix for each of the frame for a given traj.
            np.savetxt(("%s/%s.%s_%s_%s_%s.csv" % (parms.outdir, p_name, st, matrix, traj, str(fr))),
                       arr_matrix, delimiter=',', fmt="%s")
        elif separate_output_option == 'no':
            pass

        progress_bar.cli_progress(dict_average, total_frames, parms.metrics, bar_length=20)

    # print ('test', dict_average)
    return dict_average
