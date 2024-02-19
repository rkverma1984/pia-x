"""
Description:
    wrapper script to calculate the chi1 angle for the binding site residues.
    
    Dependencies:
        chi_angle_calc


Author
    Ravi kumar verma
"""

from collections import OrderedDict

import mdtraj as md
import numpy as np

import chi_angle_calc
import progress_bar


def get_residue_name(ltraj, res):
    """
    module to get the residue name from the residue number. Use mdtraj to do that.
    """
    
    indices = ltraj.topology.select('protein and (resSeq ' + str(res) + ' and name CA)')
    for i in indices:
        identifier = str(ltraj.topology.atom(i))
    return identifier[:3]


def binding_res_chi1_angle(T, p_name, st, matrix, traj, parms, dict_frames, dict_average, individual_matrix, total_frames, separate_output_option):
    """
    DESCRIPTION

        Calculate to chi1 angles.

    ARGUMENTS
    
        T = loaded trajectory
        :param dict_average: dictionary to save data for individual state for each matrix
        :param parms: parameters file
        :param individual_matrix: dictionary to store individual matrices
        :param T: loaded trajectory
        :param p_name: protein name
        :param st: state
        :param traj: name of the trajectory
        :param matrix: matrix type
        :param dict_frames: dictionary with list of frame-ids (specified by the cluster to be used)
        :param total_frames: total no of frames on which analysis is to be done
        :param separate_output_option: write individual matrices for each frame if 'yes'.
        binding_index[p_name] = Dictionary containing residue numbering for the binding site residues.
        binding_BW[p_name]  = Dictionary containing Ballesteros-Wienstein numbering for the binding site residues.
        
    Author:
        Ravi Kumar Verma

    """
    
    list_chi1_atoms = []
    list_resname = []
    list_b_res = []
    list_b_n = []
    list_chi1_atoms_ind = []
    for i in range(0, len(parms.binding_index[p_name])):
        b_res = parms.binding_index[p_name][i]
        b_n = parms.binding_BW[p_name][i]
        res_name = get_residue_name(T, b_res)
        
        list_b_res.append(b_res)
        list_b_n.append(b_n)
        list_resname.append(res_name)
        assert isinstance(chi_angle_calc.chi1_dict.get, object)
        list_chi1_atoms.append(chi_angle_calc.chi1_dict.get(res_name))
        
        if res_name != 'GLY' and res_name != 'ALA':
            ind = T.topology.select('protein and (resSeq %s and (name  %s))' % (b_res, ' or name '.join(chi_angle_calc.chi1_dict.get(res_name))))
            list_chi1_atoms_ind.append(ind)
        else:
            list_chi1_atoms_ind.append('ND')
    
    arr = md.compute_chi1(T, periodic=True, opt=True)
    for fr in dict_frames[p_name][st][traj]:
        temp = OrderedDict()
        for i in range(0, len(list_chi1_atoms_ind)):
            if list_chi1_atoms_ind[i] != 'ND':
                for j in range(0, len(arr[0])):
                    if np.array_equiv(sorted(list_chi1_atoms_ind[i]), sorted(np.array(arr[0][j]))):
                        dict_chi = ({fr: np.rad2deg(arr[1][fr - 1][j])})
                        temp.update({list_b_n[i]: dict_chi})
            elif list_chi1_atoms_ind[i] == 'ND':
                dict_chi = ({fr: 0.0})
                temp.update({list_b_n[i]: dict_chi})
        
        concatenated_matrix = []
        for key in parms.binding_BW[p_name]:
            # append all the pairwise angles for a single frame in a single array
            # print (key, temp[key].get(fr))
            data = temp[key].get(fr)
            # data = meta_dict[p_name][st][traj][matrix][key].get(fr)
            concatenated_matrix.append(data)
        
        # create a len(parms.binding_BW[p_name]) * len(parms.binding_BW[p_name]) matrix
        for x in range(0, len(concatenated_matrix)):
            if concatenated_matrix[x] < 0:
                concatenated_matrix[x] += 360
        
        arr_matrix = [concatenated_matrix[i:i + len(parms.binding_BW[p_name])] for i in
                      np.arange(0, len(concatenated_matrix), len(parms.binding_BW[p_name]))]
        # print arr_matrix
        dict_average[st][matrix].append(np.array(arr_matrix))
        individual_matrix[p_name][st][matrix][traj][fr] = np.array(arr_matrix)
        
        if separate_output_option == 'yes':
            # write matrix for each of the frame for a given traj.
            np.savetxt(("%s/%s.%s_%s_%s_%s.csv" % (parms.outdir, p_name, st, matrix, traj, str(fr))), arr_matrix,
                       delimiter=',', fmt="%s")
        elif separate_output_option == 'no':
            pass
        
        progress_bar.cli_progress(dict_average, total_frames, parms.metrics, bar_length=20)
    
    return dict_average
