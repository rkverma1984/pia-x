"""
Description:
    calculate angle between helices.
    the results are highly sensitive to the sigma_cutoff provided.

Author
    Ravi kumar verma
"""

import os
from collections import OrderedDict
import sys

current_dir = os.getcwd()
sys.path.append(current_dir + '/functions')
import multi_level_dict

"""
Create the parms.states dictionary.
"""


def get_states(states, strjs):
    dict_states = multi_level_dict.mutileveldict()
    
    for i in range(0, len(states)):
        for traj in strjs[i].keys():
            if states[i]:
                dict_states[states[i]].update({traj: {}})
    return dict_states


def get_clusters(t_data, parms):
    dict_clusters = multi_level_dict.mutileveldict()
    dict_clusters.update({parms.protein_name: {}})
    
    for i in range(0, len(parms.states)):
        if parms.states[i]:
            dict_clusters[parms.protein_name].update({parms.states[i]: {}})
            for traj in t_data.trajs[i].keys():
                dict_clusters[parms.protein_name][parms.states[i]].update({traj: t_data.trajs[i].get(traj)})
    
    return dict_clusters


def get_meta_dict(t_data, parms):
    """
    Create the Meta Dictionary (meta_dict) that store data from analysis.
    """
    meta_dict = multi_level_dict.mutileveldict()
    s_dict = get_states(parms.states, t_data.trajs)
    meta_dict[parms.protein_name].update(s_dict)
    
    for st in meta_dict[parms.protein_name].keys():
        for traj in meta_dict[parms.protein_name][st].keys():
            for matrix in parms.metrics:
                meta_dict[parms.protein_name][st][traj].update({matrix: OrderedDict()})
                if matrix == parms.matrix1 or matrix == parms.matrix2:
                    for t1 in parms.helix_names[parms.protein_name]:
                        for t2 in parms.helix_names[parms.protein_name]:
                            meta_dict[parms.protein_name][st][traj][matrix].update({t1 + '_' + t2: OrderedDict()})
                
                if matrix == parms.matrix3:
                    for b1 in parms.binding_BW[parms.protein_name]:
                        for b2 in parms.binding_BW[parms.protein_name]:
                            meta_dict[parms.protein_name][st][traj][matrix].update({b1 + '_' + b2: OrderedDict()})
                
                if matrix == parms.matrix4:
                    for b_n in parms.binding_BW[parms.protein_name]:
                        meta_dict[parms.protein_name][st][traj][matrix].update({b_n: OrderedDict()})
    return meta_dict


def get_coord_dict(t_data, parms):
    coord_dict = multi_level_dict.mutileveldict()
    s_dict = get_states(parms.states, t_data.trajs)
    coord_dict[parms.protein_name].update(s_dict)
    
    for st in coord_dict[parms.protein_name].keys():
        for traj in coord_dict[parms.protein_name][st].keys():
            for matrix in parms.metrics:
                coord_dict[parms.protein_name][st][traj].update({matrix: OrderedDict()})
                if matrix == parms.matrix1 or matrix == parms.matrix2:
                    for tm_n in parms.helix_names[parms.protein_name]:
                        coord_dict[parms.protein_name][st][traj][matrix].update({tm_n: OrderedDict()})
                
                if matrix == parms.matrix3 or matrix == parms.matrix4:
                    for b_n in parms.binding_BW[parms.protein_name]:
                        coord_dict[parms.protein_name][st][traj][matrix].update({b_n: OrderedDict()})
    
    return coord_dict


def get_frame_dict(t_data, parms, clusters_dict):
    """
    Create cluster dictionary and put required clusters and frames in it
    """
    c_data = parms.cluster_data_file
    dict_frames = multi_level_dict.mutileveldict()
    for i in range(0, len(parms.states)):
        for traj in t_data.trajs[i].keys():
            if parms.states[i]:
                dict_frames[parms.protein_name][parms.states[i]].update({traj: []})
    
    """
    Read the cluster_data_file file and accordingly append the frames required in the dict_frames dictionary for each of the traj.
    """
    fh = open(c_data, 'r').readlines()
    # fh = open(os.path.join(launch_dir, c_data), 'r').readlines()
    for line in fh:
        line = line.rstrip("\n")
        if "-" not in line[:1] and 'mclstr' not in line and line != '':
            list_line = line.split(' ')
            # print list_line
            cluster_id = list_line[0]
            t_name = list_line[1]
            fr = int(list_line[-1])
            for st in dict_frames[parms.protein_name].keys():
                for traj in dict_frames[parms.protein_name][st].keys():
                    
                    if traj == t_name:
                        for clstr in clusters_dict[parms.protein_name][st][traj]:  # check if cluster_id exists in clusters_dict.
                            # print (clstr, cluster_id)
                            # print (int(clstr), int(cluster_id))
                            if clstr == cluster_id:
                                # if cluster id exists append the frames to dict_frames dictionary.
                                dict_frames[parms.protein_name][st][traj].append(fr)
                            else:
                                pass
                    dict_frames[parms.protein_name][st][traj] = sorted(dict_frames[parms.protein_name][st][traj])
    # print ('dict_frames', dict_frames)
    return dict_frames


def get_average_dict(parms):
    """
    Create cluster dictionary and put required clusters and frames in it
    """
    
    dict_average = multi_level_dict.mutileveldict()
    for i in range(0, len(parms.states)):
        dict_average.update({parms.states[i]: OrderedDict()})
        for matrix in parms.metrics:
            dict_average[parms.states[i]].update({matrix: []})
    
    return dict_average


def get_dict_individual_frames(parms, dict_frames):
    """
    Create cluster dictionary and put required clusters and frames in it
    """
    
    ind_dict = multi_level_dict.mutileveldict()
    for parms.protein_name in dict_frames.keys():
        ind_dict.update({parms.protein_name: OrderedDict()})
        for st in dict_frames.get(parms.protein_name):
            ind_dict[parms.protein_name].update({st: OrderedDict()})
            for matrix in parms.metrics:
                ind_dict[parms.protein_name][st].update({matrix: OrderedDict()})
                for traj in dict_frames[parms.protein_name][st].keys():
                    ind_dict[parms.protein_name][st][matrix].update({traj: OrderedDict()})
                    for fr in dict_frames[parms.protein_name][st][traj]:
                        ind_dict[parms.protein_name][st][matrix][traj].update({fr: []})
    return ind_dict
