"""
Description:
    calculate angle between helices.
    the results are higly sensitive to the sigma_cutoff provided.

Author
    Ravi kumar verma
    
modification:
    The previous version of this script used to calculate helix vectors for individual frames, hence the slow speed. In the modified code the
    helix vectors for a given helical subsegment are calculated for the entire trajectory and then the needed frames are sliced.
"""

import mdtraj as md
import numpy as np
import os
from collections import OrderedDict
import sys
import math
from importlib import import_module

current_dir = os.getcwd()
sys.path.append(current_dir)

import distance_angle_calc


def helix_orientation(ltraj, start, end, sigma_cutoff):
    """
    DESCRIPTION

        Calculate to helical vector.

    ARGUMENTS
        ltraj = loaded trajectory
        start = int: start of the helix
        end   = int: end of the helix
        sigma_cutoff = float: drop outliers outside
        (standard_deviation * sigma_cutoff) {default: 1.5}

ADAPTED AND MODIFIED FROM:
    anglebetweenhelices.py

    """
    
    # find indices for the C and O atoms
    
    ind1 = ltraj.topology.select('protein and (resSeq ' + str(start) + ' to ' + str(end) + ' and (name C)' + ')')
    ind2 = ltraj.topology.select('protein and (resSeq ' + str(start) + ' to ' + str(end) + ' and (name O)' + ')')
    print(len(ind1), len(ind2))
    c1 = ltraj.xyz[:, ind1, :] * 10
    c2 = ltraj.xyz[:, ind2, :] * 10
    vec_list = c2 - c1  # create vectors using the backbone C and O for each of the residue.
    #    print "vec_list1 :\n", vec_list
    vec = np.sum(vec_list, axis=1)  # create an average vector for the helical segment.
    #    print "vector1 :", vec
    '''
    if only two residue in the helical segment take vec as it is, else remove the outliers and create a new vec.
    '''
    
    for n in range(0, ltraj.n_frames):
        # print n
        if len(vec_list[n]) > 2 and sigma_cutoff > 0:
            angle_list = [distance_angle_calc.get_angle(vec[n], x) for x in
                          vec_list[n]]  # calculate the angle between the average vec and each of the individual vec in the vec_list.
            angle_mu = np.mean(angle_list)
            angle_std = np.std(angle_list)
            # '''
            # compare each of the above calculated angle with angle_std * sigma_cutoff. if the angle value is greater than angle_std * sigma_cutoff,
            # exclude the outlier vector from the vec_list.
            # '''
            nvec_list = [vec_list[n][i] for i in range(len(vec_list[n])) if abs(angle_list[i] - angle_mu) < angle_std * sigma_cutoff]
            vec[n] = np.sum(nvec_list, axis=0)  # create an average vector for the helical segment after removing outliers
            # print "vec_list2 :\n", vec_list
            # print "vector2 :", vec
            vec[n] = vec[n] / np.linalg.norm(vec[n])
    return vec


def angle_between_helices(ltraj, start, end, frames, sigma_cutoff):
    """
    DESCRIPTION
        Calculates the angle between two helices

    ARGUMENTS
        ltraj = loaded trajectory
        start = int: start of the helix
        end   = int: end of the helix
        frames = an array that contains the frame ids for which the calculation is required

    ADAPTED AND MODIFIED FROM:
        anglebetweenhelices.py
    """
    arr_vec = OrderedDict()
    # for fr in frames:
    #     vector = helix_orientation(ltraj, start, end, fr, sigma_cutoff)
    #     arr_vec.update({fr: vector})
    
    vector = helix_orientation(ltraj, start, end, sigma_cutoff)
    for fr in frames:
        # print fr, vector[fr-1]
        arr_vec.update({fr: vector[fr-1]})
    
    return arr_vec
