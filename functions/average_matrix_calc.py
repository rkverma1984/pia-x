"""
Description:
    wrapper script to calculate average matrices
    
    Dependencies:
        None
Author
    Ravi kumar verma
"""

import os
import sys
from collections import OrderedDict
import numpy as np

current_dir = os.getcwd()
sys.path.append(current_dir)


def av_matrix(p_name, metrics, states, dict_average, outdir_av_and_diff_mat):
    """
    DESCRIPTION

        Calculate to average matrices for each state.

    ARGUMENTS
        :type p_name : protein name
        :type metrics : list of metrics
        :type states  : list of states
        :type dict_average : dictionary with {st:{matrix:{calculated matrix table}}}
        :type outdir_av_and_diff_mat : dir where output average metrics to be stored:
    Author:
        Ravi Kumar Verma
    """
    
    average_matrix = OrderedDict()
    print("\n\n")
    for st in states:
        if st:
            average_matrix.update({st: OrderedDict()})
            cumulative_matrix = []
            for matrix in metrics:
                # print matrix
                average_matrix[st].update({matrix: []})
                for i in range(0, len(dict_average[st][matrix])):
                    data = dict_average[st][matrix][i]
                    # print data
                    if i == 0:
                        cumulative_matrix = data
                    else:
                        cumulative_matrix = (cumulative_matrix + data)
                average_matrix[st][matrix] = cumulative_matrix / len(dict_average[st][matrix])
                print(st, matrix, ", frames used :", len(dict_average[st][matrix]))
                np.savetxt(outdir_av_and_diff_mat + p_name + '.' + st + '_' + matrix + '.csv', average_matrix[st][matrix],
                           delimiter=',', fmt="%s")
