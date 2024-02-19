"""
Description:
    wrapper script to calculate differential matrices
    
    Dependencies:
        None


Author
    Ravi kumar verma
"""

import os
import sys
import numpy as np

current_dir = os.getcwd()
sys.path.append(current_dir)


def diff_matrix(p_name, metrics, outdir_av_and_diff_mat, st1, st2):
    """
    DESCRIPTION

        Calculate to average matrices for each state.

    ARGUMENTS
        :type p_name : protein name
        :type metrics : list of metrics
        :type outdir_av_and_diff_mat : dir where output average metrics to be stored:
        :type st1: name of first state
        :type st2: name of second state
    Author:
        Ravi Kumar Verma
    """
    
    for matrix in metrics:
        matrix1 = np.loadtxt(("%s/%s.%s_%s.csv" % (outdir_av_and_diff_mat, p_name, st1, matrix)), delimiter=",",
                             dtype=float)
        matrix2 = np.loadtxt(("%s/%s.%s_%s.csv" % (outdir_av_and_diff_mat, p_name, st2, matrix)), delimiter=",",
                             dtype=float)
        
        differential_matrix = matrix1 - matrix2
        np.savetxt(outdir_av_and_diff_mat + p_name + '.' + st1 + '_vs_' + st2 + '_' + matrix + '.csv',
                   differential_matrix, delimiter=',', fmt="%s")