import glob
import numpy as np
from collections import OrderedDict


def report_differences(matrix, fr, difference, cutoff):
    """
    function to calculate values above cutoff.
    """
    check_matrix = difference[np.where(difference > cutoff)]
    
    if len(check_matrix) == 0:
        print("%s    \t-->\t%s: %s\t-->\t%s" % (matrix, 'frame_id', fr, 'no significant difference'))
    else:
        print("ERROR :%s\t-->\t%s, %s, %s" % (matrix, 'frame_id', fr, check_matrix))


def report_differences_average_matrix(matrix, diff, cutoff):
    """
    function to calculate values above cutoff.
    """
    check_matrix = diff[np.where(diff > cutoff)]
    
    if len(check_matrix) == 0:
        print("%s    \t-->\t%s" % (matrix, 'no significant difference'))
    else:
        print("ERROR :%s\t-->\t%s" % (matrix, check_matrix))


def compare_individual_per_frame_matrix(p_name, dict_frames, parameters):
    """
    module to calculate differences between individual metrics generated using pymol and mdtraj
    """
    arr_dir_indi = [parameters.outdir, parameters.outdir_pymol]
    
    for matrix in parameters.metrics:
        arr = OrderedDict()
        for st in dict_frames[p_name].keys():
            print("\n#######", 'individual frame:', st, "#########")
            for traj in dict_frames[p_name][st].keys():
                for fr in dict_frames[p_name][st][traj]:
                    arr.update({fr: []})
                    for dir_name in arr_dir_indi:
                        arr_csv_files = glob.glob("%s/*.csv" % (dir_name))
                        for a in arr_csv_files:
                            if p_name in a and st in a and '_' + str(fr) + '.csv' in a and matrix in a:
                                arr[fr].append(a)
            for fr in arr.keys():
                difference = []
                for i in range(0, len(arr[fr])):
                    data = np.loadtxt(arr[fr][i], delimiter=",", dtype=float)
                    if i == 0:
                        difference = data
                    elif i == 1:
                        difference -= data
                difference = np.around(difference, decimals=3)
                
                if matrix == parameters.matrix2 or matrix == parameters.matrix3:
                    # if difference in distance or carbon matrix is greater than 0.001 ang, report it.
                    report_differences(matrix, fr, difference, cutoff=0.001)
                if matrix == parameters.matrix1 or matrix == parameters.matrix4:
                    # if difference in angles or chi1 matrix is greater than 0.1 degree, report it.
                    report_differences(matrix, fr, difference, cutoff=0.1)


def compare_average_matrix(p_name, dict_frames, parms):
    """
    module to compare difference between average metrics generated using  pymol and mdtraj
    """
    
    arr_dir_av_diff = [parms.outdir_av_and_diff_mat, parms.outdir_av_and_diff_mat_pymol]
    
    for st in dict_frames[p_name].keys():
        print("\n#######", 'average matrix:', st, "#########")
        for matrix in parms.metrics:
            # print st, matrix
            arr = []
            
            for dir_name in arr_dir_av_diff:
                arr_csv_files = glob.glob("%s/%s*%s*%s.csv" % (dir_name, p_name, st, matrix))
                for fn in arr_csv_files:
                    if '_vs_' not in fn:
                        arr.append(fn)
            difference = np.around(np.loadtxt(arr[0], delimiter=",", dtype=float) - np.loadtxt(arr[1], delimiter=",",
                                                                                               dtype=float), decimals=3)
            if matrix == parms.matrix2 or matrix == parms.matrix3:
                # if difference in distance or carbon matrix is greater than 0.001 ang, report it.
                report_differences_average_matrix(matrix, difference, cutoff=0.001)
            if matrix == parms.matrix1 or matrix == parms.matrix4:
                # if difference in angles or chi1 matrix is greater than 0.1 degree, report it.
                report_differences_average_matrix(matrix, difference, cutoff=0.1)


def compare_differential_matrix(parms):
    """
    module to compare difference between differential metrics generated using  pymol and mdtraj
    """
    arr_dir_av_diff = [parms.outdir_av_and_diff_mat, parms.outdir_av_and_diff_mat_pymol]
    
    print("\n#######", 'differential matrix:', "#########")
    
    for matrix in parms.metrics:
        arr = []
        for dir_name in arr_dir_av_diff:
            arr_csv_files = glob.glob("%s/*vs*%s.csv" % (dir_name, matrix))
            arr.append(arr_csv_files[0])
        difference = np.around(np.loadtxt(arr[0], delimiter=",", dtype=float) - np.loadtxt(arr[1], delimiter=",", dtype=float),
                               decimals=3)
        
        if matrix == parms.matrix2 or matrix == parms.matrix3:
            # if difference in distance or carbon matrix is greater than 0.001 ang, report it.
            report_differences_average_matrix(matrix, difference, cutoff=0.001)
        if matrix == parms.matrix1 or matrix == parms.matrix4:
            # if difference in angles or chi1 matrix is greater than 0.1 degree, report it.
            report_differences_average_matrix(matrix, difference, cutoff=0.1)

# TODO: Write more documentation
