from __future__ import print_function

"""
PIA_GPCR_mdtraj script

DESCRIPTION:
    perform pairwise interaction analysis on trajectories using mdtraj
    
    Calculate four pairwise metrics.
        Distances between the helical sub-segments.
        Angle between the helical sub-segments.
        Distances between the CA atoms of the binding site residues.
        Chi1 angles for the binding site residues.
    Dependencies:
        mdtraj
        angle_sub_seg_wrapper    : in functions dir
        distance_sub_seg_wrapper : in functions dir
        chi1_angle_wrapper       : in functions dir
        ca_distance_wrapper      : in functions dir
        
     Inputs:
        parameters_file          : a file containing all the parameters needed for the Pairwise interaction calculations.
                                   Visit sample_parameters.py to have a look on the format.
        traj_data_file           : a file consisting information about trajectory names and cluster_ids needed to be included in analysis.
                                   Visit sample_traj_data.py to have a look on the format.
    
    Results:
        average(state1) - average(state2)
            
Author:
    Ravi Kumar Verma

New features:
    added functionality to calculate center of mass. Use mass=1 to consider center of masses. Else, if just centroid is needed use mass=None
    added functionality to provide atom selection for COM calculation in parameters file.
    optimized the calculation modules to speed up calculations
    added functionality to specify if the intermediate files are needed or not.
    unit test is integrated in the same script and can be accessed using:
        pia_gpcr -d unit_test [Arguments]...
    integrated help functionality in the code. Use pia_gpcr -h to see the usage and arguments.
    

    speed improvement:
        angle calculation: 2x
        distance calculation: 2x
        chi angle calculation: 5x
        carbon calculations: 2x

Modification date:
    12/24/2017
 
** New **
    1. Added functionality to create a logfile
    2. Now work
    3. Added functionality to perform calculations on pdb files as well.
        use trajectory_suffix='pdb' in parameters file to use this.
        use "PDB_ch_0" format in traj_data_file.py. number denotes the Chain ID index.
    4: include in parameters file
        mass=1 if you want to calculate center of mass, or
        mass=0 if centroid is to be calculated.
        
        if "atomselection_com_TM" in parameters file is CA both mass=1 or mass=0 gives same results.
    


    


Modification date:
    01/20/2018
"""
import warnings

import os
import sys
import glob
import getopt
import mdtraj as md
import numpy as np

# import analysis modules from modules dir.

import angle_sub_seg_wrapper
import distance_sub_seg_wrapper
import chi1_angle_wrapper
import ca_distance_wrapper
import get_dictionary
import average_matrix_calc
import differential_matrix_calc
import matrix_comparison
import graph_plotter

warnings.filterwarnings('ignore')

os.system('date')

launch_dir = os.getcwd()
print('\nlaunch_dir =', launch_dir)


# get parms and trajectory data file from system arguments #


def available_colors():
    """
    list of available colors
    :return: list
    """
    av_colors = ["Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r, GnBu",
                 "GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1",
                 "Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples",
                 "Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1, Set1_r",
                 "Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r",
                 "YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r, brg, brg_r, bwr, bwr_r, cool",
                 "cool_r, coolwarm, coolwarm_r, copper, copper_r, cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r",
                 "gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r, gist_stern",
                 "gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r",
                 "inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, ocean, ocean_r, pink, pink_r, plasma",
                 "plasma_r, prism, prism_r, rainbow, rainbow_r, seismic, seismic_r, spectral, spectral_r, spring, spring_r, summer, summer_r",
                 "terrain, terrain_r, viridis, viridis_r, winter, winter_r"]
    
    return av_colors


def show_help():
    """
    function to print out help information
    """
    print("\nUsage format:\n\tpia_gpcr -p sample_parameters.py -t sample_traj_data.py [Arguments]...", '\n')
    print("To run unit-test do:\n\tpia_gpcr -d yes [Arguments]...", '\n')
    print("Author:\tRavi Kumar Verma")
    print("Arguments")
    print("\t-p, --pfile                                Parameters file name")
    print("\t-t, --tfile                                Trajectory file name")
    print("\t-d, --do=yes|no                            Run unit test")
    print("\t-w, --write-inter=yes|no, default yes      If \"yes\" write the metrics into .npy files.")
    print("\t                                           The .npy files stores traj and frame")
    print("\t                                           ids for each of the individual matrix.")
    print("\t-o, --out-plt=yes|no|plots-only            Produce differential graphs")
    print("\t              , default yes                if \"only\", graphs are generated using previously")
    print("\t                                           calculated differential metrics")
    print("\t-c, --color, default RdBu                  specify which colors to use in graphs")
    print("\t-l, --logname, default run.log             specify logfile name")
    print("\nAvailable colors:\n\n", ", ".join(available_colors()))
    print("\n")


def main(argv):
    if argv:
        parameters_file = ''
        so = 'no'
        wo = 'yes'
        gpl = 'yes'
        clr = 'RdBu'
        
        traj_data_file = ''
        do_ut = ''
        path = ''
        log = 'run.log'
        
        try:
            opts, argument = getopt.getopt(argv, "hp:t:d:w:o:c:l:", ["pfile=", "tfile=", "do=", "write-inter=", "out-plt=", "color=", "logname"])
        except getopt.GetoptError:
            print('\nERROR:\n    Check your arguments', '\n')
            show_help()
            sys.exit(2)
        for opt, arg in opts:
            if opt == '-h':
                show_help()
                sys.exit()
            elif opt in ("-p", "--ifile"):
                parameters_file = arg
                path = launch_dir
            elif opt in ("-t", "--ofile"):
                traj_data_file = arg
                path = launch_dir
            elif opt in ("-d", "--do"):
                os.system("mkdir -p unit_set/gpcr_mdtraj unit_set/plots unit_set/traj_mdtraj")
                parameters_file = 'parameters_unit_test.py'
                traj_data_file = 'traj_unit_test.py'
                log = 'unit_test_run.log'
                path = os.path.join(launch_dir, 'unit_set/inputs/')
                do_ut = arg
                so = 'yes'
            elif opt in ("-w", "--write-inter"):
                wo = arg
            elif opt in ("-o", "--out-plt"):
                gpl = arg
            elif opt in ("-c", "--color"):
                clr = arg
            elif opt in ("-l", "--logname"):
                log = arg
        
        print(clr)
        print('\nInput file :', parameters_file)
        print('Output file is :', traj_data_file, '\n')
        return [parameters_file, traj_data_file, so, wo, do_ut, gpl, clr, log, path]
    
    else:
        show_help()
        sys.exit()


def get_total_frames(dict_frames):
    tot_fr = 0
    for p_name in dict_frames:
        for st in dict_frames[p_name]:
            for traj in dict_frames[p_name][st]:
                tot_fr += len(dict_frames[p_name][st].get(traj))
    return tot_fr


def write_intermediate_outputs(ind_mat, parms):
    for p_name in ind_mat:
        for st in ind_mat[p_name]:
            for matrix_name in ind_mat[p_name][st]:
                # write matrix into a single file
                np.save(("%s/%s.%s_%s_%s.csv" % (parms.outdir, p_name, st, matrix_name, 'all')), ind_mat[p_name][st][matrix_name])


def perform_analysis(coord_dict, dict_frames, dict_average, individual_matrix, total_fr, parms, ut_option, sep_file_option,
                     wr_int_files_option, plot_option, color_id, outf):
    
    # create output directories if not exists.
    if not os.path.isdir(parms.outdir):
        os.mkdir(parms.outdir)
    
    if not os.path.isdir(parms.outdir_av_and_diff_mat):
        os.mkdir(parms.outdir_av_and_diff_mat)

    if not os.path.isdir(parms.plot_dir):
        os.mkdir(parms.plot_dir)

    if plot_option != 'plots-only':
        """
        Perform pairwise analysis
        """
        for st in parms.states:
            for p_name in dict_frames.keys():
                for traj in dict_frames[p_name][st].keys():
                    if os.path.exists((os.path.join(parms.metadata_path, traj + parms.suffix))):
                        os.chdir(os.path.join(parms.metadata_path, traj + parms.suffix))
                    else:
                        os.chdir(os.path.join(parms.metadata_path))
                    
                    # os.system('ls')
                    if parameters.trajectory_suffix != 'pdb':
                        arr_pdb = glob.glob(traj + '*' + parms.topology_suffix)  # list pdb files
                        arr_traj = glob.glob(traj + '*' + parms.trajectory_suffix)  # list of traj files
                        # print (arr_traj[0])
                        loaded_traj = md.load(arr_traj[0], stride=parms.stride, top=arr_pdb[0])  # load the trajectories
                        
                        indices = loaded_traj.topology.select('protein and not resname ACE and not resname NME')
                        loaded_traj = loaded_traj.atom_slice(indices, inplace=False)
                    
                    else:
                        pdb_name, chainid = traj.split('_ch_')
                        
                        arr_pdb = glob.glob(pdb_name + '*' + parms.topology_suffix)  # list pdb files
                        full_traj = md.load(arr_pdb[0], standard_names=False)  # load pdb files
                        
                        chain_indices = full_traj.topology.select('protein and chainid ' + chainid)
                        loaded_traj = full_traj.atom_slice(chain_indices, inplace=False)
                        
                        print('\n', arr_pdb[0], "frames in traj:", loaded_traj.n_frames, "frames used:", len(dict_frames[p_name][st][traj]))
                        outf.write("%s\t%s\t%s\n" % (arr_pdb[0], loaded_traj.n_frames, len(dict_frames[p_name][st][traj])))
                    
                    for matrix_name in parms.metrics:
                        # Angle calculation between sub segments
                        if matrix_name == parms.matrix1:
                            angle_sub_seg_wrapper.sub_seg_angle(loaded_traj, p_name, st, matrix_name, traj, parms, coord_dict, dict_frames,
                                                                dict_average, individual_matrix, total_fr, sep_file_option)
                        
                        # Distance calculation between COM of sub segments
                        if matrix_name == parms.matrix2:
                            distance_sub_seg_wrapper.sub_seg_distance(loaded_traj, p_name, st, traj, parms, matrix_name, coord_dict, dict_frames,
                                                                      dict_average, individual_matrix, total_fr, sep_file_option, mass=parms.mass)
                        
                        # Distance calculation between CA atoms in binding site residues
                        if matrix_name == parms.matrix3:
                            ca_distance_wrapper.binding_res_ca_distance(loaded_traj, p_name, st, traj, matrix_name, parms, coord_dict, dict_frames,
                                                                        dict_average, individual_matrix, total_fr, sep_file_option, mass=parms.mass)
                            
                            # Chi1 dihedral angle calculation for residues for binding site residues
                        if matrix_name == parms.matrix4:
                            chi1_angle_wrapper.binding_res_chi1_angle(loaded_traj, p_name, st, matrix_name, traj, parms, dict_frames, dict_average,
                                                                      individual_matrix, total_fr, sep_file_option)
        
        """
        module to calculate and write average and differential metrics.
        """
        average_matrix_calc.av_matrix(parms.protein_name, parms.metrics, parms.states, dict_average, parms.outdir_av_and_diff_mat)
        
        if parms.states[0] and parms.states[1]:
            differential_matrix_calc.diff_matrix(parms.protein_name, parms.metrics, parms.outdir_av_and_diff_mat, parms.state1,
                                                 parms.state2)
        
        if plot_option == 'yes':
            print("\n")
            graph_plotter.plot_pia_matrix(parms.metrics, parms, color_id)
            print("\n")
        
        """
        module to write intermediate files.
        """
        if wr_int_files_option == 'yes':
            write_intermediate_outputs(individual_matrix, parms)
        
        """
        module to compare output with pymol generated outputs.
        """
        if ut_option == 'yes':
            matrix_comparison.compare_individual_per_frame_matrix(parms.protein_name, dict_frames, parms)
            matrix_comparison.compare_average_matrix(parms.protein_name, dict_frames, parms)
            matrix_comparison.compare_differential_matrix(parms)
    
    elif plot_option == 'plots-only':
        differential_matrix_calc.diff_matrix(parms.protein_name, parms.metrics, parms.outdir_av_and_diff_mat, parms.state1,
                                             parms.state2)
        graph_plotter.plot_pia_matrix(parms.metrics, parms, color_id)


if __name__ == "__main__":
    args = main(sys.argv[1:])
    # print args[0], args[1]
    sys.path.append(args[8])
    parameters = __import__(args[0].replace('.py', ''))
    traj_data = __import__(args[1].replace('.py', ''))
    separate_output_option = args[2]
    write_intermediate_file_option = args[3]
    perform_unit_test = args[4]
    graph_plot_option = args[5]
    colors = args[6]
    log_f_name = args[7]
    
    color_options = available_colors()
    
    out = open(log_f_name, 'w')
    f_h1 = open(os.path.join(args[8], args[0]), 'r').readlines()
    f_h2 = open(os.path.join(args[8], args[1]), 'r').readlines()
    out.write("##########################PARAMETERS USED#############################%s\n\n" % '')
    
    for li in f_h1:
        out.write("%s" % li)
    
    out.write('##########################TRAJDATA USED###############################%s\n\n' % '')
    
    for li in f_h2:
        out.write("%s" % li)
    out.write("######################################################################%s\n\n" % '')
    out.write("%s\t%s\t%s\n\n" % ('trajectory', "number of frames in traj", "number of frames used"))
    
    if colors in ','.join(color_options):
        assert isinstance(traj_data, object)
        
        # create local variables from parms and traj_data files
        for k in dir(parameters):
            locals()[k] = getattr(parameters, k)
        for l in dir(traj_data):
            locals()[l] = getattr(traj_data, l)
        
        # dictionary to store coordinate data.
        coord_dictionary = get_dictionary.get_coord_dict(traj_data, parameters)
        print('coord_dictionary done')
        
        # dictionary to store cluster_id information.
        clusters_dictionary = get_dictionary.get_clusters(traj_data, parameters)
        print('clusters_dictionary done')

        print(clusters_dictionary['d3r'].keys())
        # dictionary to store frames to be included in the analysis. It is populated using clusters_dict.
        frames_dictionary = get_dictionary.get_frame_dict(traj_data, parameters, clusters_dictionary)
        print('frames_dictionary done')
        
        average_dictionary = get_dictionary.get_average_dict(parameters)
        print('average_dictionary done')
        
        individual_matrix_dictionary = get_dictionary.get_dict_individual_frames(parameters, frames_dictionary)
        print('individual_matrix_dictionary done')
        
        total_frames = get_total_frames(frames_dictionary)
        
        assert isinstance(coord_dictionary, dict)
        assert isinstance(frames_dictionary, dict)
        assert isinstance(average_dictionary, dict)
        assert isinstance(individual_matrix_dictionary, dict)
        
        perform_analysis(coord_dictionary, frames_dictionary, average_dictionary, individual_matrix_dictionary, total_frames, parameters,
                         perform_unit_test, separate_output_option, write_intermediate_file_option, graph_plot_option, colors, out)
    
    else:
        print("Colormap " + colors + " is not recognized. Use one of following values for -c:\n")
        
        for seg in color_options:
            print(seg)

print("\n")
os.system('date')
