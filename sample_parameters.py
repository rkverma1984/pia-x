"""

Description:

Contains all the parameters needed for the Pairwise interaction calculations.

Change following parameters:
    metadata_path: Path where the wrapped and concatenated traj are saved. metadata dir must have a dcd and a pdb file with
    'protein.cpptraj.wrapped.dcd' and 'cms.wrap.protein.pdb', respectively.
    outdir: Directory where outputs are directed.
    outdir_av_and_diff_mat: Directory where differential matrices are directed.
    suffix: change suffix for the dir, where traj and the pdb are saved

    stride: stride for the traj reading.
    states: provide name of the states
    
    protein_name: protein_name type

    modify the binding tuple to include binding site residues. key used for binding tuple must match protein_name type enteries.
        binding_tuple[protein_name] : dictionary containing binding site residues for a given protein_name.
        binding_BW[protein_name] : dictionary containing Ballesteros-Wienstein numbering for binding site residues.
        helix_names[protein_name] :dictionary containing names of the helical sub segments.
        helix_index[protein_name] : dictionary containing helix start and end positions.
        
    Read further comments in the code and change accordingly.

AUTHOR
    Ravi Kumar Verma
    Version date: Nov 16, 2016
    Modified on:
        Dec 24, 2016
"""

import numpy as np
import os
import sys
from collections import OrderedDict

user = os.path.expanduser('~')

current_dir = os.getcwd()
sys.path.append(current_dir)
sys.path.append(current_dir + '/modules')
sys.path.append(current_dir + '/functions')

binding_res = OrderedDict()
helix_names = OrderedDict()
h_index = OrderedDict()

################## Modify following variables ################################

metadata_path = user + '/Research/gpcr/desmond_traj_analysis/'

current_dir = os.getcwd()
outdir = current_dir + '/results_unit_set/gpcr/'  # dir to save data for each of the frame.
outdir_av_and_diff_mat = current_dir + '/results_unit_set/traj/'  # dir to save average matrices for each of the state as well as the differential
# matrices.
plot_dir = current_dir + '/results_unit_set/plots/'

suffix = '_analysis/traj/'  # suffix for the dir, where traj and the pdb are saved
stride = 2  # skip traj every 2 frames
mass = 1  # use 0 if want to calculate centroid, else fpr center of mass use 1

protein_name = 'd3r'
trajectory_suffix = 'protein*.h5'  # change *protein*.dcd if using dcd files
topology_suffix = 'cms.wrap.protein.pdb'  # suffix for topology file. topology file must be in the same dir as dcd file.

cluster_data_file = 'sample_cluster.dat'  # name of the file containing PIA derived cluster information for frames from a given trajectory.
matrix1 = 'angle'
matrix2 = 'distance'
matrix3 = 'carbon'
matrix4 = 'chi1'

metrics = [matrix1, matrix2, matrix3, matrix4]  # tuple of analysis to be carried out.

state1 = 'active'
state2 = 'inactive'

states = [state1, state2]  # states to be considered, results would be state1 -state2, here state1 - state2

sigma_cutoff = 0  # cutoff used to select helical vectors. default cutoff in pymol script is 1.5

"""
provide atom selection for calculating COM of TM sub-segments
    modes available:
        all : consider all atoms.
        all_heavy : only heavy atoms.
        CA : only CA atoms
"""

atomselection_com_TM = 'CA'

"""
Ballesteros- Weinstein numbering for binding site residues
"""
binding_BW = OrderedDict()
binding_BW[protein_name] = (
    "2.61", "2.64", "2.65", "EL1.50", "3.28", "3.29", "3.32", "3.33", "3.36", "3.37", "3.40", "4.57", "EL2.52", "5.38",
    "5.39", "5.42", "5.43", "5.46", "5.47", "6.44", "6.48", "6.51", "6.52", "6.55", "6.56", "6.58", "6.59", "7.32", "7.35",
    "7.36", "7.39", "7.40", "7.42", "7.43")

"""
Definition of binding site residues.
Required for histograms (residue wise dGs against Ballesteros- Weinstein numbering)
"""
binding_index = OrderedDict()
binding_index[protein_name] = np.array(
        [86, 89, 90, 96, 106, 107, 110, 111, 114, 115, 118, 165, 183, 188, 189, 192, 193, 196, 197, 338, 342, 345, 346, 349, 350, 352, 353, 362, 365,
         366, 369, 370, 372, 373])

"""
name of the helical sub segments.
"""

helix_names[protein_name] = ['2e', '3e', '4e', '5e', '6e', '7e', '1m', '2m', '3m', '4m', '5m', '6m', '7m']

"""
start and end positions of the helical segments.
"""

helix_index = OrderedDict()
helix_index[protein_name] = np.array(
        [[87, 91], [100, 109], [164, 170], [186, 194], [349, 354], [362, 368], [34, 42], [77, 86], [110, 116], [158, 163], [195, 202], [341, 348],
         [369, 376]])
