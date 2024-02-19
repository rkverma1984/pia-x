"""

Description:

Contains all the parameters needed for the Pairwise interaction calculations.

Change following parameters:
    metadata_path: Path where the wrapped and concatenated traj are saved. metadata dir must have a dcd and a pdb file with
    'protein.cpptraj.wrapped.dcd' and 'cms.wrap.protein.pdb', respectively.
    outdir: Directory where outputs are directed.
    outdir_av_and_diff_mat: Directory where diiferential matrices are directed.
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
        Dec 27, 2016
"""

import mdtraj as md
import numpy as np
import os
import sys
from collections import OrderedDict

from os.path import expanduser

user = os.path.expanduser('~')
current_dir = os.getcwd()

binding_res = OrderedDict()
helix_names = OrderedDict()
h_index = OrderedDict()

################## Modify following variables ################################


metadata_path = current_dir +  '/inputs/traj/'

current_dir = os.getcwd()
outdir = current_dir + '/pia_frames/'  				   			# dir to save data for each of the frame.
outdir_av_and_diff_mat = current_dir + '/pia_average/'  		# dir to save average matrices for each of the state as well as the differential matrices.
#outdir_pymol = current_dir + '/pymol_frames/'  					# dir to save data for each of the frame.
#outdir_av_and_diff_mat_pymol = current_dir + '/pymol_average/' 	# dir to save average matrices for each of the state as well as the differential matrices.
outdir_pymol = current_dir + '/vmd_frames/'  					# dir to save data for each of the frame.
outdir_av_and_diff_mat_pymol = current_dir + '/vmd_average/' 	# dir to save average matrices for each of the state as well as the differential matrices.
plot_dir = current_dir + '/pia_plots/'

suffix = '_analysis/traj/'  # suffix for the dir, where traj and the pdb are saved
stride = 1  # skip traj every 2 frames
mass = 1
protein_name = 'dat'
trajectory_suffix = 'dcd'  								# change *protein*.dcd if using dcd files
topology_suffix = 'pdb'  							# suffix for topology file. topology file must be in the same dir as dcd file.

cluster_data_file = current_dir + '/inputs/cluster_unit_test.dat'	# name of the file contaioning PIA derived cluster information for frames from a given trajectory.

matrix1 = 'angle'
matrix2 = 'distance'
matrix3 = 'carbon'
matrix4 = 'chi1'

# metrics = [matrix3, matrix4]
metrics = [matrix1, matrix2, matrix3, matrix4]  # tuple of analysis to be carried out.

state1 = 'jjc8091'
state2 = 'jjc8088'

states = [state1, state2]  # states to be considered, results would be state1 -state2, here state1 - state2

sigma_cutoff = 0  # cutoff used to select helical vectors. 1.5 value taked from pymol

"""
provide atom selection for calculating COM of TM subsegments
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
binding_BW[protein_name] = ('r0', 'r1')

"""
Definition of binding site residues.
Required for histograms (residue wise dGs against Ballesteros- Weinstein numbering)
"""
binding_index = OrderedDict()
binding_index[protein_name] = np.array([ 76,  77])

"""
name of the helical sub segments.
"""
helix_names[protein_name] = ['NT', '1i']

"""
start and end positions of the helical segments.
"""
helix_index = OrderedDict()
helix_index[protein_name] = np.array([[58, 64],[65, 74]])
