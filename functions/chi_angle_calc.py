"""
Description:
    wrapper script to calculate the angle between the helical segments.
    
    Dependencies:
        distance_angle_calculation.py
        parameters.py
        helix_angle.py
        
Dependencies:
    dihedral

Author
    Ravi kumar verma
"""

import os
import sys
from collections import OrderedDict


current_dir = os.getcwd()
sys.path.append(current_dir)

import dihedral_calc

chi1_dict = OrderedDict()
chi1_dict = {
    'ARG': ['N', 'CA', 'CB', 'CG'],
    'ASN': ['N', 'CA', 'CB', 'CG'],
    'ASP': ['N', 'CA', 'CB', 'CG'],
    'CYS': ['N', 'CA', 'CB', 'SG'],
    'GLN': ['N', 'CA', 'CB', 'CG'],
    'GLU': ['N', 'CA', 'CB', 'CG'],
    'HIS': ['N', 'CA', 'CB', 'CG'],
    'ILE': ['N', 'CA', 'CB', 'CG1'],
    'LEU': ['N', 'CA', 'CB', 'CG'],
    'LYS': ['N', 'CA', 'CB', 'CG'],
    'MET': ['N', 'CA', 'CB', 'CG'],
    'PHE': ['N', 'CA', 'CB', 'CG'],
    'PRO': ['N', 'CA', 'CB', 'CG'],
    'SER': ['N', 'CA', 'CB', 'OG'],
    'THR': ['N', 'CA', 'CB', 'OG1'],
    'TRP': ['N', 'CA', 'CB', 'CG'],
    'TYR': ['N', 'CA', 'CB', 'CG'],
    'VAL': ['N', 'CA', 'CB', 'CG1']
}

chi2_dict = OrderedDict()
chi2_dict = {
    'ARG': ['CA', 'CB', 'CG', 'CD'], 
    'ASN': ['CA', 'CB', 'CG', 'OD1'], 
    'ASP': ['CA', 'CB', 'CG', 'OD1'], 
    'GLN': ['CA', 'CB', 'CG', 'CD'], 
    'GLU': ['CA', 'CB', 'CG', 'CD'], 
    'HIS': ['CA', 'CB', 'CG', 'ND1'], 
    'ILE': ['CA', 'CB', 'CG1', 'CD'], 
    'LEU': ['CA', 'CB', 'CG', 'CD1'], 
    'LYS': ['CA', 'CB', 'CG', 'CD'], 
    'MET': ['CA', 'CB', 'CG', 'SD'], 
    'PHE': ['CA', 'CB', 'CG', 'CD1'], 
    'PRO': ['CA', 'CB', 'CG', 'CD'], 
    'TRP': ['CA', 'CB', 'CG', 'CD1'], 
    'TYR': ['CA', 'CB', 'CG', 'CD1'],
}


def get_chi1(T, res, res_n, atom_names, frames):
    """
    Description:
        ch1 calculation module. return a dictionary with frame number as the key.
    Arguments:
        T = loaded trajectory
        res = residue id
        res_n = residue name
        atom_names =  list consisting atom names for the for atoms required for chi1 angle calculation. names are \
        derived after matching residue name with chi1_dict keys.
    """

    arr_dih = OrderedDict()
    if res_n != 'GLY' and res_n != 'ALA':
        
        ind = T.topology.select('protein and (resSeq ' + str(res) + ' and (name ' + ' or name '.join(atom_names) + '))')

        # extract coordinates for the each of the atom
        for fr in frames:
            r = T.xyz[fr - 1, ind, :]
            
            arr_dih.update({fr: dihedral_calc.calc_dihedral_0_360(r[0], r[1], r[2], r[3])})

    else:
        for fr in frames:
            arr_dih.update({fr: 0.0})
    return arr_dih
