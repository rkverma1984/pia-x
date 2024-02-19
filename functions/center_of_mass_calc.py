"""	
Author
    Ravi kumar verma
"""

from collections import OrderedDict

import numpy as np

import atomic_mass  # Dictionary consisting atomic masses.


def calculate_com(ltraj, frames, indices, mass=None):
    """
    Description:
        module to calculate COMs. return an dictionary of COM with frame_id as the key.
    Arguments:
        ltraj = loaded trajectory
        Dict_frames = list of frame_ids (specified by the cluster to be used)
        indices = indices for the atoms
        :rtype: object

    """
    
    arr_com = OrderedDict()
    if mass is not None:
        for fr in frames:
            coord = np.zeros([1, 3])
            total_mass = 0
            for i in indices:
                atom_name = ltraj.topology.atom(i)
                for key in atomic_mass.dict_mass.keys():
                    if key in str(atom_name).split('-')[-1]:
                        # get the coordinates and multiply them with the atomic mass for the given atom
                        coord += ltraj.xyz[fr - 1, i, :] * 10 * atomic_mass.dict_mass.get(key)
                        total_mass += atomic_mass.dict_mass.get(key)
                    #            print fr, coord / total_mass , np.sum(ltraj.xyz[fr - 1, indices, :], axis=0 )* 10
            arr_com.update({fr: coord / total_mass})
    
    else:
        # coord = np.zeros([ltraj.n_frames, 3])
        coord = np.sum(ltraj.xyz[:, indices, :], axis=1) * 10
        
        for fr in frames:
            arr_com.update({fr: coord[fr - 1] / len(indices)})
    return arr_com
