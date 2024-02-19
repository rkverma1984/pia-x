"""
Description:
    Function calculates the dihedral angle between.

Arguments:
    Input coordinates of four atoms p0, p1, p2, p3
    
Author
    Ravi kumar verma
"""

import numpy as np


def calc_dihedral_0_360(p0, p1, p2, p3):
    v01 = np.array(p0) - np.array(p1)
    v32 = np.array(p3) - np.array(p2)
    v12 = np.array(p1) - np.array(p2)

    # calculate vectors perpendicular to the bonds
    v0 = np.cross(v12, v01)
    v3 = np.cross(v12, v32)

    cosa = np.dot(np.squeeze(np.asarray(v0)), np.squeeze(np.asarray(v3))) / np.linalg.norm(v0) / np.linalg.norm(v3)
    torsion = np.rad2deg(np.arccos(np.clip(cosa, -1, 1)))

    if np.dot(np.cross(np.squeeze(np.asarray(v0)), np.squeeze(np.asarray(v3))), np.squeeze(np.asarray(v12))) > 0:
        torsion = -torsion

        if torsion < 0:
            torsion = torsion + 360  # convert the torsional angle to 0 to 360 scale

    return torsion
