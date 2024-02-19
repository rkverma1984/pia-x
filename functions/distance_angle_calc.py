"""
Description:
    functions to calculate distances and angles between vectors.
    
Author
    Ravi kumar verma
"""

import numpy as np
from collections import OrderedDict
import math


def get_distance(u, v, f):
    """
    Description:
        distance calculation module. return a dictionary with frame number as the key.
    Arguments:
        vectors u and v
        frame id
    """
    # print u, v
    dict_dist = OrderedDict()

    dist = np.linalg.norm(u - v)
    dict_dist.update({f: dist})

    return dict_dist


def get_angle(u, v):
    """
    Description:
        angle calculation module. angles are given in radians. return a dictionary with frame number as the key.
    Arguments:
        vectors a and b
        frame id
    """
    u = np.squeeze(np.asarray(u))
    v = np.squeeze(np.asarray(v))

    cosa = np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v)  # -> cosine of the angle
    angle = np.arccos(np.clip(cosa, -1, 1))
    return angle


def convert_angle(angle):
    """
    Description:
        module to convert radians to degrees.
    Arguments:
        vectors a and b
        frame id
    """
    angle_deg = math.degrees(float(angle))

    if angle_deg < 0:
        angle_deg += 360
        return angle_deg
    else:
        return angle_deg
