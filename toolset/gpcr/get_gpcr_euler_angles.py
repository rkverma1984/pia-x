import MDAnalysis as mda
import pandas as pd
from sys import argv

# define enviroment
import sys, os
from pathlib import Path
home = str(Path.home())
toolset_dir = home+'/repositories/pia-x/toolset'
toolset_common = toolset_dir+"/common"
toolset_gpcr = toolset_dir+"/gpcr"
sys.path.insert(0, toolset_dir)
sys.path.insert(0, toolset_common)
sys.path.insert(0, toolset_gpcr)
# end define enviroment

from check import *
from common_gpcr import *
from math_math import *


euler_convention='xyx'

def get_gpcr_euler_angles(u, rec_seg_in, ga_seg_in, rec_chainid, ga_chainid):
    '''
    :param u: mda object
    :param rec_seg_in:
    :param ga_seg_in:
    :param rec_chainid:
    :param ga_chainid:
    :return:
    '''

    sel_R = ""
    for i in range(len(rec_seg_in)):
        if i == 0:
            sel_R = sel_R+"resid "+rec_seg_in[i]
        else:
            sel_R = sel_R+" or resid "+rec_seg_in[i]
    print(sel_R)
    sel_A = ""
    for i in range(len(ga_seg_in)):
        if i == 0:
            sel_A = sel_A+"resid "+ga_seg_in[i]
        else:
            sel_A = sel_A+" or resid "+ga_seg_in[i]

    euler_series=[]
    for ts in u.trajectory:
        print(ts)
    # for ts in range(1000):

        # atomselection
        segR_CA = u.select_atoms('protein and segid PRO'+rec_chainid+' and backbone and '+sel_R)
        # moment of inertia
        # segR_I = segR_CA.moment_of_inertia()
        # principal_axes
        segR_UT = segR_CA.principal_axes()

        # atomselection
        segA_CA = u.select_atoms('protein and segid PRO'+ga_chainid+' and backbone and '+sel_A)
        # moment of inertia
        # segA_I = segA_CA.moment_of_inertia()
        # principal_axes
        segA_UT = segA_CA.principal_axes()

        Va = segR_UT[0] + segR_UT[1] + segR_UT[2]
        Vb = segA_UT[0] + segA_UT[1] + segA_UT[2]
        rotM = rot(Va, Vb)
        rotM2 = R.from_matrix(rotM)
        euler_series.append(rotM2.as_euler(euler_convention, degrees=True).tolist())

    df = pd.DataFrame(euler_series, columns=['alpha', 'beta', 'gamma'])

    return df


if __name__ == '__main__':

    # argvs
    pdbin = argv[1]
    dcdin = argv[2]
    gpcrin = argv[3]
    gain = argv[4]

    # check output
    output_dir = "output"
    check_output_dir(output_dir)

    # check gpcrin
    if gpcrin == "d3": bwtable = convert_Residue2TMs(gpcrin="DRD3")
    elif gpcrin == "d2": bwtable = convert_Residue2TMs(gpcrin="DRD2")
    elif gpcrin == "d2short": bwtable = convert_Residue2TMs(gpcrin="D2R_short")
    elif gpcrin == "d1": bwtable = convert_Residue2TMs(gpcrin="D1R")
    elif gpcrin == "d4": bwtable = convert_Residue2TMs(gpcrin="D4R")
    elif gpcrin == "b2": bwtable = convert_Residue2TMs(gpcrin="B2")
    elif gpcrin == "m2": bwtable = convert_Residue2TMs(gpcrin="M2")
    else: print("gpcrin needs")


    u = mda.Universe(pdbin, dcdin)

    #rec_segments2 = ['68:96', '104:138', '147:172', '187:219']
    # TM2.50/3.50/5.60/6.50/7.50
    #rec_segments2 = ['75', '128', '200', '344', '380']
    #rec_segments2 = ['75', '200', '344']

    gi_segments = ['320:340']
    go_segments = gi_segments

    rec_chainid = 'R'
    ga_chainid = 'A'

    if gain == "gi":
        _df = get_gpcr_euler_angles(u, rec_segments2, gi_segments, rec_chainid, ga_chainid)
        _df.to_hdf(output_dir + '/df_euler_angle.hdf', 'df', mode='w')
    elif gain == "go":
        _df = get_gpcr_euler_angles(u, rec_segments2, go_segments, rec_chainid, ga_chainid)
        _df.to_hdf(output_dir + '/df_euler_angle.hdf', 'df', mode='w')
    else:
        print("grcp_input needs")

