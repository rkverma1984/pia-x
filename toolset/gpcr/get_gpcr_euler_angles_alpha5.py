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

def get_gpcr_euler_angles(u, rec_segments, rec_chainid, ga_chainid):

    sel_R = ""
    for i in range(len(rec_segments)):
        if i==0:
            sel_R=sel_R+"resid "+rec_segments[i]
        else:
            sel_R=sel_R+" or resid "+rec_segments[i]
    #print(sel)
    sel_alpha5 = ""
    for i in range(len(ga_segments)):
        if i==0:
            sel_alpha5=sel_alpha5+"resid "+ga_segments[i]
        else:
            sel_alpha5=sel_alpha5+" or resid "+ga_segments[i]

    euler_series=[]
    for ts in u.trajectory:
    #for ts in range(1000):

        # atomselection
        segR_CA = u.select_atoms('protein and segid C'+rec_chainid+' and name CA and '+sel_R)
        # moment of inertia
        #segR_I = segR_CA.moment_of_inertia()
        # principal_axes
        segR_UT = segR_CA.principal_axes()

        # atomselection
        segA_CA = u.select_atoms('protein and segid C'+ga_chainid+' and name CA and '+sel_alpha5)
        # moment of inertia
        #segA_I = segA_CA.moment_of_inertia()
        # principal_axes
        segA_UT = segA_CA.principal_axes()

        Va=segR_UT[0]+segR_UT[1]+segR_UT[2]
        Vb=segA_UT[0]+segA_UT[1]+segA_UT[2]
        rotM = rot(Va,Vb)
        rotM2=R.from_matrix(rotM)
        euler_series.append(rotM2.as_euler(euler_convention, degrees=True).tolist())

    df = pd.DataFrame(euler_series, columns=['alpha', 'beta', 'gamma'])

    return df


if __name__ == '__main__':

    # argvs
    pdbin=argv[1]
    dcdin=argv[2]
    gpcrin = argv[3]

    # check output
    output_dir="output"
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

    _df=get_gpcr_euler_angles(u, rec_segments, rec_chainid, ga_chainid)
    _df.to_hdf(output_dir+'/gpcr_euler_alpha5.hdf','df', mode='w')
