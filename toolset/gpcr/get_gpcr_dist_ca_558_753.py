from sys import argv
from vmd import atomsel
from vmd import molecule
import pandas as pd
import numpy as np
import MDAnalysis as mda
import pickle

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

from common_gpcr import *

import warnings
warnings.filterwarnings('ignore')


## distance
def _cal_distance(u, rec_df):
    pairin = ["5.58", "7.53"]
    pairin_resid=convert_bwidx_to_resid(pairin, rec_df)
    list_dist=[]
    for ts in u.trajectory:                # iterate through all frames
        selA = u.select_atoms('protein and segid PROR and resid '+str(pairin_resid[0])+' and name CA')
        selB  = u.select_atoms('protein and segid PROR and resid '+str(pairin_resid[1])+' and name CA')
        r = selA[0].position - selB[0].position  # distance vector from atom positions
        list_dist.append(np.linalg.norm(r))
    return list_dist


if __name__ == '__main__':

    # argvs
    pdbin = argv[1]
    dcdin = argv[2]
    gpcrin = argv[3]

    # check output
    output_dir = "output"
    check_output_dir(output_dir)

    #
    if gpcrin == "d3": bwtable = convert_Residue2TMs(gpcrin="DRD3")
    elif gpcrin == "d2": bwtable = convert_Residue2TMs(gpcrin="DRD2")
    elif gpcrin == "d1": bwtable = convert_Residue2TMs(gpcrin="D1R")
    elif gpcrin == "d4": bwtable = convert_Residue2TMs(gpcrin="D4R")
    elif gpcrin == "b2": bwtable = convert_Residue2TMs(gpcrin="B2")
    elif gpcrin == "m2": bwtable = convert_Residue2TMs(gpcrin="M2")
    else: print("gpcrin needs")


    u = mda.Universe(pdbin, dcdin)
    molid = molecule.load("pdb", pdbin)
    rec_df = prep_receptor_df(molid, bwtable)

    list_dist = _cal_distance(u, rec_df)
    print(list_dist)

    # save to files
    with open(output_dir + '/gpcr_dist_ca_5.58_7.53.p', 'wb') as fp:
        pickle.dump(list_dist, fp)



