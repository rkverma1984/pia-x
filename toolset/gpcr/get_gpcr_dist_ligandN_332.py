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
    pairin = ["3.32"]
    pairin_resid=convert_bwidx_to_resid(pairin, rec_df)
    list_dist=[]
    for ts in u.trajectory:                # iterate through all frames
        _mindist = 999
        sel_lig_N = u.select_atoms('not protein and not resname NMA CYT and type N')
        selA_OD1 = u.select_atoms('protein and segid PROR and resid ' + str(pairin_resid[0]) + ' and name OD1')
        selA_OD2 = u.select_atoms('protein and segid PROR and resid ' + str(pairin_resid[0]) + ' and name OD2')
        for i in range(len(sel_lig_N)):
            r1=sel_lig_N[i].position - selA_OD1[0].position
            r2=sel_lig_N[i].position - selA_OD2[0].position
            _dist = min(np.linalg.norm(r1),np.linalg.norm(r2))
            _mindist = min(_mindist,_dist)
        list_dist.append(_mindist)
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

    # save to files
    with open(output_dir + '/gpcr_dist_3.32_ligandN.p', 'wb') as fp:
        pickle.dump(list_dist, fp)



