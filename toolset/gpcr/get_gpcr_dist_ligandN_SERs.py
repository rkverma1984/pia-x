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
def _cal_distance(u, rec_df, bwidxin):
    pairin = [bwidxin]
    pairin_resid=convert_bwidx_to_resid(pairin, rec_df)
    list_dist = []
    for ts in u.trajectory:  # iterate through all frames
        _mindist = 999
        sel_lig_N = u.select_atoms('not protein and not resname NMA CYT and type N')
        sel_OG = u.select_atoms('protein and segid PROR and resid ' + str(pairin_resid[0]) + ' and name OG')
        for i in range(len(sel_lig_N)):
            r1 = sel_lig_N[i].position - sel_OG[0].position
            _dist = np.linalg.norm(r1)
            _mindist = min(_mindist, _dist)
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

    bwidxin="5.42"
    list_dist = _cal_distance(u, rec_df, bwidxin)
    print(list_dist)
    # save to files
    with open(output_dir + '/gpcr_dist_'+str(bwidxin)+'OG_ligandN.p', 'wb') as fp:
        pickle.dump(list_dist, fp)

    bwidxin="5.43"
    list_dist = _cal_distance(u, rec_df, bwidxin)
    print(list_dist)
    # save to files
    with open(output_dir + '/gpcr_dist_'+str(bwidxin)+'OG_ligandN.p', 'wb') as fp:
        pickle.dump(list_dist, fp)

    bwidxin="5.46"
    list_dist = _cal_distance(u, rec_df, bwidxin)
    print(list_dist)
    # save to files
    with open(output_dir + '/gpcr_dist_'+str(bwidxin)+'OG_ligandN.p', 'wb') as fp:
        pickle.dump(list_dist, fp)

