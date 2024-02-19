from sys import argv
from vmd import atomsel
from vmd import molecule
import pandas as pd
import numpy as np
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

from check import *
from common_gpcr import *


if __name__ == '__main__':

    # argvs
    pdbin = argv[1]
    dcdin = argv[2]
    refpdbin = argv[3]
    gpcrin = argv[4]

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


    # load pdb/dcd
    molid = molecule.load("pdb", pdbin)
    rec_df = prep_receptor_df(molid, bwtable)

    NPxxY_bwidx = ["7.49", "7.50", "7.51", "7.52", "7.53"]
    NPxxY_resid = convert_bwidx_to_resid(NPxxY_bwidx, rec_df)

    _str = ""
    for i in range(len(NPxxY_resid)):
        _str = _str + str(NPxxY_resid[i]) + " "

    molid = molecule.load("pdb", refpdbin)
    NPxxY_ref = atomsel('protein and chain R and backbone and resid '+_str)

    molid = molecule.load("pdb", pdbin)
    molecule.read(molid, "dcd", dcdin, stride=1, waitfor=-1)
    nframes = molecule.numframes(molid)  # number of frames (PDB+DCD)
    sel_NPxxY = atomsel('protein and chain R and backbone and resid '+_str)
    sel_all = atomsel('all')
    _list=[]
    for t in range(1, nframes):
        sel_NPxxY.frame=t
        sel_NPxxY.update()
        sel_all.frame=t
        sel_all.update()
        M = sel_NPxxY.fit(selection=NPxxY_ref)
        sel_all.move(M)
        _list.append(sel_NPxxY.rmsd(NPxxY_ref))

    # save to files
    with open(output_dir + '/gpcr_rmsdbb_NPxxY.p', 'wb') as fp:
        pickle.dump(_list, fp)
