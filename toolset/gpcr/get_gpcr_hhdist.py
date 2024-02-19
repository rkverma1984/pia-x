from natsort import natsorted
## define enviroment
import sys, os
from pathlib import Path

home = str(Path.home())
toolset_dir = home + '/repositories/pia-x/toolset'
toolset_common = toolset_dir + "/common"
toolset_gpcr = toolset_dir + "/gpcr"
sys.path.insert(0, toolset_dir)
sys.path.insert(0, toolset_common)
sys.path.insert(0, toolset_gpcr)
## end define enviroment
from gpcr_contact import *
# from mkplots import *
from plot_master import *

from sys import argv

# sys.path.append('/home/khlee/repositories/pia-x/coupler/scripts')
# from run import *
# from _heatmap import *

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
# import pandas as pd
import numpy as np
# import pickle
# import glob
# import time

# import multiprocessing
from multiprocessing import Process, Pool

# from natsort import natsorted, ns

import warnings

warnings.filterwarnings('ignore')

# from check import *
from common_gpcr import *


# conf inputs
# conf_dir='/home/khlee/repositories/dopamine/toolset/anal_d2d3gigo'
# sys.path.insert(0, conf_dir)
# from common_d2d3gigo import *

def _mp_calc_hhdist(_selA, _selB):
    timeseries = []
    for ts in u.trajectory:
        dist = round(np.amin(contacts.distance_array(_selA.positions, _selB.positions)), 3)
        timeseries.append(dist)
    return timeseries


def get_hhdist(u, _ls_resid_A, _chainid_A, _ls_resid_B, _chainid_B):
    # atom selection array for selA
    _ls_selA = []
    for i in range(len(_ls_resid_A)):
        # _sel = str("protein and segid C"+_chainid_A+" and not (name H*) and resid " + str(_ls_resid_A[i]))
        _sel = str("protein and segid PRO" + _chainid_A + " and not (name H*) and resid " + str(_ls_resid_A[i]))
        _ls_selA.append(u.select_atoms(_sel))

    # atom selection array for selB
    _ls_selB = []
    for i in range(len(_ls_resid_B)):
        # _sel = str("protein and segid C"+_chainid_B+" and not (name H*) and resid " + str(_ls_resid_B[i]))
        _sel = str("protein and segid PRO" + _chainid_B + " and not (name H*) and resid " + str(_ls_resid_B[i]))
        _ls_selB.append(u.select_atoms(_sel))

    # start = time.time()

    _ls = []
    for i in range(len(_ls_selA)): [_ls.append((_ls_selA[i], _ls_selB[j])) for j in range(len(_ls_selB))]
    pool = Pool(64)
    _results = pool.starmap(_mp_calc_hhdist, _ls)
    pool.close()
    pool.join()

    k = 0
    dict_dist_mp = {}
    # for i in range(len(_ls_selA)):
    for i in _ls_resid_A:
        # for j in range(len(_ls_selB)):
        for j in _ls_resid_B:
            name = "hhdist_" + str(_chainid_A) + "_" + str(i) + "_" + str(_chainid_B) + "_" + str(j)
            dict_dist_mp[name] = _results[k]
            k = k + 1

    # end = time.time()
    # print(end - start)

    return pd.DataFrame.from_dict(dict_dist_mp)


def get_cadist(u, _ls_resid_A, _chainid_A, _ls_resid_B, _chainid_B):
    # atom selection array for selA
    _ls_selA = []
    for i in range(len(_ls_resid_A)):
        # _sel = str("protein and segid C"+_chainid_A+" and not (name H*) and resid " + str(_ls_resid_A[i]))
        _sel = str("protein and segid PRO" + _chainid_A + " and name CA and resid " + str(_ls_resid_A[i]))
        _ls_selA.append(u.select_atoms(_sel))

    # atom selection array for selB
    _ls_selB = []
    for i in range(len(_ls_resid_B)):
        # _sel = str("protein and segid C"+_chainid_B+" and not (name H*) and resid " + str(_ls_resid_B[i]))
        _sel = str("protein and segid PRO" + _chainid_B + " and name CA and resid " + str(_ls_resid_B[i]))
        _ls_selB.append(u.select_atoms(_sel))

    # start = time.time()

    _ls = []
    for i in range(len(_ls_selA)): [_ls.append((_ls_selA[i], _ls_selB[j])) for j in range(len(_ls_selB))]
    pool = Pool(64)
    _results = pool.starmap(_mp_calc_hhdist, _ls)
    pool.close()
    pool.join()

    k = 0
    dict_dist_mp = {}
    # for i in range(len(_ls_selA)):
    for i in _ls_resid_A:
        # for j in range(len(_ls_selB)):
        for j in _ls_resid_B:
            name = "cadist_" + str(_chainid_A) + "_" + str(i) + "_" + str(_chainid_B) + "_" + str(j)
            dict_dist_mp[name] = _results[k]
            k = k + 1

    # end = time.time()
    # print(end - start)

    return pd.DataFrame.from_dict(dict_dist_mp)


if __name__ == '__main__':

    # argv
    PDBin = argv[1]
    DCDin = argv[2]
    gpcrin = argv[3]

    # load pdb/dcd
    u = mda.Universe(PDBin, DCDin)
    # check output directory
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


    # print(common_R)
    # print(common_A)
    # print(common_B)
    # print(common_R_il2)
    # print(common_A_il2)

    # inputs
    # _ls_resid_R = common_R
    # _ls_resid_A = common_A
    # _ls_resid_B = common_B
    # _chainid_R = "R"
    # _chainid_A = "A"
    # _chainid_B = "B"
    # # _df = get_hhdist(u,_ls_resid_R,_chainid_R,_ls_resid_A,_chainid_A)
    # # _df.to_hdf(output_dir+'/df_hhdist_RA.hdf', key='df', mode='w')
    # # _df = get_hhdist(u,_ls_resid_R,_chainid_R,_ls_resid_B,_chainid_B)
    # # _df.to_hdf(output_dir+'/df_hhdist_RB.hdf', key='df', mode='w')
    #
    # _df = get_cadist(u, _ls_resid_R, _chainid_R, _ls_resid_A, _chainid_A)
    # _df.to_hdf(output_dir + '/df_cadist_RA.hdf', key='df', mode='w')
    # _df = get_cadist(u, _ls_resid_R, _chainid_R, _ls_resid_B, _chainid_B)
    # _df.to_hdf(output_dir + '/df_cadist_RB.hdf', key='df', mode='w')

    # inputs
    # _ls_resid_R = common_R_il2
    # _ls_resid_A = common_A_il2
    # _chainid_R = "R"
    # _chainid_A = "A"
    # # _df = get_hhdist(u,_ls_resid_R,_chainid_R,_ls_resid_A,_chainid_A)
    # # _df.to_hdf(output_dir+'/df_hhdist_RA_il2.hdf', key='df', mode='w')
    #
    # _df = get_cadist(u, _ls_resid_R, _chainid_R, _ls_resid_A, _chainid_A)
    # _df.to_hdf(output_dir + '/df_cadist_RA_il2.hdf', key='df', mode='w')


    # pre-scanning
    ls_R_resids = []
    for ts in u.trajectory:
        sel = u.select_atoms('(around 5 (segid PROA) and not type H) and segid PROR')
        ls_R_resids = ls_R_resids + sel.residues.resids.tolist()
        ls_R_resids = natsorted(list(set(ls_R_resids)))
    #print(len(ls_R_resids))
    #ls_R_resids_around_A = ls_R_resids

    # scaning the trajectory for chain A that near RA interface
    ls_A_resids = []
    for ts in u.trajectory:
        sel = u.select_atoms('(around 5 (segid PROR) and not type H) and segid PROA')
        ls_A_resids = ls_A_resids + sel.residues.resids.tolist()
        ls_A_resids = natsorted(list(set(ls_A_resids)))
    #print(len(ls_A_resids))

    _ls_resid_R = ls_R_resids
    _ls_resid_A = ls_A_resids
    _chainid_R = "R"
    _chainid_A = "A"
    # _df = get_hhdist(u,_ls_resid_R,_chainid_R,_ls_resid_A,_chainid_A)
    # _df.to_hdf(output_dir+'/df_hhdist_RA_il2.hdf', key='df', mode='w')

    _df = get_cadist(u, _ls_resid_R, _chainid_R, _ls_resid_A, _chainid_A)
    _df.to_hdf(output_dir + '/df_cadist_RA.hdf', key='df', mode='w')
    _df = get_hhdist(u,_ls_resid_R,_chainid_R,_ls_resid_A,_chainid_A)
    _df.to_hdf(output_dir+'/df_hhdist_RA.hdf', key='df', mode='w')
