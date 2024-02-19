from sys import argv
# from vmd import atomsel
from vmd import molecule
# import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
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

# distance
def _cal_distance_hhdist(u, bwtable):
    pairin = ["6.50", "7.38"]
    # pairin_resid=convert_bwidx_to_resid(pairin, rec_df)
    # convert_bwidx_to_resid(pairin, rec_df)
    pairin_resid = []
    for ipairin in pairin:
        print(ipairin, get_key(bwtable,"TM"+ipairin))
        pairin_resid.append(get_key(bwtable,"TM"+ipairin))
    list_dist=[]
    for ts in u.trajectory:                # iterate through all frames
        selA = u.select_atoms('protein and segid PROR and resid '+str(pairin_resid[0])+' and not (name H*)')
        print(selA)
        print(pairin_resid[1])
        selB = u.select_atoms('protein and segid PROR and resid '+str(pairin_resid[1])+' and not (name H*)')
        print(selB)
        dist = round(np.amin(contacts.distance_array(selA.positions, selB.positions)), 2)
        list_dist.append(dist)
    return list_dist


def _cal_distance_cadist(u, rec_df):
    pairin = ["6.50", "7.38"]
    # pairin_resid=convert_bwidx_to_resid(pairin, rec_df)
    # convert_bwidx_to_resid(pairin, rec_df)
    pairin_resid = []
    for ipairin in pairin:
        pairin_resid.append(get_key(bwtable,"TM"+ipairin))
    list_dist=[]
    for ts in u.trajectory:                # iterate through all frames
        selA = u.select_atoms('protein and segid PROR and resid '+str(pairin_resid[0])+' and name CA')
        selB = u.select_atoms('protein and segid PROR and resid '+str(pairin_resid[1])+' and name CA')
        dist = round(np.amin(contacts.distance_array(selA.positions, selB.positions)), 2)
        list_dist.append(dist)
    return list_dist

def _cal_distance_prodist(u, rec_df):
    pairin = ["6.50", "7.38"]
    # pairin_resid=convert_bwidx_to_resid(pairin, rec_df)
    # convert_bwidx_to_resid(pairin, rec_df)
    pairin_resid = []
    for ipairin in pairin:
        pairin_resid.append(get_key(bwtable,"TM"+ipairin))
    list_dist=[]
    for ts in u.trajectory:                # iterate through all frames
        selA = u.select_atoms('protein and segid PROR and resid '+str(pairin_resid[0])+' and (name CA or name O or name N or name C)')
        selB = u.select_atoms('protein and segid PROR and resid '+str(pairin_resid[1])+' and not (name H*)')
        dist = round(np.amin(contacts.distance_array(selA.positions, selB.positions)), 2)
        list_dist.append(dist)
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
    #rec_df = prep_receptor_df(molid, bwtable)


    #list_dist = _cal_distance_hhdist(u, rec_df)
    list_dist = _cal_distance_hhdist(u, bwtable)
    print(list_dist)

    # save to files
    with open(output_dir + '/gpcr_hhdist_6.50_7.38.p', 'wb') as fp:
        pickle.dump(list_dist, fp)


    # list_dist = _cal_distance_cadist(u, rec_df)
    list_dist = _cal_distance_cadist(u, bwtable)
    print(list_dist)

    # save to files
    with open(output_dir + '/gpcr_cadist_6.50_7.38.p', 'wb') as fp:
        pickle.dump(list_dist, fp)


    #list_dist = _cal_distance_prodist(u, rec_df)
    list_dist = _cal_distance_prodist(u, bwtable)
    print(list_dist)

    # save to files
    with open(output_dir + '/gpcr_prodist_6.50_7.38.p', 'wb') as fp:
        pickle.dump(list_dist, fp)

