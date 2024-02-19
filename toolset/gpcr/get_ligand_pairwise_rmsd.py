from sys import argv
from vmd import atomsel
from vmd import molecule
import pandas as pd
import numpy as np
import multiprocessing
from multiprocessing import Process, Pool

# define enviroment
import sys, os
from pathlib import Path
home = str(Path.home())
toolset_dir = home + '/repositories/pia-x/toolset'
toolset_common = toolset_dir + "/common"
toolset_gpcr = toolset_dir + "/gpcr"
sys.path.insert(0, toolset_dir)
sys.path.insert(0, toolset_common)
sys.path.insert(0, toolset_gpcr)
from common_gpcr import *
import pickle
import warnings
warnings.filterwarnings('ignore')


def _calc_pairwise_rmsd(frame_t):
    '''
    to calculate ligand pairwise RMSD
    :param frame_t: input time/frame
    :return: pariwise rmsd for the ligands
    '''
    _ligand = vmd_ligand
    ls_ligand_pairwise_rmsd = []
    ref_ligand = atomsel(_ligand, frame=frame_t)

    for t in range(frame_t + 1, nframes):
        # for t in range(1,10):
        sel_allatom.frame = t
        sel_allatom.update()

        sel_prot_sel_backbone.frame = t
        sel_prot_sel_backbone.update()

        sel_prot_R_tm_backbone.frame = t
        sel_prot_R_tm_backbone.update()

        sel_prot_A_sel_backbone.frame = t
        sel_prot_A_sel_backbone.update()

        sel_OBS_backbone.frame = t
        sel_OBS_backbone.update()

        sel_ligand.frame = t
        sel_ligand.update()

        sel_backbone.frame = t
        sel_backbone.update()

        M = sel_OBS_backbone.fit(selection=ref_OBS_backbone)
        sel_allatom.move(M)
        ls_ligand_pairwise_rmsd.append(sel_ligand.rmsd(ref_ligand))
    return ls_ligand_pairwise_rmsd


if __name__ == '__main__':

    # pdbin="/home/khlee/desmond/output/d2gi/qnp/d2gi_qnp.f8/mda_protein.pdb"
    # dcdin="/home/khlee/desmond/output/d2gi/qnp/d2gi_qnp.f8/d2gi_qnp.f8_e0-1800.protein.wrapped.dcd"
    # gpcrin="d2"
    pdbin = argv[1]
    dcdin = argv[2]
    gpcrin = argv[3]

    molid = molecule.load("pdb", pdbin)
    molecule.read(molid, "dcd", dcdin, stride=1, waitfor=-1)
    nframes = molecule.numframes(molid)  # number of frames (PDB+DCD)

    if gpcrin == "d3": bwtable = convert_Residue2TMs(gpcrin="DRD3")
    elif gpcrin == "d2": bwtable = convert_Residue2TMs(gpcrin="DRD2")
    elif gpcrin == "d1": bwtable = convert_Residue2TMs(gpcrin="D1R")
    elif gpcrin == "d4": bwtable = convert_Residue2TMs(gpcrin="D4R")
    elif gpcrin == "b2": bwtable = convert_Residue2TMs(gpcrin="B2")
    elif gpcrin == "m2": bwtable = convert_Residue2TMs(gpcrin="M2")
    else: print("gpcrin needs")

    rec_df = prep_receptor_df(molid, bwtable)
    # print(rec_df)

    # OBS_df = prep_bindingsite_df(OBS, rec_df)
    # SBP_df = prep_bindingsite_df(SBP, rec_df)
    # OBS_resid = list(OBS_df['resid'])

    #print(OBS)
    # get OBS resid
    OBS_resid = []
    for i in OBS:
        if "EL2.52" not in i:
            #print(get_key(bwtable, "TM"+str(i)))
            OBS_resid.append(get_key(bwtable, "TM"+str(i)))

    sel_allatom, sel_backbone, ref_allatom, ref_backbone, \
    sel_prot_sel, sel_prot_sel_backbone, ref_prot_sel, ref_prot_sel_backbone, \
    sel_prot_R_tm, sel_prot_R_tm_backbone, ref_prot_R_tm, ref_prot_R_tm_backbone, \
    sel_prot_A_sel, sel_prot_A_sel_backbone, ref_prot_A_sel, ref_prot_A_sel_backbone, \
    sel_OBS, sel_OBS_backbone, ref_OBS, ref_OBS_backbone, \
    sel_ligand, ref_ligand, ref_ligand_last = \
        prep_vmd_atomsel(vmd_prot_R, vmd_prot_A, vmd_prot_B, vmd_prot_C,
                         vmd_ignore_prot_A_sel,
                         vmd_add_backbone, vmd_add_heavyatoms, vmd_add_CA,
                         vmd_add_rec_tm, OBS_resid, vmd_ligand, nframes)

    # create list of variables for multiprocesses
    _ls = []
    for i in range(1, nframes):
        _ls.append(i)
    # len(_ls)
    pool = Pool(40)
    _ligand_pairwise_rmsd = pool.map(_calc_pairwise_rmsd, _ls)
    pool.close()

    # put everything into list
    _ligand_pairwise_rmsd2 = []
    for i in range(len(_ligand_pairwise_rmsd)):
        _ligand_pairwise_rmsd2 = _ligand_pairwise_rmsd2 + _ligand_pairwise_rmsd[i]

    # dump pairwise RMSD into array
    output_dir = "output"
    check_output_dir(output_dir)
    with open(output_dir + '/ligand_pairwise_rmsd.p', 'wb') as fp:
        pickle.dump(_ligand_pairwise_rmsd2, fp)
