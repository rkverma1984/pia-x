# define enviroment
import sys
# import pandas as pd
from pathlib import Path
home = str(Path.home())
toolset_dir = home + '/repositories/pia-x/toolset'
toolset_common = toolset_dir + "/common"
toolset_gpcr = toolset_dir + "/gpcr"
sys.path.insert(0, toolset_dir)
sys.path.insert(0, toolset_common)
sys.path.insert(0, toolset_gpcr)
# from check import *
from common_gpcr import *
import MDAnalysis as mda
from sys import argv
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from get_gpcr_chi1 import get_gpcr_rec_chi1


if __name__ == '__main__':

    # check
    output_dir = "output"
    check_output_dir(output_dir)

    # argv
    pdbin = argv[1]
    dcdin = argv[2]
    gpcrin = argv[3]

    bwtable = convert_Residue2TMs(gpcrin)

    u = mda.Universe(pdbin, dcdin)
    df = pd.DataFrame()

    SBPin = SBP


    # SBP chi1
    SBPin_resids = []
    for i in range(len(SBPin)):
        if "EL" in SBPin[i]:
            resid = get_key(bwtable, SBPin[i])
        else:
            bwidxin = 'TM' + SBPin[i]
            resid = get_key(bwtable, bwidxin)
            SBPin_resids.append(resid)
            mda_chi1 = get_gpcr_rec_chi1(u, resid)
            df["sbp_chi1_" + SBPin[i]] = mda_chi1

    # SBP cadist
    # for i in range(len(SBPin)):
    #     for j in range(i + 1, len(SBPin)):
    #         pairin = [SBPin[i], SBPin[j]]
    #         list_dist = []
    #         for ts in u.trajectory:  # iterate through all frames
    #             selA = u.select_atoms('protein and segid PROR and resid ' + str(SBPin_resids[i]) + ' and name CA')
    #             selB = u.select_atoms('protein and segid PROR and resid ' + str(SBPin_resids[j]) + ' and name CA')
    #             r = selA[0].position - selB[0].position  # distance vector from atom positions
    #             list_dist.append(np.linalg.norm(r))
    #         df["sbp_cadist_" + SBPin[i] + "_" + SBPin[j]] = list_dist



    print("pre-calculate CA coor")
    all_seg_vec1 = []
    for i in range(len(SBPin)):
        _selseg1 = u.select_atoms( 'protein and segid PROR and resid ' + str(SBPin_resids[i]) + ' and name CA')
        _seg_vec = []
        for ts in u.trajectory:
            _a = _selseg1.center_of_mass()
            _seg_vec.append(_a)
        all_seg_vec1.append(_seg_vec)

    all_seg_vec2 = all_seg_vec1

    print("pairwise distance")
    all_seg_dist = []
    colnames = []
    #for i in range(len(OBSin_resids)):
    for i in range(len(SBPin)):
        _va = all_seg_vec1[i]
        _len_va = len(_va)

        #for j in range(len(OBSin_resids)):
        for j in range(i + 1, len(SBPin)):
            _vb = all_seg_vec2[j]
            _len_vb = len(_vb)
            assert _len_va == _len_vb

            list_dist = []
            for k in range(_len_va):
                list_dist.append(np.linalg.norm(_va[k]- _vb[k]))

            df["sbp_cadist_" + SBPin[i] + "_" + SBPin[j]] = list_dist








    # SBP hhdist
    # for i in range(len(SBPin)):
    #     for j in range(i + 1, len(SBPin)):
    #         pairin = [SBPin[i], SBPin[j]]
    #         list_dist = []
    #         for ts in u.trajectory:  # iterate through all frames
    #             selA = u.select_atoms('protein and segid PROR and resid ' + str(SBPin_resids[i]) + ' and name CA')
    #             selB = u.select_atoms('protein and segid PROR and resid ' + str(SBPin_resids[j]) + ' and name CA')
    #             r = selA[0].position - selB[0].position  # distance vector from atom positions
    #             list_dist.append(np.linalg.norm(r))
    #         df["sbp_hhdist_" + SBPin[i] + "_" + SBPin[j]] = list_dist

    df.to_hdf(output_dir + '/df_sbp.hdf', key='df', mode='w')
