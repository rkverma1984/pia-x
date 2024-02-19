# define enviroment
import sys
import MDAnalysis as mda
from sys import argv
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
home = str(Path.home())
toolset_dir = home + '/repositories/pia-x/toolset'
toolset_common = toolset_dir + "/common"
toolset_gpcr = toolset_dir + "/gpcr"
sys.path.insert(0, toolset_dir)
sys.path.insert(0, toolset_common)
sys.path.insert(0, toolset_gpcr)
from get_gpcr_chi1 import get_gpcr_rec_chi1
from anal_rawdata import *
from anal_interface_contact import *


def get_pia_binding_site(u, bwtable, residue_bw_indexes, sidechain=False):
    '''
    prepare dataframe for PIA binding site residues

    :param u: given u from mda
    :param residue_bw_indexes: given a binding site residues
    :return: rawdata in dataframe and PIA infromation in dictionary
    '''

    df = pd.DataFrame()

    # chi1 measurement
    print("chi1 measurement")
    residue_ids = []
    for i in range(len(residue_bw_indexes)):
        if "EL" in residue_bw_indexes[i]:
            resid = get_key(bwtable, residue_bw_indexes[i])
        else:
            bwidxin = 'TM' + residue_bw_indexes[i]
            resid = get_key(bwtable, bwidxin)
        residue_ids.append(resid)
        mda_chi1 = get_gpcr_rec_chi1(u, resid)
        df["piabs_chi1_" + residue_bw_indexes[i]] = mda_chi1
    chi1 = get_pia_bs_chi1(df, residue_bw_indexes)

    # PIA distance measurement
    print("Pre-calculate CA coordinates")
    all_seg_vec1 = []
    for i in range(len(residue_bw_indexes)):
        if sidechain:
            # try:
            #     _selseg1 = u.select_atoms('protein and segid PROR and resid ' + str(residue_ids[i]) + ' and not backbone and not type H')
            # except:
            #     _selseg1 = u.select_atoms('protein and segid PROR and resid ' + str(residue_ids[i]) + ' and name CA')
            _selseg1 = u.select_atoms('protein and segid PROR and resid ' + str(residue_ids[i]) + ' and ((not backbone and not type H) or (name CA))')
            # if len(_selseg1)==0:
            #     _selseg1 = u.select_atoms('protein and segid PROR and resid ' + str(residue_ids[i]) + ' and name CA')
        else:
            _selseg1 = u.select_atoms('protein and segid PROR and resid ' + str(residue_ids[i]) + ' and name CA')
        _seg_vec = []
        for ts in u.trajectory:
            _a = _selseg1.center_of_mass()
            #print(i, len(_selseg1),_a)
            _seg_vec.append(_a)
        all_seg_vec1.append(_seg_vec)
    all_seg_vec2 = all_seg_vec1

    print("pairwise distance between binding site residues")
    for i in range(len(residue_bw_indexes)):
        _va = all_seg_vec1[i]
        _len_va = len(_va)
        for j in range(i + 1, len(residue_bw_indexes)):
            _vb = all_seg_vec2[j]
            _len_vb = len(_vb)
            assert _len_va == _len_vb
            list_dist = []
            for k in range(_len_va):
                list_dist.append(np.linalg.norm(_va[k] - _vb[k]))
            df["piabs_cadist_" + residue_bw_indexes[i] + "_" + residue_bw_indexes[j]] = list_dist
            #print(i,j,list_dist)

    distance_matrix = get_obs_dist_matrix(df, residue_bw_indexes)

    dict_binding_site = {"label": residue_bw_indexes, "matrix": distance_matrix, "chi1": chi1}
    return df, dict_binding_site

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
    ls_observable = [OBS, OBSplus, SBP, BSI]
    ls_out_filename = ["pia_bs_obs", "pia_bs_obsplus", "pia_bs_sbp", "pia_bs_bsi"]
    for i in range(len(ls_observable)):
        residue_bw_indexes = ls_observable[i]
        out_filename = output_dir + '/' + ls_out_filename[i]
        df, dict_bs = get_pia_binding_site(u, bwtable, residue_bw_indexes, sidechain=False)
        write_df_hdf(df, out_filename)
        write_dict_pickle(dict_bs, out_filename)

    #ls_observable = [OBS]
    #ls_out_filename = ["pia_bs_obs_scm"]
    ls_observable = [OBS, OBSplus, SBP, BSI]
    ls_out_filename = ["pia_bs_obs_COM", "pia_bs_obsplus_COM", "pia_bs_sbp_COM", "pia_bs_bsi_COM"]
    # ls_observable = [BSI]
    # ls_out_filename = ["pia_bs_bsi_COM"]
    for i in range(len(ls_observable)):
        residue_bw_indexes = ls_observable[i]
        out_filename = output_dir + '/' + ls_out_filename[i]
        df, dict_bs = get_pia_binding_site(u, bwtable, residue_bw_indexes, sidechain=True)
        write_df_hdf(df, out_filename)
        write_dict_pickle(dict_bs, out_filename)