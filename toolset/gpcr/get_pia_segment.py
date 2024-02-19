from sys import argv
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
# import pickle
# import multiprocessing
import MDAnalysis as mda
import numpy as np
import math
import pandas as pd

# from pylab import cross,dot,inv
# # https://stackoverflow.com/questions/36409140/create-a-rotation-matrix-from-2-normals
# from scipy.spatial.transform import Rotation as R
# # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html
import warnings
warnings.filterwarnings('ignore')


def get_segments_distance(u, segmanet_list1, segmanet_chainid1, labels1, segmanet_list2, segmanet_chainid2, labels2, backbone_CA=True):


    if backbone_CA:  # backbone CA
        print("pre-calculate helix distance for segment1")
        all_seg_vec1 = []
        for i in range(len(segmanet_list1)):
            _selseg1 = u.select_atoms( 'protein and segid PRO' + segmanet_chainid1 + ' and resid ' + segmanet_list1[i] + ' and (not type H) and name CA')
            _seg_vec = []
            for ts in u.trajectory:
                _a = _selseg1.center_of_mass()
                _seg_vec.append(_a)
            all_seg_vec1.append(_seg_vec)

        print("pre-calculate helix distance for segment2")
        all_seg_vec2 = []
        for i in range(len(segmanet_list2)):
            _selseg1 = u.select_atoms('protein and segid PRO' + segmanet_chainid2 + ' and resid ' + segmanet_list2[i] + ' and (not type H) and name CA')
            _seg_vec = []
            for ts in u.trajectory:
                _a = _selseg1.center_of_mass()
                _seg_vec.append(_a)
            all_seg_vec2.append(_seg_vec)
    else: # side-chain COM, if not work, then using backbone CA
        print("pre-calculate helix distance for segment1")
        all_seg_vec1 = []
        for i in range(len(segmanet_list1)):
            # try:
            #     _selseg1 = u.select_atoms('protein and segid PRO' + segmanet_chainid1 + ' and resid ' + segmanet_list1[i] + ' and (not type H) and (not backbone)')
            # except:
            #     _selseg1 = u.select_atoms('protein and segid PRO' + segmanet_chainid1 + ' and resid ' + segmanet_list1[i] + ' and (not type H) and name CA')
            _selseg1 = u.select_atoms('protein and segid PRO' + segmanet_chainid1 + ' and resid ' + segmanet_list1[i] + ' and (not type H)')
            _seg_vec = []
            for ts in u.trajectory:
                _a = _selseg1.center_of_mass()
                _seg_vec.append(_a)
            all_seg_vec1.append(_seg_vec)

        print("pre-calculate helix distance for segment2")
        all_seg_vec2 = []
        for i in range(len(segmanet_list2)):
            # try:
            #     _selseg1 = u.select_atoms('protein and segid PRO' + segmanet_chainid2 + ' and resid ' + segmanet_list2[i] + ' and (not type H) and (not backbone)')
            # except:
            #     _selseg1 = u.select_atoms('protein and segid PRO' + segmanet_chainid2 + ' and resid ' + segmanet_list2[i] + ' and (not type H) and name CA')
            _selseg1 = u.select_atoms('protein and segid PRO' + segmanet_chainid2 + ' and resid ' + segmanet_list2[i] + ' and (not type H)')
            _seg_vec = []
            for ts in u.trajectory:
                _a = _selseg1.center_of_mass()
                _seg_vec.append(_a)
            all_seg_vec2.append(_seg_vec)

    print("pairwise distance")
    all_seg_dist = []
    colnames = []
    for i in range(len(segmanet_list1)):
        _va = all_seg_vec1[i]
        _len_va = len(_va)

        for j in range(len(segmanet_list2)):
            _vb = all_seg_vec2[j]
            _len_vb = len(_vb)
            assert _len_va == _len_vb

            _ls_diff_dist = []
            for k in range(_len_va):
                print(np.linalg.norm(_va[k]- _vb[k]))
                _ls_diff_dist.append(np.linalg.norm(_va[k]- _vb[k]))

            colnames.append("seg_dist_"+labels1[i] + "_" + labels2[j])
            all_seg_dist.append(_ls_diff_dist)

    df = pd.DataFrame(all_seg_dist).T
    df.columns = colnames

    return df


def get_segments_angle(u, segmanet_list1, segmanet_chainid1, labels1, segmanet_list2, segmanet_chainid2, labels2):

    # pre-calculate helix vector
    print("pre-calculate helix vector for segment1")
    all_seg_vec1 = []
    for i in range(len(segmanet_list1)):
        indCi = u.select_atoms('protein and segid PRO' + segmanet_chainid1 + ' and resid ' + segmanet_list1[i] + ' and name C')
        indOi = u.select_atoms('protein and segid PRO' + segmanet_chainid1 + ' and resid ' + segmanet_list1[i] + ' and name O')
        # print("debug1------->",len(indCi),len(indOi))
        # print("debug1------->",indCi)
        # print(segmanet_list1[i])
        # if len(indCi)==len(indOi):
        _seg_vec = []
        for ts in u.trajectory:
            c1i = indCi.positions
            c2i = indOi.positions
            _vi = cal_vec(c1i, c2i)
            _seg_vec.append(_vi)
        all_seg_vec1.append(_seg_vec)

    # pre-calculate helix vector
    print("pre-calculate helix vector for segment2")
    all_seg_vec2 = []
    for i in range(len(segmanet_list2)):
        indCi = u.select_atoms('protein and segid PRO' + segmanet_chainid2 + ' and resid ' + segmanet_list2[i] + ' and name C')
        indOi = u.select_atoms('protein and segid PRO' + segmanet_chainid2 + ' and resid ' + segmanet_list2[i] + ' and name O')
        # print("debug2------->",len(indCi),len(indOi))
        # if len(indCi)==len(indOi):
        _seg_vec = []
        for ts in u.trajectory:
            c1i = indCi.positions
            c2i = indOi.positions
            _vi = cal_vec(c1i, c2i)
            _seg_vec.append(_vi)
        all_seg_vec2.append(_seg_vec)


    print("pairwise angle")
    all_angle = []
    colnames = []
    for i in range(len(segmanet_list1)):
        _va = all_seg_vec1[i]
        _len_va = len(_va)

        for j in range(len(segmanet_list2)):
            _vb = all_seg_vec2[j]
            _len_vb = len(_vb)
            assert _len_va==_len_vb

            _ls_angle = []
            for k in range(_len_va):
                _ls_angle.append(get_angle(_va[k],_vb[k]))
                #print(i,j, get_angle(_va[k],_vb[k]))

            colnames.append("seg_angle_"+labels1[i] + "_" + labels2[j])
            all_angle.append(_ls_angle)

    df = pd.DataFrame(all_angle).T
    #df = df.T
    df.columns = colnames

    return df




def cal_vec(c1, c2):
    vecs = np.array([])
    # c1i = indCi.positions
    # c2i = indOi.positions
    vec_list = c2 - c1
    vec = np.sum(vec_list, axis=0)

    angle_list = [get_angle(vec, x) for x in vec_list]
    # calculate the angle between the average vec and each of the individual vec in the vec_list.
    angle_mu = np.mean(angle_list)
    angle_std = np.std(angle_list)
    # '''
    # compare each of the above calculated angle with angle_std * sigma_cutoff. if the angle value is greater than angle_std * sigma_cutoff,
    # exclude the outlier vector from the vec_list.
    # '''
    sigma_cutoff = 1.5
    nvec_list = [vec_list[i] for i in range(len(vec_list)) if
                 abs(angle_list[i] - angle_mu) < angle_std * sigma_cutoff]
    vec = np.sum(nvec_list,
                 axis=0)  # create an average vector for the helical segment after removing outliers
    # print "vec_list2 :\n", vec_list
    # print "vector2 :", vec
    vec = vec / np.linalg.norm(vec)
    return vec


def get_angle(u, v):
    """
    Description:
        angle calculation module. angles are given in radians. return a dictionary with frame number as the key.
    Arguments:
        vectors a and b
        frame id
    """
    u = np.squeeze(np.asarray(u))
    v = np.squeeze(np.asarray(v))

    cosa = np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v)  # -> cosine of the angle
    angle = np.arccos(np.clip(cosa, -1, 1))

    angle = angle / (2 * math.pi ) * 360
    return angle


def convert_angle(angle):
    """
    Description:
        module to convert radians to degrees.
    Arguments:
        vectors a and b
        frame id
    """
    angle_deg = math.degrees(float(angle))

    if angle_deg < 0:
        angle_deg += 360
        return angle_deg
    else:
        return angle_deg


def helix_orientation(u, segid, resid, sigma_cutoff=1.5):
    # ADAPTED AND MODIFIED FROM:
    #     anglebetweenhelices.py
    #     helix_angle_calc.py

    ind1 = u.select_atoms('protein and segid PRO'+segid+' and resid '+resid+' and name C')
    ind2 = u.select_atoms('protein and segid PRO'+segid+' and resid '+resid+' and name O')

    vecs = np.array([])
    n=0
    for ts in u.trajectory:
        c1 = ind1.positions
        c2 = ind2.positions
        vec_list = c2 - c1
        vec = np.sum(vec_list, axis=0)

        angle_list = [get_angle(vec, x) for x in vec_list]
        # calculate the angle between the average vec and each of the individual vec in the vec_list.
        angle_mu = np.mean(angle_list)
        angle_std = np.std(angle_list)
        # '''
        # compare each of the above calculated angle with angle_std * sigma_cutoff. if the angle value is greater than angle_std * sigma_cutoff,
        # exclude the outlier vector from the vec_list.
        # '''
        nvec_list = [vec_list[i] for i in range(len(vec_list)) if abs(angle_list[i] - angle_mu) < angle_std * sigma_cutoff]
        vec = np.sum(nvec_list, axis=0)  # create an average vector for the helical segment after removing outliers
        # print "vec_list2 :\n", vec_list
        # print "vector2 :", vec
        vec = vec / np.linalg.norm(vec)
        print(vec)
        # print(n)
    return vec





def get_rec_info2(bwtable):
    rec_segments_pia=[str(get_key(bwtable, "TM1.30"))+":"+str(get_key(bwtable, "TM1.33")),
                      str(get_key(bwtable, "TM2.63"))+":"+str(get_key(bwtable, "TM2.66")),
                      str(get_key(bwtable, "TM3.22"))+":"+str(get_key(bwtable, "TM3.25")),
                      str(get_key(bwtable, "TM4.59"))+":"+str(get_key(bwtable, "TM4.62")),
                      str(get_key(bwtable, "EL2.52")),
                      str(get_key(bwtable, "TM5.36"))+":"+str(get_key(bwtable, "TM5.39")),
                      str(get_key(bwtable, "TM6.57"))+":"+str(get_key(bwtable, "TM6.60")),
                      str(get_key(bwtable, "TM7.32"))+":"+str(get_key(bwtable, "TM7.35"))]
    rec_labels_pia=["1ee", "2ee",
                    "3ee", "4ee",'EL2',
                    "5ee", "6ee",
                    "7ee"]
    rec_chainid="R"
    return rec_segments_pia, rec_labels_pia, rec_chainid



if __name__ == '__main__':

    pdbin=argv[1]
    dcdin=argv[2]
    gpcrin=argv[3]
    gain = argv[4]

    u = mda.Universe(pdbin, dcdin)
    output_dir="output"
    check_output_dir(output_dir)

    bwtable = convert_Residue2TMs(gpcrin)

    # check gpcrin
    if gain == "go":
        ga_segments = gi_segments
    elif gain == "gi":
        ga_segments = go_segments
    else:
        print("ga_segments needs")

    rec_segments_pia, rec_labels_pia, rec_chainid, = get_rec_info(bwtable)

    # segments distance: rec vs ga
    _df=get_segments_distance(u, rec_segments_pia, rec_chainid, rec_labels_pia, ga_segments, ga_chainid, ga_labels)
    _df.to_hdf(output_dir + '/df_seg_dist_rec_ga.hdf', key='df', mode='w')

    # segments distance: rec vs rec (PIA)
    backbone_CA_in = True
    _df=get_segments_distance(u, rec_segments_pia, rec_chainid, rec_labels_pia, rec_segments_pia, rec_chainid, rec_labels_pia, backbone_CA_in)
    _df.to_hdf(output_dir + '/df_seg_dist_rec_rec.hdf', key='df', mode='w')

    # segments distance: rec vs rec (PIA); side-chain COM
    backbone_CA_in = False
    _df=get_segments_distance(u, rec_segments_pia, rec_chainid, rec_labels_pia, rec_segments_pia, rec_chainid, rec_labels_pia, backbone_CA_in)
    _df.to_hdf(output_dir + '/df_seg_dist_rec_rec_COM.hdf', key='df', mode='w')

    # segments angle: rec vs ga
    _df=get_segments_angle(u, rec_segments_pia, rec_chainid, rec_labels_pia, ga_segments, ga_chainid, ga_labels)
    _df.to_hdf(output_dir+'/df_seg_angle_rec_ga.hdf', key='df', mode='w')

    # segments angle: rec vs rec (PIA)
    _df=get_segments_angle(u, rec_segments_pia, rec_chainid, rec_labels_pia, rec_segments_pia, rec_chainid, rec_labels_pia)
    _df.to_hdf(output_dir+'/df_seg_angle_rec_rec.hdf', key='df', mode='w')



    # end matrix
    rec_segments_pia, rec_labels_pia, rec_chainid, = get_rec_info2(bwtable)

    # segments distance: rec vs rec (PIA)
    backbone_CA_in = True
    _df=get_segments_distance(u, rec_segments_pia, rec_chainid, rec_labels_pia, rec_segments_pia, rec_chainid, rec_labels_pia, backbone_CA_in)
    _df.to_hdf(output_dir + '/df_seg_dist_rec_rec_end.hdf', key='df', mode='w')

    # segments distance: rec vs rec (PIA); side-chain COM
    backbone_CA_in = False
    _df=get_segments_distance(u, rec_segments_pia, rec_chainid, rec_labels_pia, rec_segments_pia, rec_chainid, rec_labels_pia, backbone_CA_in)
    _df.to_hdf(output_dir + '/df_seg_dist_rec_rec_end_COM.hdf', key='df', mode='w')
