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
import numpy as np
import math
import pandas as pd
import MDAnalysis as mda
from IPython.display import display
from itertools import combinations
import re
import matplotlib.pyplot as plt



# predefine counter-clockwise extracellular turns order (from the top view)
dict_cc = {
    'x2ee': ['x3ee','x4ee','x5ee','x6ee','x7ee'],
    'x3ee': ['x4ee','x5ee','x6ee','x7ee','x2ee'],
    'x4ee': ['x5ee','x6ee','x7ee','x2ee','x3ee'],
    'x5ee': ['x6ee','x7ee','x2ee','x3ee','x4ee'],
    'x6ee': ['x7ee','x2ee','x3ee','x4ee','x5ee'],
    'x7ee': ['x2ee','x3ee','x4ee','x5ee','x6ee']
    }


# list of extracellular turns order
list_xee = ['x2ee','x3ee','x4ee','x5ee','x6ee','x7ee']


def get_ee_definition(bwtable):
    # ex to get residue ID
    print("1e", "1.30 - 1.33", get_key(bwtable, "TM1.30"), get_key(bwtable, "TM1.33"))
    print("2e", "2.63 - 2.66", get_key(bwtable, "TM2.63"), get_key(bwtable, "TM2.66"))
    print("3e", "3.22 - 3.25", get_key(bwtable, "TM3.22"), get_key(bwtable, "TM3.25"))
    print("4e", "4.59 - 4.62", get_key(bwtable, "TM4.59"), get_key(bwtable, "TM4.62"))
    print("5e", "5.36 - 5.39", get_key(bwtable, "TM5.36"), get_key(bwtable, "TM5.39"))
    print("6e", "6.57 - 6.60", get_key(bwtable, "TM6.57"), get_key(bwtable, "TM6.60"))
    print("7e", "7.32 - 7.35", get_key(bwtable, "TM7.32"), get_key(bwtable, "TM7.35"))
    print("EL2", "EL2.52", get_key(bwtable, "EL2.52"))


def prep_resid_for_xee(bwtable):
    '''
    prepare residue ids for given extracellular ends
    :param bwtable:
    :return:
    '''
    # residue ids for the extracellular turns
    x2ee = str(get_key(bwtable, "TM2.63"))+':'+str(get_key(bwtable, "TM2.66"))
    x3ee = str(get_key(bwtable, "TM3.22"))+':'+str(get_key(bwtable, "TM3.25"))
    x4ee = str(get_key(bwtable, "TM4.59"))+':'+str(get_key(bwtable, "TM4.62"))
    x5ee = str(get_key(bwtable, "TM5.36"))+':'+str(get_key(bwtable, "TM5.39"))
    x6ee = str(get_key(bwtable, "TM6.57"))+':'+str(get_key(bwtable, "TM6.60"))
    x7ee = str(get_key(bwtable, "TM7.32"))+':'+str(get_key(bwtable, "TM7.35"))
    return [x2ee, x3ee, x4ee, x5ee, x6ee, x7ee]


def get_area(a,b,c):
    '''given 3 points & return area'''
    v1 = b-a
    v2 = a-c
    return round(np.linalg.norm(np.cross(v1, v2))/2,2)


def get_COM_by_xXee(u, xxee):
    '''provide u & xxee, return COM for xxee
    xxee is any extracecullar turns
    '''
    return u.select_atoms('protein and segid PROR and resid '+xxee+' and name CA').center_of_mass()



def get_area_by_sel(xee_a,xee_b,xee_c, u):
    '''given 3 points & return area'''
    a = get_COM_by_xXee(u, xee_a)
    b = get_COM_by_xXee(u, xee_b)
    c = get_COM_by_xXee(u, xee_c)
    v1 = b-a
    v2 = a-c
    return round(np.linalg.norm(np.cross(v1, v2))/2,2)



def get_pdb_classes(subtype, xdict):
    # get keys & values from dictionary
    val_list = list(xdict.values())
    key_list = list(xdict.keys())

    # active pdbs
    act_pdbs = []
    for i in range(len(val_list)):
        val = val_list[i]
        key = key_list[i]
        if val.replace(subtype + '.', "") == "a":
            act_pdbs.append(key)

    # inactive pdbs
    inact_pdbs = sorted(key_list)
    for i_act_pdb in act_pdbs:
        if i_act_pdb in inact_pdbs:
            inact_pdbs.remove(i_act_pdb)
    inact_pdbs = sorted(inact_pdbs)

    # all pdbs
    all_pdbs = sorted(key_list)
    return act_pdbs, inact_pdbs, all_pdbs


def get_common_trend(df, act_pdbs, inact_pdbs, subtotal_only=False, do_return=False):
    dfin = df.reset_index(drop=False)
    if subtotal_only:
        dfin = dfin[["subtotal" in i for i in dfin['labels']]].set_index('labels')
    else:
        dfin = dfin[["subtotal" not in i for i in dfin['labels']]].set_index('labels')
    num_inact = len(inact_pdbs)
    num_act = len(act_pdbs)
    num_total = num_inact*num_act

    print("number of inactive: " + str(num_inact)+" ==> "+str(inact_pdbs))
    print("number of active: " +   str(num_act)+"   ==> "+str(act_pdbs))

    # inactive > active
    for i_act in act_pdbs:
        _dfout = pd.DataFrame()
        for i_inact in inact_pdbs:
            i_inact = i_inact.lower()
            _dfout[i_inact] = dfin[i_inact] > dfin[i_act]
            _dfout[i_inact] = _dfout[i_inact].astype(int)
        _dfout['sum'] = _dfout.sum(axis=1)

        dfout = pd.DataFrame()
        dfout = pd.concat([dfout, _dfout[_dfout['sum'] == 0]])
        print("common [inactive > active]: ", len(_dfout[_dfout['sum'] == num_inact]))
        dfout = pd.concat([dfout, _dfout[_dfout['sum'] == num_inact]])
        dfout = dfout[dfout['sum']==num_total]
        dfout_biger_inact = dfout

    # active > inactive
    for i_act in act_pdbs:
        _dfout = pd.DataFrame()
        for i_inact in inact_pdbs:
            i_inact = i_inact.lower()
            _dfout[i_inact] = dfin[i_inact] < dfin[i_act]
            _dfout[i_inact] = _dfout[i_inact].astype(int)

        _dfout['sum'] = _dfout.sum(axis=1)

        dfout = pd.DataFrame()
        dfout = pd.concat([dfout, _dfout[_dfout['sum'] == 0]])
        print("common [inactive < active]: ", len(_dfout[_dfout['sum'] == num_inact]))
        dfout = pd.concat([dfout, _dfout[_dfout['sum'] == num_inact]])
        dfout = dfout[dfout['sum']==num_total]
        dfout_biger_act = dfout

    if do_return:
        return dfout_biger_inact, dfout_biger_act


def reload_area(subtype, all_pdbs):
    '''reload the pre-analyzed data'''

    df = pd.DataFrame()
    for ipdb in all_pdbs:
        _df = pd.read_hdf('../xtal_rawdata/'+subtype+'/'+ipdb+'_R_output/df_area.hdf')
        df = pd.concat([df,_df],axis=1)
    df.columns = [re.sub("_ligand","", i) for i in list(df.columns)]
    df.columns = [re.sub("_R","", i) for i in list(df.columns)]
    return df

def get_cc_area(df, dict_cc):

    df_ref = df.reset_index(drop=False)  # all 20 combination is used as reference
    dfcc = pd.DataFrame()  # dataframe for counter-clockwise extracellular turns
    verbose = False
    for i_cc in list(dict_cc.keys()):
        #i_cc = 'x2ee' # find the top from pre-define counter-clockwise extracellular turns order
        ls_cc = dict_cc[i_cc] # the list of counter-clockwise extracellular turns order
        if verbose: print(i_cc, ls_cc)
        _df_icc = pd.DataFrame()
        ls_ori_str = []
        for i in range(len(ls_cc)-1):
            ls_ls = list([i_cc, ls_cc[i],ls_cc[i+1]]) # list of order for 3 points
            ori_str = str(ls_ls).replace("[","").replace("]","").replace("'","") # convert 3 points to order
            if verbose: print(ori_str)
            key_str = str(sorted(ls_ls)).replace("[","").replace("]","").replace("'","") # convert 3 points to order
            _df_key_str = df_ref[df_ref['labels'] == key_str]
            _df_icc = pd.concat([_df_icc,_df_key_str],axis=0).reset_index(drop=True)
            ls_ori_str.append(ori_str)
        _df_icc['labels'] = ls_ori_str
        _df_icc.loc[len(_df_icc),:] = _df_icc.sum(axis=0)
        _df_icc['labels'][len(_df_icc)-1] = list(_df_icc['labels'])[-1].split(',')[0]+" subtotal"
        if verbose: display(_df_icc)
        dfcc = pd.concat([dfcc,_df_icc],axis=0).reset_index(drop=True)
    dfcc = dfcc.set_index('labels')
    return dfcc


def get_4points_area(dfall):
    '''
    dfall is the look up table
    consider the combination and take area from the look up table
    then return the final df
    '''
    df_ref = dfall.reset_index(drop=False)  # all 20 combination is used as reference

    dfout = pd.DataFrame()
    list_combinations = list()
    list_combinations += list(combinations(list_xee, 4))
    newlabels = []
    for i_comb in list_combinations:
        # print(i_comb)
        _list = list(i_comb)
        sublist_comb = list()
        sublist_comb += list(combinations(_list, 3))
        _sub_dfout = pd.DataFrame()
        for i_sub_comb in sublist_comb:
            _newlabel = str(i_comb) + "; " + str(i_sub_comb)
            newlabels.append((str(i_comb) + "____" + str(i_sub_comb)).replace('(', '').replace(')', '').replace("'", ''))
            key_str = str(i_sub_comb).replace('(', '').replace(')', '').replace("'", '')
            _df_key_str = df_ref[df_ref['labels'] == key_str]
            _df_key_str['labels'][0] = _newlabel
            # display(_df_key_str)
            _sub_dfout = pd.concat([_sub_dfout, _df_key_str], axis=0).reset_index(drop=True)
        #_df_icc['labels'] = ls_ori_str
        _sub_dfout.loc[len(_sub_dfout),:] = _sub_dfout.sum(axis=0)
        _sub_dfout['labels'][len(_sub_dfout)-1] = str(i_comb)+" subtotal"
        newlabels.append(str(i_comb)+" subtotal")

        dfout = pd.concat([dfout, _sub_dfout], axis=0).reset_index(drop=True)
    dfout['labels'] = newlabels
    dfout = dfout.set_index('labels')

    return dfout

def get_4p_2tri(dfall):
    '''
    dfall is the look up table
    consider the combination and take area from the look up table
    then return the final df
    '''
    list_combinations = list()
    list_combinations += list(combinations(list_xee, 4))
    n_combinations = len(list_combinations)


    dfref = dfall.copy()
    dfref = dfref.reset_index(drop=False)

    _names = []
    _labels_sorted = []
    df_4p_2tri = pd.DataFrame()

    for i in range(n_combinations):
        _list_sel = list_combinations[i]
        list_idx = [] # prepare index for convience coverting name in list_xee
        [list_idx.append(list_xee.index(i)) for i in _list_sel]
        # print(_list_sel)    # name in list_xee
        # print(list_idx) # index in list_xee

        _df_2tri_sum = pd.DataFrame()
        for ii in list_idx:
            xx = list_idx.copy()
            xx.remove(ii)
            for iidx in range(len(xx)-1):
                #print(_list_sel, "____", list_xee[ii], list_xee[xx[iidx]], list_xee[xx[iidx+1]])
                _name = list_xee[ii], list_xee[xx[iidx]], list_xee[xx[iidx+1]]
                _name = str(_name).replace("(","").replace(")","").replace("'","")
                _name_4p = _list_sel
                _name_4p = str(_name_4p).replace("(","").replace(")","").replace("'","")

                _label_sorted = sorted([list_xee[ii], list_xee[xx[iidx]], list_xee[xx[iidx+1]]])
                _label_sorted = str(_label_sorted).replace("[","").replace("]","").replace("'","")
    #            print(_name_4p+"___"+_name)
                _names.append(_name_4p+"___"+_name)
                _labels_sorted.append(_label_sorted)
                _df = dfref[dfref['labels'] == _label_sorted]
                _df_2tri_sum = pd.concat([_df_2tri_sum, _df], axis=0).reset_index(drop=True)
        # sum up subset
        _df_2tri_sum.loc[len(_df_2tri_sum),:] = _df_2tri_sum.sum(axis=0)
        _df_2tri_sum['labels'][len(_df_2tri_sum)-1] = str(_list_sel)+" subtotal"
        _names.append(str(_list_sel)+" subtotal")
        df_4p_2tri = pd.concat([df_4p_2tri, _df_2tri_sum], axis=0).reset_index(drop=True)

    df_4p_2tri['labels'] = _names
    df_4p_2tri = df_4p_2tri.set_index('labels')
    return df_4p_2tri

if __name__ == '__main__':

    pdbin=argv[1]
    gpcrin=argv[2]

    u = mda.Universe(pdbin)
    output_dir="output"
    check_output_dir(output_dir)

    bwtable = convert_Residue2TMs(gpcrin)

    #print(list_xee)
    list_xee_resid = prep_resid_for_xee(bwtable)


    list_combinations = list()
    list_combinations += list(combinations(list_xee, 3))
    list_resids = prep_resid_for_xee(bwtable)

    #print(list_combinations)
    labels = []
    areas = []
    for i in range(len(list_combinations)):
        labels.append(str(list_combinations[i]).replace(")","").replace("(","").replace("'",""))

        areas.append(get_area_by_sel(list_xee_resid[list_xee.index(list_combinations[i][0])],
                                     list_xee_resid[list_xee.index(list_combinations[i][1])],
                                     list_xee_resid[list_xee.index(list_combinations[i][2])],u))

    #print(labels)
    #print(areas)
    df = pd.DataFrame()
    df['labels'] = labels
    df[pdbin.replace(".pdb","")] = areas
    df = df.set_index('labels', drop=True)
    #display(df)
    df.to_hdf(output_dir + '/df_area.hdf', key='df', mode='w')




# x2ee, x3ee, x4ee, x5ee, x6ee, x7ee
# list(list_com[0])


# dict_coor = {
#     'x2ee': [8, 6],
#     'x3ee': [6, 8],
#     'x4ee': [2, 6],
#     'x5ee': [1, 4],
#     'x6ee': [4, 2],
#     'x7ee': [6, 3]
#     }



def get_xx(_df):
    xx0 = list(_df.index)
    if "subtotal" in xx0[0]:
        do_4points = True
        do_3points = False
    else:
        do_4points = False
        do_3points = True
    xx=[]
    for i in range(len(xx0)):
        if do_4points:
            xx.append(xx0[i].replace(") subtotal","").replace("(","").replace("'","").replace('"','').split(', '))
        elif do_3points:
            xx.append(xx0[i].split(', '))
    return xx

def plot_area(_df1,_df2,intitle,dict_coor):

    if len(_df1) > 0:
        xx1 = get_xx(_df1)
        print("list1", xx1)
    if len(_df2) > 0:
        xx2 = get_xx(_df2)
        print("list2", xx2)

    circle2 = plt.Circle(dict_coor['x2ee'], 2, color='lightgrey', clip_on=False)
    circle3 = plt.Circle(dict_coor['x3ee'], 2, color='lightgrey', clip_on=False)
    circle4 = plt.Circle(dict_coor['x4ee'], 2, color='lightgrey', clip_on=False)
    circle5 = plt.Circle(dict_coor['x5ee'], 2, color='lightgrey', clip_on=False)
    circle6 = plt.Circle(dict_coor['x6ee'], 2, color='lightgrey', clip_on=False)
    circle7 = plt.Circle(dict_coor['x7ee'], 2, color='lightgrey', clip_on=False)

    figwidth=6
    figheight=6
    fig, ax = plt.subplots(figsize=(figwidth,figheight))

    plt.text(dict_coor['x2ee'][0]-1, dict_coor['x2ee'][1]+2, "TM2", fontsize=16)
    plt.text(dict_coor['x3ee'][0]+2, dict_coor['x3ee'][1], "TM3", fontsize=16)
    plt.text(dict_coor['x4ee'][0]-1, dict_coor['x4ee'][1]-3, "TM4", fontsize=16)
    plt.text(dict_coor['x5ee'][0]-1, dict_coor['x5ee'][1]-3, "TM5", fontsize=16)
    plt.text(dict_coor['x6ee'][0]-4, dict_coor['x6ee'][1], "TM6", fontsize=16)
    plt.text(dict_coor['x7ee'][0]-1, dict_coor['x7ee'][1]+2, "TM7", fontsize=16)

    coor_shift = 0.1
    if len(_df1) > 0:
        for i in range(len(xx1)):
            list_combinations = list()
            list_combinations += list(combinations(xx1[i], 2))
            n_combinations = len(list_combinations)
            for j in range(len(list_combinations)):
                _coorA = dict_coor[list_combinations[j][0]]
                _coorB = dict_coor[list_combinations[j][1]]
                coorA = [_coorA[0]-coor_shift,_coorB[0]-coor_shift]
                coorB = [_coorA[1]-coor_shift,_coorB[1]-coor_shift]
                ax.plot(coorA, coorB, 'b', linestyle="-")

    if len(_df2) > 0:
        for i in range(len(xx2)):
            list_combinations = list()
            list_combinations += list(combinations(xx2[i], 2))
            n_combinations = len(list_combinations)
            for j in range(len(list_combinations)):
                _coorA = dict_coor[list_combinations[j][0]]
                _coorB = dict_coor[list_combinations[j][1]]
                coorA = [_coorA[0]+coor_shift,_coorB[0]+coor_shift]
                coorB = [_coorA[1]+coor_shift,_coorB[1]+coor_shift]
                ax.plot(coorA, coorB, 'r', linestyle="-")


    ax.add_patch(circle2)
    ax.add_patch(circle3)
    ax.add_patch(circle4)
    ax.add_patch(circle5)
    ax.add_patch(circle6)
    ax.add_patch(circle7)
    #ax.set_xlim(0, 10)
    #ax.set_ylim(0, 10)
    ax.set_title(intitle,fontsize=20)
    ax.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
    plt.box(False)
