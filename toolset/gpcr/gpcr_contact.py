"""
Kuo Hao Lee
protein-protein association analyzer (PPAA)
the program is used to analyse protein-protein association
"""


import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
from MDAnalysis.analysis import contacts
from statistics import mean, stdev
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# %matplotlib inline
import time
import glob


def _contacts_within_cutoff(u, random_list, sel_a, sel_b, radius=5):
    timeseries = []
    # for ts in u.trajectory:
    for ts in u.trajectory[random_list]:
        # calculate distances between sel_a (selection a) and sel_b (selection b)
        dist = contacts.distance_array(sel_a.positions, sel_b.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts, int(n_contacts > 0)])
    return np.array(timeseries)


def _contacts_frequency(u, sel2,
                        sel_1_resids,
                        sel_1_resnames,
                        sel1_atoms,
                        contacts_cutoff=5):

    # bootstrapping 10 frames
    random_list = []
    np.random.seed(313)
    for a in range(0, 10):
        random_list.append(np.random.randint(0, len(u.trajectory)))
    # print(random_list)

    lig_neighbor_contacts = []
    for i in range(len(sel1_atoms)):
        ca = _contacts_within_cutoff(u, random_list, sel2, sel1_atoms[i], radius=contacts_cutoff)
        # ca.shape
        ca_df = pd.DataFrame(ca, columns=['Frame', '# Contacts', 'Contact'])

        lig_neighbor_contacts.append([sel_1_resids[i],
                                      sel_1_resnames[i],
                                      ca_df['# Contacts'].sum(),
                                      ca_df['Contact'].sum(),
                                      # len(u.trajectory),
                                      # format(ca_df['Contact'].sum()/len(u.trajectory)*100, '.2f')
                                      len(random_list),
                                      format(ca_df['Contact'].sum() / len(random_list) * 100, '.2f')
                                      ])

    ligand_contacts = pd.DataFrame(lig_neighbor_contacts,
                                   columns=['resids', 'resnames', '# atom contacts', '# residue contacts',
                                            'total frames', 'freq (%)'])
    return ligand_contacts


def get_df_contacts_frequency(u, chain1="R", chain2="A"):
    t1 = time.time()

    # only heavy atoms considered
    sel_chain1_atoms = u.select_atoms('protein and not type H and segid C' + chain1)
    sel_chain1_resids = list(sel_chain1_atoms.residues.resids)
    sel_chain1_resns = list(sel_chain1_atoms.residues.resnames)
    sel_chain2_atoms = u.select_atoms('protein and not type H and segid C' + chain2)

    # atom selection array for ligand neighbor resids
    lig_neighbor_selection = []
    for i in range(len(sel_chain1_resids)):
        _sel = str(
            "protein and not type H and segid C" + chain1 + " and resid " + str(list(sel_chain1_atoms.residues.resids)[i]))
        lig_neighbor_selection.append(u.select_atoms(_sel))

    df = _contacts_frequency(u,
                             sel_chain2_atoms,
                             sel_chain1_resids,
                             sel_chain1_resns,
                             lig_neighbor_selection,
                             contacts_cutoff=5)

    t2 = time.time()
    print("TIMER: Function completed time is: %.5f" % (t2 - t1))
    return df



def get_union(in_df1, in_df2):
    # given two residue ids and returen the union

    # find residue with contact frequency bigger than 0
    _df_1 = in_df1[in_df1['freq (%)'] != '0.00'].reset_index()
    _df_2 = in_df2[in_df2['freq (%)'] != '0.00'].reset_index()

    # find the union of set between Recptor-proteinG-complex
    union_set = list(set(_df_1['resids']) | set(_df_2['resids']))
    union_set.sort()
    # print(union_set)

    return union_set


def get_union_df(in_df, in_resids):
    # given dataset and a range of resids
    # retrun the dataset with selection of resids

    df = pd.concat([in_df[in_df['resids'] == in_resids[i]] for i in range(len(in_resids))]).reset_index()

    return df





def get_distance_matrix(PDB, DCDs, common_rec, common_ga):
    t1 = time.time()

    sel_resid_REC = "resid " + " or resid ".join(str(x) for x in common_rec)
    sel_resid_Ga = "resid " + " or resid ".join(str(x) for x in common_ga)

    u = mda.Universe(PDB, DCDs)

    selA = u.select_atoms('protein and segid CR and (' + sel_resid_REC + ') and name CA')
    selB = u.select_atoms('protein and segid CA and (' + sel_resid_Ga + ') and name CA')

    # atom selection array for selA
    ls_selA = []
    for i in range(len(list(selA.residues.resids))):
        _sel = str("protein and segid CR and name CA and resid " + str(list(selA.residues.resids)[i]))
        ls_selA.append(u.select_atoms(_sel))
    # ls_selA

    ## atom selection array for selB
    ls_selB = []
    for i in range(len(list(selB.residues.resids))):
        _sel = str("protein and segid CA and name CA and resid " + str(list(selB.residues.resids)[i]))
        ls_selB.append(u.select_atoms(_sel))
    # ls_selB

    # create matrix
    matrix = np.zeros((len(ls_selA), len(ls_selB)))

    # define two sets
    _selA = ls_selA
    _selB = ls_selB

    # TODO: move this part into bootstrapping
    # define a list of random
    random_list = []
    np.random.seed(313)
    for a in range(0, 100):
        random_list.append(np.random.randint(0, len(u.trajectory)))
    # print(random_list)

    # calculate distance matrix
    for i in range(0, len(_selA)):
        for j in range(0, len(_selB)):
            timeseries = []
            # for ts in u.trajectory:
            for ts in u.trajectory[random_list]:
                # dist = contacts.distance_array(ls_selA[i].positions, ls_selB[j].positions).tolist()[0][0]
                dist = contacts.distance_array(_selA[i].positions, _selB[j].positions).tolist()[0][0]
                timeseries.append(dist)
            matrix[i][j] = mean(timeseries)

    # matrix

    t2 = time.time()
    print("TIMER: Function completed time is: %.5f" % (t2 - t1))
    return (matrix)


def get_dist_freq_matrix(u, common_resids1, chain1, common_resids2,chain2):
    t1 = time.time()

    # atom selection array for selA
    ls_selA = []
    for i in range(len(common_resids1)):
        _sel = str('protein and segid C'+chain1+' and not type H and resid ' + str(common_resids1[i]))
        ls_selA.append(u.select_atoms(_sel))

    # atom selection array for selB
    ls_selB = []
    for i in range(len(common_resids2)):
        _sel = str('protein and segid C'+chain2+' and not type H and resid ' + str(common_resids2[i]))
        ls_selB.append(u.select_atoms(_sel))

    # create matrix
    matrix = np.zeros((len(ls_selA), len(ls_selB)))
    matrix_sd = np.zeros((len(ls_selA), len(ls_selB)))
    matrix_freq = np.zeros((len(ls_selA), len(ls_selB)))

    # define two sets
    _selA = ls_selA
    _selB = ls_selB

    # calculate distance matrix
    for i in range(0, len(_selA)):
        for j in range(0, len(_selB)):
            timeseries = []
            counts = 0
            for ts in u.trajectory:
                dist = np.amin(contacts.distance_array(_selA[i].positions, _selB[j].positions))
                timeseries.append(dist)
                if dist <= 5:
                    counts = counts + 1
            # TODO: add more standard here to skip far-away pairs
            matrix[i][j] = mean(timeseries)
            matrix_sd[i][j] = stdev(timeseries)
            matrix_freq[i][j] = counts / len(u.trajectory)

    t2 = time.time()
    print("TIMER: Function completed time is: %.5f" % (t2 - t1))
    return matrix, matrix_sd, matrix_freq

def get_matrix_diff(matrix_freq_1, matrix_freq_2):
    matrix_i = matrix_freq_1.shape[0]
    matrix_j = matrix_freq_1.shape[1]
    matrix_diff = np.zeros((matrix_i, matrix_j))

    # calculate distance matrix
    for i in range(matrix_i):
        for j in range(matrix_j):
            matrix_diff[i][j] = matrix_freq_1[i][j] - matrix_freq_2[i][j]
    return matrix_diff


def plot_matrix(u, mymatrix, common_rec, common_ga):
    sel_resid_REC = "resid " + " or resid ".join(str(x) for x in common_rec)
    sel_resid_Ga = "resid " + " or resid ".join(str(x) for x in common_ga)

    selA = u.select_atoms('protein and segid CR and (' + sel_resid_REC + ') and name CA')
    selB = u.select_atoms('protein and segid CA and (' + sel_resid_Ga + ') and name CA')

    dfA = {"resid": list(selA.residues.resids), "resname": list(selA.residues.resnames)}
    dfA = pd.DataFrame(dfA)
    dfA["label"] = dfA["resid"].astype('str') + "_" + dfA["resname"].astype('str')

    dfB = {"resid": list(selB.residues.resids), "resname": list(selB.residues.resnames)}
    dfB = pd.DataFrame(dfB)
    dfB["label"] = dfB["resid"].astype('str') + "_" + dfB["resname"].astype('str')

    mynamesA = dfA["label"]
    mynamesB = dfB["label"]
    # mynamesA=dfA["resid"]
    # mynamesB=dfB["resid"]

    # cbticks=np.arange(-npa,npa+1,npa)
    # cbticks = np.arange(3, 8, 0.5)
    cbticks = np.arange(0, 1, 0.1)

    plt.figure(figsize=(18, 18))
    sz = 20
    plt.imshow(mymatrix, interpolation='nearest', cmap=plt.cm.hot)
    cbar = plt.colorbar(orientation='horizontal', shrink=0.2, pad=0.1, aspect=4, ticks=cbticks)
    cbar.ax.tick_params(labelsize=sz)
    l = max(np.hstack(mymatrix).min(), np.hstack(mymatrix).max(), key=abs)

    plt.clim(cbticks[0], cbticks[-1])
    ax = plt.gca()
    plt.xticks(np.arange(0, len(mynamesB), 1), size=sz, rotation=90)
    ax.set_xticklabels(mynamesB)
    plt.yticks(np.arange(0, len(mynamesA), 1), size=sz)
    ax.set_yticklabels(mynamesA)
    ax.yaxis.labelpad = 100
    ax.grid(True)
    plt.tight_layout()

    plt.xlabel('Ga', fontsize=30)
    plt.ylabel('REC', fontsize=30)
    plt.title('CA(Ga)-CA(REC) Distance', fontsize=30)
    return




def plot_matrix_diff(u, mymatrix, common_rec, common_ga):
    sel_resid_REC = "resid " + " or resid ".join(str(x) for x in common_rec)
    sel_resid_Ga = "resid " + " or resid ".join(str(x) for x in common_ga)

    selA = u.select_atoms('protein and segid CR and (' + sel_resid_REC + ') and name CA')
    selB = u.select_atoms('protein and segid CA and (' + sel_resid_Ga + ') and name CA')

    dfA = {"resid": list(selA.residues.resids), "resname": list(selA.residues.resnames)}
    dfA = pd.DataFrame(dfA)
    dfA["label"] = dfA["resid"].astype('str') + "_" + dfA["resname"].astype('str')

    dfB = {"resid": list(selB.residues.resids), "resname": list(selB.residues.resnames)}
    dfB = pd.DataFrame(dfB)
    dfB["label"] = dfB["resid"].astype('str') + "_" + dfB["resname"].astype('str')

    mynamesA = dfA["label"]
    mynamesB = dfB["label"]
    # mynamesA=dfA["resid"]
    # mynamesB=dfB["resid"]

    # cbticks=np.arange(-npa,npa+1,npa)
    cbticks = np.arange(-8, 8, 0.5)

    plt.figure(figsize=(22, 22))
    sz = 20
    plt.imshow(mymatrix, interpolation='nearest', cmap=plt.cm.hot)
    cbar = plt.colorbar(orientation='horizontal', shrink=0.2, pad=0.1, aspect=4, ticks=cbticks)
    cbar.ax.tick_params(labelsize=sz)
    l = max(np.hstack(mymatrix).min(), np.hstack(mymatrix).max(), key=abs)

    plt.clim(cbticks[0], cbticks[-1])
    ax = plt.gca()
    plt.xticks(np.arange(0, len(mynamesB), 1), size=sz, rotation=90)
    ax.set_xticklabels(mynamesB)
    plt.yticks(np.arange(0, len(mynamesA), 1), size=sz)
    ax.set_yticklabels(mynamesA)
    ax.yaxis.labelpad = 100
    ax.grid(True)
    plt.tight_layout()

    plt.xlabel('Ga', fontsize=30)
    plt.ylabel('REC', fontsize=30)
    plt.title('CA(Ga)-CA(REC) Distance Difference', fontsize=30)
    return


def get_labels(u, common_rec, common_ga):
    sel_resid_REC = "resid " + " or resid ".join(str(x) for x in common_rec)
    sel_resid_Ga = "resid " + " or resid ".join(str(x) for x in common_ga)

    selA = u.select_atoms('protein and segid CR and (' + sel_resid_REC + ') and name CA')
    selB = u.select_atoms('protein and segid CA and (' + sel_resid_Ga + ') and name CA')

    dfA = {"resid": list(selA.residues.resids), "resname": list(selA.residues.resnames)}
    dfA = pd.DataFrame(dfA)
    dfA["label"] = dfA["resid"].astype('str') + "_" + dfA["resname"].astype('str')

    dfB = {"resid": list(selB.residues.resids), "resname": list(selB.residues.resnames)}
    dfB = pd.DataFrame(dfB)
    dfB["label"] = dfB["resid"].astype('str') + "_" + dfB["resname"].astype('str')

    xlabels = dfA["label"]
    ylabels = dfB["label"]
    y_axis_label = "Rec"
    x_axis_label = "Ga"
    return xlabels, ylabels, x_axis_label, y_axis_label


def get_labels2(u, resids, chainid):
    sel_resid_REC = "resid " + " or resid ".join(str(x) for x in resids)

    selA = u.select_atoms('protein and segid C'+chainid+' and (' + sel_resid_REC + ') and name CA')

    dfA = {"resid": list(selA.residues.resids), "resname": list(selA.residues.resnames)}
    dfA = pd.DataFrame(dfA)
    dfA["label"] = dfA["resid"].astype('str') + "_" + dfA["resname"].astype('str')

    xlabels = dfA["label"]
    return xlabels