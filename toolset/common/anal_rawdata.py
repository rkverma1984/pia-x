import warnings
warnings.filterwarnings('ignore')
import glob
import pickle
# import seaborn as sns
import pandas as pd
import numpy as np
from natsort import natsorted
import re
import statistics
import multiprocessing
from multiprocessing import Process, Pool
ncpus = multiprocessing.cpu_count()

import glob
from natsort import natsorted
import pandas as pd


def get_bootstrap_data(receptor, ligand, bsdatadir, ls_includes_hdf, ls_includes_p):
    ##
    _filein = glob.glob(bsdatadir + "/" + receptor + "_" + ligand + "/bootstrap_[0-9]*")
    _filein = natsorted(_filein)

    df_out = pd.DataFrame()
    for i in range(len(_filein)):
        bsidx = _filein[i].split("/")[-1].replace("bootstrap_", "bs")
        df = pd.DataFrame()

        for ihdf in ls_includes_hdf:
            df2 = pd.read_hdf(_filein[i] + '/' + ihdf + '.hdf', 'df')
            df2 = df2.reset_index()
            del df2['index']
            df = pd.concat([df, df2], axis=1, ignore_index=False)

        for inc in ls_includes_p:
            pin = pickle.load(open(_filein[i] + '/' + inc + '.p', "rb"))
            df[inc] = pin

        df.insert(0, "Bootstrap", bsidx)
        df.insert(0, "Jobname", receptor + "_" + ligand)
        df.insert(0, "Ligand", ligand)
        df.insert(0, "Receptor", receptor)
        df_out = pd.concat([df_out, df])
    return df_out


def get_bootstrap_data_append(df, receptors, ligands, bsdatadir, ls_includes_hdf, ls_includes_p):
    for receptor in receptors:
        for ligand in ligands:
            _df = get_bootstrap_data(receptor, ligand, bsdatadir, ls_includes_hdf, ls_includes_p)
            df = pd.concat([df, _df], ignore_index=True)
    return df


# def read_data(receptor, ligand, time_limit=0, Verbose=False, exclude=[], include=[]):
#     """ Returns a dataframe containing all data
#
#     :param receptor: (str) string representing receptor name ex. 'd2gi'
#     :param ligand: (str) string representing ligand name ex. 'pd' or 'prm'
#     :param time_limit: (int) only include data where time is >= time_limit
#     :param Verbose: (bool) prints path of files to be read in if True
#     :param exclude: (sequence of str) jobnames to exclude
#     :param include: (sequence of str) list of pickle files to include
#     :return: dataframe containing all data from specified hdf files
#
#     """
#     _filein = glob.glob("/home/khlee/work/desmond/output/" + receptor + "/" + ligand + "/%s_*/" % receptor)
#     _filein = natsorted(_filein, key=lambda y: y.lower())
#     if Verbose:
#         print(_filein)
#     df_out = pd.DataFrame()
#     for i in range(len(_filein)):
#         if _filein[i].split("/")[-2] in exclude:
#             continue
#         try:
#             t = np.loadtxt(_filein[i] + '/ene/time.dat')
#             frame = list(range(1, len(t) + 1))
#             obs = pickle.load(open(_filein[i] + '/output/gpcr_interface.p', "rb"))
#             if (len(t)!=len(obs)): print("obs:"+str(len(obs)))
#             obs_a = pickle.load(open(_filein[i] + '/output/gpcr_interface_RA.p', "rb"))
#             if (len(t)!=len(obs_a)): print("obs_a:" + str(len(obs_a)))
#             obs_b = pickle.load(open(_filein[i] + '/output/gpcr_interface_RB.p', "rb"))
#             if (len(t)!=len(obs_b)): print("obs_b:" + str(len(obs_b)))
#
#         except:
#             print("Job not loaded because a file was not found:", _filein[i].split("/")[-2])
#             print("\n")
#             continue
#
#         try:
#             assert (len(t) == len(obs)) & (len(obs) == len(obs_a)) & (len(obs) == len(obs_b))
#             rec = [receptor for a in range(len(t))]
#             lig = [ligand for a in range(len(t))]
#             jobs = [_filein[i].split("/")[-2] for a in range(len(t))]
#             df = pd.DataFrame({"Receptor": rec,
#                                "Ligand": lig,
#                                "Jobname": jobs,
#                                "Time": t, "Frame": frame,
#                                "Interface": obs,
#                                "Interface_A": obs_a,
#                                "Interface_B": obs_b})
#             if len(include) > 0:
#                 df2 = pd.DataFrame()
#                 _files = glob.glob(_filein[i]+"/output/*.p")
#                 for inc in include:
#                     for _file in _files:
#                         if inc in str(_file):
#                             df2[inc] = pickle.load(open(_file, 'rb'))
#                             if (len(t)!=len(df2[inc])): print(str(inc)+": "+ str(len(df2[inc])))
#                 assert (len(df) == len(df2))
#                 df = pd.concat([df,df2], axis=1)
#         except:
#             print("Job not loaded because number of datapoints mismatch: ",
#                   _filein[i].split("/")[-2])
#             print("Time:", len(t), "Interface:", len(obs), "Interface_A:", len(obs_a), "Interface_B:", len(obs_b))
#             print("\n")
#             continue
#
#         try:
#             df1 = pd.read_hdf(_filein[i] + '/output/df_seg_dist_rec_ga.hdf', 'df')
#             if (df1.isnull().values.any()): print("check: df_seg_dist_rec_ga.hdf")
#             #print(len(df1))
#             df2 = pd.read_hdf(_filein[i] + '/output/df_seg_dist_rec_rec.hdf', 'df')
#             if (df2.isnull().values.any()): print("check: df_seg_dist_rec_rec.hdf")
#             #print(len(df2))
#             df3 = pd.read_hdf(_filein[i] + '/output/df_seg_angle_rec_ga.hdf', 'df')
#             if (df3.isnull().values.any()): print("check: df_seg_angle_rec_ga.hdf")
#             #print(len(df3))
#             df4 = pd.read_hdf(_filein[i] + '/output/df_seg_angle_rec_rec.hdf', 'df')
#             if (df4.isnull().values.any()): print("check: df_seg_angle_rec_rec.hdf")
#             #print(len(df4))
#             df5 = pd.read_hdf(_filein[i] + '/output/df_hhdist_RA.hdf', 'df')
#             if (df5.isnull().values.any()): print("check: df_hhdist_RA.hdf")
#             #print(len(df5))
#             df6 = pd.read_hdf(_filein[i] + '/output/df_hhdist_RB.hdf', 'df')
#             if (df6.isnull().values.any()): print("check: df_hhdist_RB.hdf")
#             #print(len(df6))
#             #             df5 = pd.read_hdf(_filein[i]+'/output/rmsdbb_alin_OBS_backbone.hdf','df')
#             #             df5.columns =["alin_obs_"+df5.columns[i] for i in range(len(df5.columns))]
#             #             df6 = pd.read_hdf(_filein[i]+'/output/rmsdbb_alin_prot_A_backbone.hdf','df')
#             #             df6.columns = ["alin_gA_"+_df6.columns[i] for i in range(len(df6.columns))]
#             #             df7 = pd.read_hdf(_filein[i]+'/output/rmsdbb_alin_prot_R_backbone.hdf','df')
#             #             df7.columns = ["alin_gpcr_"+_df7.columns[i] for i in range(len(df7.columns))]
#             #             df = pd.concat([df,df1,df2,df3,df4,df5,df6,df7], axis=1, ignore_index=False)
#             df = pd.concat([df, df1, df2, df3, df4, df5, df6], axis=1, ignore_index=False)
#
#         except:
#             print("Job not loaded because dataframe file was not found:",
#                   _filein[i].split("/")[-2])
#             print("\n")
#             continue
#
#         try:
#             assert (df.isnull().values.any() == False)
#             df_out = pd.concat([df_out, df])
#
#         except:
#             print("Job not loaded because null entries:",
#                   _filein[i].split("/")[-2])
#             print("\n")
#             continue
#
#     if not df_out.empty:
#         df_out[df_out["Time"] >= time_limit].reset_index(drop=True)
#
#     return df_out

# def read_data_hdf(receptor, ligand, time_limit=0, Verbose=False, exclude=[], include=[]):
#     """ Read in only hdf files
#
#     :param receptor: (str) string representing receptor name ex. 'd2gi'
#     :param ligand: (str) string representing ligand name ex. 'pd' or 'prm'
#     :param time_limit: (int) only include data where time is >= time_limit
#     :param Verbose: (bool) prints path of files to be read in if True
#     :param exclude: (sequence of str) jobnames to exclude
#     :param include: (sequence of str) list of hdf files to include
#     :return: dataframe containing all data from specified hdf files
#     """
#     _filein = glob.glob("/home/khlee/work/desmond/output/" + receptor + "/" + ligand + "/%s_*/" % receptor)
#     _filein = natsorted(_filein, key=lambda y: y.lower())
#     if Verbose:
#         print(_filein)
#     df_out = pd.DataFrame()
#     for i in range(len(_filein)):
#         if _filein[i].split("/")[-2] in exclude:
#             continue
#         try:
#             t = np.loadtxt(_filein[i] + '/ene/time.dat')
#             frame = list(range(1, len(t) + 1))
#
#         except:
#             print("Job not loaded because a file was not found:", _filein[i].split("/")[-2])
#             print("\n")
#             continue
#
#         try:
#             rec = [receptor for a in range(len(t))]
#             lig = [ligand for a in range(len(t))]
#             jobs = [_filein[i].split("/")[-2] for a in range(len(t))]
#             df = pd.DataFrame({"Receptor": rec,
#                                "Ligand": lig,
#                                "Jobname": jobs,
#                                "Time": t, "Frame": frame})
#
#
#
#         except:
#             print("Job not loaded because number of datapoints mismatch: ", _filein[i].split("/")[-2])
#             print("Time:", len(t))
#             print("\n")
#             continue
#
#         try:
#             if len(include) > 0:
#                 for inc in include:
#                     df2 = pd.read_hdf(_filein[i] + '/output/'+inc+'.hdf', 'df')
#                     if (len(t)!=len(df2)):
#                         print(str(inc)+": "+ str(len(df2)))
#                     df = pd.concat([df,df2], axis=1, ignore_index=False)
#
#         except:
#             print("Job not loaded because dataframe file was not found:",
#                   _filein[i].split("/")[-2])
#             print("\n")
#             continue
#
#         try:
#             assert (df.isnull().values.any() == False)
#             df_out = pd.concat([df_out, df])
#
#         except:
#             print("Job not loaded because null entries:",
#                   _filein[i].split("/")[-2])
#             print("\n")
#             continue
#
#     if not df_out.empty:
#         df_out[df_out["Time"] >= time_limit].reset_index(drop=True)
#
#     return df_out


def read_data2(receptor, ligand, time_limit=0, Verbose=False, exclude_ng=[], include_hdf=[], include_p=[]):
    """ Read in only hdf and pickle files

    :param receptor: (str) string representing receptor name ex. 'd2gi'
    :param ligand: (str) string representing ligand name ex. 'pd' or 'prm'
    :param time_limit: (int) only include data where time is >= time_limit
    :param Verbose: (bool) prints path of files to be read in if True
    :param exclude_ng: (sequence of str) jobnames to exclude
    :param include_hdf: (sequence of str) list of hdf files to include
    :param include_p: (sequence of str) list of hdf files to include
    :return: dataframe containing all data from specified hdf files
    """
    _filein = glob.glob("/home/khlee/work/desmond/output/" + receptor + "/" + ligand + "/%s_*/" % receptor)
    _filein = natsorted(_filein, key=lambda y: y.lower())
    if Verbose:
        print(_filein)
    df_out = pd.DataFrame()


    for i in range(len(_filein)):

        if _filein[i].split("/")[-2] in exclude_ng:
            continue

        try:
            t = np.loadtxt(_filein[i] + '/ene/time.dat')
            frame = list(range(1, len(t) + 1))

        except:
            print("Job not loaded because a file was not found:", _filein[i].split("/")[-2])
            print("\n")
            continue

        try:
            rec = [receptor for a in range(len(t))]
            lig = [ligand for a in range(len(t))]
            jobs = [_filein[i].split("/")[-2] for a in range(len(t))]
            df = pd.DataFrame({"Receptor": rec,
                               "Ligand": lig,
                               "Jobname": jobs,
                               "Time": t, "Frame": frame})

        except:
            print("Job not loaded because number of datapoints mismatch: ", _filein[i].split("/")[-2])
            print("Time:", len(t))
            print("\n")
            continue

        try:
            if len(include_p) > 0:
                for inc in include_p:
                    pin = pickle.load(open(_filein[i] + '/output/' + inc + '.p', "rb"))
                    df[inc] = pin

        except:
            print("Job not loaded because dataframe file was not found:",
                  _filein[i].split("/")[-2])
            print("\n")
            continue

        try:
            if len(include_hdf) > 0:
                for inc in include_hdf:
                    df2 = pd.read_hdf(_filein[i] + '/output/'+inc+'.hdf', 'df')
                    df2 = df2.reset_index()
                    del df2['index']
                    if (len(t)!=len(df2)):
                        print(str(inc)+": "+ str(len(df2)))
                    df = pd.concat([df,df2], axis=1, ignore_index=False)

        except:
            print("Job not loaded because dataframe file was not found:",
                  _filein[i].split("/")[-2])
            print("\n")
            continue



        try:
            assert (df.isnull().values.any() == False)
            df_out = pd.concat([df_out, df])

        except:
            print("Job not loaded because null entries:", _filein[i].split("/")[-2])
            print("\n")
            continue

    if not df_out.empty:
        df_out = df_out[df_out["Time"] >= time_limit].reset_index(drop=True)

    return df_out



def read_data_append(df, receptors, ligands, ls_exclude_ng, ls_includes_p, ls_includes_hdf, ctrl_time_limit):
    """ read dataframe df and append to existing df

    """
    for receptor in receptors:
        for ligand in ligands:
            _df = read_data2(receptor, ligand,
                             time_limit=ctrl_time_limit,
                             exclude_ng=ls_exclude_ng,
                             include_p=ls_includes_p,
                             include_hdf=ls_includes_hdf)
            df = pd.concat([df, _df], ignore_index=True)
    return df

# def read_data_hdf(receptor, ligand, time_limit=0, Verbose=False, exclude=[], include=[]):
#     """ Read in only hdf files
#
#     :param receptor: (str) string representing receptor name ex. 'd2gi'
#     :param ligand: (str) string representing ligand name ex. 'pd' or 'prm'
#     :param time_limit: (int) only include data where time is >= time_limit
#     :param Verbose: (bool) prints path of files to be read in if True
#     :param exclude: (sequence of str) jobnames to exclude
#     :param include: (sequence of str) list of hdf files to include
#     :return: dataframe containing all data from specified hdf files
#     """
#     _filein = glob.glob("/home/khlee/work/desmond/output/" + receptor + "/" + ligand + "/%s_*/" % receptor)
#     _filein = natsorted(_filein, key=lambda y: y.lower())
#     if Verbose:
#         print(_filein)
#     df_out = pd.DataFrame()
#     for i in range(len(_filein)):
#         if _filein[i].split("/")[-2] in exclude:
#             continue
#         try:
#             t = np.loadtxt(_filein[i] + '/ene/time.dat')
#             frame = list(range(1, len(t) + 1))
#
#         except:
#             print("Job not loaded because a file was not found:", _filein[i].split("/")[-2])
#             print("\n")
#             continue
#
#         try:
#             rec = [receptor for a in range(len(t))]
#             lig = [ligand for a in range(len(t))]
#             jobs = [_filein[i].split("/")[-2] for a in range(len(t))]
#             df = pd.DataFrame({"Receptor": rec,
#                                "Ligand": lig,
#                                "Jobname": jobs,
#                                "Time": t, "Frame": frame})
#
#
#
#         except:
#             print("Job not loaded because number of datapoints mismatch: ", _filein[i].split("/")[-2])
#             print("Time:", len(t))
#             print("\n")
#             continue
#
#         try:
#             if len(include) > 0:
#                 for inc in include:
#                     df2 = pd.read_hdf(_filein[i] + '/output/'+inc+'.hdf', 'df')
#                     if (len(t)!=len(df2)):
#                         print(str(inc)+": "+ str(len(df2)))
#                     df = pd.concat([df,df2], axis=1, ignore_index=False)
#
#         except:
#             print("Job not loaded because dataframe file was not found:",
#                   _filein[i].split("/")[-2])
#             print("\n")
#             continue
#
#         try:
#             assert (df.isnull().values.any() == False)
#             df_out = pd.concat([df_out, df])
#
#         except:
#             print("Job not loaded because null entries:",
#                   _filein[i].split("/")[-2])
#             print("\n")
#             continue
#
#     if not df_out.empty:
#         df_out[df_out["Time"] >= time_limit].reset_index(drop=True)
#
#     return df_out
#
#
# def read_data_p(receptor, ligand, time_limit=0, Verbose=False, exclude=[], include=[]):
#     """ Read in only pickle files
#
#     :param receptor: (str) string representing receptor name ex. 'd2gi'
#     :param ligand: (str) string representing ligand name ex. 'pd' or 'prm'
#     :param time_limit: (int) only include data where time is >= time_limit
#     :param Verbose: (bool) prints path of files to be read in if True
#     :param exclude: (sequence of str) jobnames to exclude
#     :param include: (sequence of str) list of hdf files to include
#     :return: dataframe containing all data from specified hdf files
#     """
#     _filein = glob.glob("/home/khlee/work/desmond/output/" + receptor + "/" + ligand + "/%s_*/" % receptor)
#     _filein = natsorted(_filein, key=lambda y: y.lower())
#     if Verbose:
#         print(_filein)
#     df_out = pd.DataFrame()
#     for i in range(len(_filein)):
#         if _filein[i].split("/")[-2] in exclude:
#             continue
#         try:
#             t = np.loadtxt(_filein[i] + '/ene/time.dat')
#             frame = list(range(1, len(t) + 1))
#
#         except:
#             print("Job not loaded because a file was not found:", _filein[i].split("/")[-2])
#             print("\n")
#             continue
#
#         try:
#             rec = [receptor for a in range(len(t))]
#             lig = [ligand for a in range(len(t))]
#             jobs = [_filein[i].split("/")[-2] for a in range(len(t))]
#             df = pd.DataFrame({"Receptor": rec,
#                                "Ligand": lig,
#                                "Jobname": jobs,
#                                "Time": t, "Frame": frame})
#
#
#
#         except:
#             print("Job not loaded because number of datapoints mismatch: ", _filein[i].split("/")[-2])
#             print("Time:", len(t))
#             print("\n")
#             continue
#
#         try:
#             if len(include) > 0:
#                 for inc in include:
#                     # df2 = pd.read_hdf(_filein[i] + '/output/'+inc+'.hdf', 'df')
#                     pin = pickle.load(open(_filein[i] + '/output/' + inc + '.p', "rb"))
#                     df[inc] = pin
#
#         except:
#             print("Job not loaded because dataframe file was not found:",
#                   _filein[i].split("/")[-2])
#             print("\n")
#             continue
#
#         try:
#             assert (df.isnull().values.any() == False)
#             df_out = pd.concat([df_out, df])
#
#         except:
#             print("Job not loaded because null entries:",
#                   _filein[i].split("/")[-2])
#             print("\n")
#             continue
#
#     if not df_out.empty:
#         df_out[df_out["Time"] >= time_limit].reset_index(drop=True)
#
#     return df_out


def rename_df(df):
    # extra copy of original Receptor, Ligand, and Jobname
    # df['Receptor']=df['Receptor_Ori']
    # df['Ligand']=df['Ligand_Ori']
    # df['Jobname']=df['Jobname_Ori']

    # extra copy of original Receptor, Ligand, and Jobname
    df['Receptor_Ori'] = df['Receptor']
    df['Ligand_Ori'] = df['Ligand']
    df['Jobname_Ori'] = df['Jobname']

    # rename drd3 -> d3gi
    # rename drd3ao -> d3go
    # rename prm7,prm3 -> prm & pd2 -> pd
    df.loc[df['Receptor_Ori'] == 'drd3', 'Receptor'] = 'd3gi'
    df.loc[df['Receptor_Ori'] == 'drd3ao', 'Receptor'] = 'd3go'
    df.loc[df['Ligand_Ori'] == 'prm7', 'Ligand'] = 'prm'
    df.loc[df['Ligand_Ori'] == 'prm3', 'Ligand'] = 'prm'
    df.loc[df['Ligand_Ori'] == 'pd2', 'Ligand'] = 'pd'

    # rename jobnames
    # ls_xx = [xx for xx in set(df['Jobname_Ori']) if xx.startswith('drd3')]
    # ls_xx2 = [xx.replace("drd3ao","d3go").replace("drd3","d3gi").replace("prm7","prm").replace("prm3","prm").replace("pd2","pd")
    #           + ".3e" for xx in ls_xx]
    ls_xx = [xx for xx in set(df['Jobname_Ori']) if xx.startswith('drd3')]
    ls_xx2 = [
        xx.replace("drd3ao", "d3go").replace("drd3", "d3gi").replace("prm7", "prm").replace("prm3", "prm").replace(
            "pd2", "pd")
        + ".3e" for xx in ls_xx]
    for i in range(len(ls_xx)):
        df.loc[df[df['Jobname_Ori'].str.contains(ls_xx[i])].index, 'Jobname'] = ls_xx2[i]
    return df


def do_extra_calculation_df(df):
    # extra calculation for interface ratio
    df['ratio'] = round(df['gpcr_interface_nopolyg'] / df['gpcr_interface'] * 100, 2)
    # extra calculation for interaction evergy
    df['interaction_all'] = df['elec_all'] + df['vdw_all']
    df['interaction_no_il3'] = df['elec_no_il3'] + df['vdw_no_il3']
    df['interaction_il1'] = df['elec_il1'] + df['vdw_il1']
    df['interaction_il2'] = df['elec_il2'] + df['vdw_il2']
    df['interaction_il3'] = df['elec_il3'] + df['vdw_il3']
    df['interaction_tm1i'] = df['elec_tm1i'] + df['vdw_tm1i']
    df['interaction_tm2i'] = df['elec_tm2i'] + df['vdw_tm2i']
    df['interaction_tm3i'] = df['elec_tm3i'] + df['vdw_tm3i']
    df['interaction_tm4i'] = df['elec_tm4i'] + df['vdw_tm4i']
    df['interaction_tm5i'] = df['elec_tm5i'] + df['vdw_tm5i']
    df['interaction_tm6i'] = df['elec_tm6i'] + df['vdw_tm6i']
    df['interaction_tm7i-H8'] = df['elec_tm7i-H8'] + df['vdw_tm7i-H8']
    return df




def get_status(df):
    """
    a quick summary of number of simulation status using information from dataframe df
    """
    ls_sim = []
    ls_length = []
    _ls_jobnames = natsorted(list(set(df['Jobname'])), key=lambda y: y.lower(), reverse=False)
    for ii in _ls_jobnames:
        ls_sim.append(ii)
        ls_length.append(max(df[df['Jobname'] == ii]['Time']))
    df_sum = pd.DataFrame({'Jobname': ls_sim, 'Length': ls_length})

    print("d3gi")
    #print(df_sum.loc[df_sum['Jobname'].str.contains("d3gi_", case=False)].sort_values('Jobname').to_string(index=False))
    print(df_sum.loc[df_sum['Jobname'].str.contains("d3gi_", case=False)].to_string(index=False))
    print("\n")
    print("d3go")
    print(df_sum.loc[df_sum['Jobname'].str.contains("d3go_", case=False)].to_string(index=False))
    print("\n")
    print("d2gi")
    #print(df_sum.loc[df_sum['Jobname'].str.contains("d2gi_", case=False)].sort_values('Jobname').to_string(index=False))
    print(df_sum.loc[df_sum['Jobname'].str.contains("d2gi_", case=False)].to_string(index=False))
    print("\n")
    print("d2go")
    print(df_sum.loc[df_sum['Jobname'].str.contains("d2go_", case=False)].to_string(index=False))
    # print(len(set(df['Jobname'])))


def split_prodrun(df):
    """ Split dataframe into production and equilibrium sets

    :param df: (pandas dataframe) dataframe containing all input data obtained using read_data or read_caller
    :return: 2 dataframes -- equilibrium data, production data
    """
    df_prod = df[df['Jobname'].str.contains(".f\d")]
    i = list(df_prod.index)
    df = df.drop(i)
    return df.reset_index(drop=True), df_prod.reset_index(drop=True)


def check_interface(df):
    """ Perform analysis on interface surface area difference and append results to dataframe

    :param df: (pandas dataframe) dataframe with interface data
    :return: pandas dataframe with 3 additional columns
    """
    df["RAB-AB"] = abs(df["Interface"] - (df["Interface_A"] + df["Interface_B"]))
    df["AB"] = df["Interface_A"] + df["Interface_B"]
    df["%Diff"] = abs(df["RAB-AB"] / df["Interface"] * 100)
    return df


def read_caller(keywords, exclude=[], include=[], Verbose=False):
    """ Used to call read_data function for certain pairs of receptor/ligand

    :param keywords: (sequence of tuples of str) list of pairs of receptor and ligands to read in ex. [(d3gi, pd)]
    :param exclude: (sequence of str) list of jobnames to exclude
    :param include: (sequence of str) list of pickle files to include
    :param Verbose: (bool) prints path of files to be read in if True
    :return: dataframe containing all input data
    """
    df = pd.DataFrame()
    for pair in keywords:
        df = pd.concat([df, read_data(pair[0], pair[1], exclude=exclude, include=include, Verbose=Verbose)])
    return df.reset_index(drop=True)


def gen_random_picks(control_seed=2020, num_splits=50, max_value=10000, Verbose=False):
    """ Generate random numbers

    :param control_seed: (int) random seed
    :param num_splits: (int) number of random numbers to generate
    :param max_value: (int) max limit of random number to generate
    :param Verbose: (bool) print random numbers generated if True
    :return: list of generated random numbers
    """
    random_splits = []
    np.random.seed(control_seed)

    for a in range(0, num_splits):
        new_random = np.random.randint(0, max_value)
        # this part should be removed. we want them have chance to be repeatedfrom matplotlib import gridspec
        # while (new_random in random_splits_old) or (new_random in random_splits):
        # while (new_random in random_splits):
        #     # ensure no duplicate
        #     new_random = np.random.randint(0, max_value)
        random_splits.append(new_random)

    # random_splits
    if Verbose:
        for random_split in random_splits:
            print(random_split)
    return random_splits


def bootstrap(indf, splits=1, n=1000, seed=2020):
    """

    :param indf: (pandas dataframe) input dataframe to perform bootstrapping
    :param splits: (int) number of times perform data selection
    :param n: (int) number of entries to randomly select, passed to gen_random_picks function
    :param seed: (int) random seed
    :return: bootstrapped dataframe
    """
    recs = list(set(indf["Receptor"]))
    sel = pd.DataFrame()
    _splits = gen_random_picks(seed, num_splits=splits)
    for rec in recs:
        _df_rec = indf[indf["Receptor"] == rec].reset_index(drop=True)
        ligs = list(set(_df_rec["Ligand"]))
        for lig in ligs:
            _df_rec_lig = _df_rec[_df_rec["Ligand"] == lig].reset_index(drop=True)
            for i in range(splits):
                r = gen_random_picks(_splits[i], num_splits=n, max_value=len(_df_rec_lig))
                sel1 = _df_rec_lig.iloc[r].reset_index(drop=True)
                sel1["Split"] = [i + 1 for j in range(len(sel1))]
                sel = pd.concat([sel, sel1], ignore_index=True)
    return sel


def bootstrap_stats(indf, col="Interface", by_receptor=False):
    """ Get the statistics from bootstrapped dataframe

    :param indf:
    :param col:
    :return:
    """
    recs = list(set(indf["Receptor"]))
    data = {"Split": [], "Jobname": [], col: []}
    if by_receptor==True:
        for rec in recs:
            temp = indf[indf["Receptor"] == recs[0]].reset_index(drop=True)
            splits = list(set(temp["Split"]))
            for split in splits:
                df = temp[temp["Split"] == split]
                data["Split"].append(split)
                data["Jobname"].append(rec)
                data[col].append(statistics.mean(df[col]))
    else:
        for rec in recs:
            temp = indf[indf["Receptor"] == rec].reset_index(drop=True)
            ligs = list(set(temp["Ligand"]))
            for lig in ligs:
                temp2 = temp[temp["Ligand"] == lig].reset_index(drop=True)
                splits = list(set(temp2["Split"]))
                for split in splits:
                    df = temp2[temp2["Split"] == split]
                    data["Split"].append(split)
                    data["Jobname"].append(rec + "_" + lig)
                    data[col].append(statistics.mean(df[col]))
    return pd.DataFrame(data)



def get_seg_lists(df):
    ls_seg = []
    for i in range(len(df.columns)):
        if "seg_dist_" in df.columns[i]:
            ls_seg.append(df.columns[i])
        elif "seg_angle_" in df.columns[i]:
            ls_seg.append(df.columns[i])

    segi = []
    segj = []
    for i in range(len(ls_seg)):
        segi.append(ls_seg[i].split("_")[2])
        segj.append(ls_seg[i].split("_")[3])
    segj = list(set(segj))
    segi = list(set(segi))
    segi = natsorted(segi, key=lambda y: y.lower())
    segj = natsorted(segj, key=lambda y: y.lower())

    return segi, segj


def get_matrix(df, receptor, ligand, sel, by_receptor=False):
#def get_matrix(df, receptor, ligand, sel):
    # average of all
    segi, segj = get_seg_lists(df)
    matrix = np.zeros((len(segi), len(segj)))
    if sel == "distance":
        sel = "seg_dist_"
    elif sel == "angle":
        sel = "seg_angle_"

    df = df[df["Receptor"] == receptor]
    if by_receptor==False:
        df = df[df["Ligand"] == ligand]
    for i in range(len(segi)):
        for j in range(len(segj)):
            matrix[i][j] = statistics.mean(df[sel + segi[i] + "_" + segj[j]])

    return matrix


def get_matrix_dictionary(dfbs, sels=["seg_dist_", "seg_angle_"]):
    matrices = {}
    for sel in sels:
        for receptor in set(dfbs["Receptor"]):
            dfbs_sel = dfbs[dfbs["Receptor"] == receptor]
            for ligand in set(dfbs_sel["Ligand"]):
                matrices[receptor + "_" + ligand] = get_matrix(dfbs, receptor, ligand, sel)
    return matrices


def get_matrix_diff(matrix_pair1, matrix_pair2):
    '''
    :param matrix_pair1: input matrix #1
    :param matrix_pair2: input matrix #2
    :return:  matrix difference
    '''
    nshape0 = matrix_pair1.shape[0]
    nshape1 = matrix_pair1.shape[1]
    matrix_diff = np.zeros((nshape0, nshape1))
    # calculate difference
    for i in range(nshape0):
        for j in range(nshape1):
            matrix_diff[i][j] = matrix_pair1[i][j] - matrix_pair2[i][j]
    return matrix_diff

def get_binary_matrix(matrix, cutoff=0):
    '''
    convert matrix into binary matrix
    :param matrix:
    :param cutoff:
    :return:
    '''
    nshape0 = matrix.shape[0]
    nshape1 = matrix.shape[1]
    matrixout = np.zeros((nshape0, nshape1))
    # based on cutoff to convert matrix into three classes (-1,0,1)
    for i in range(nshape0):
        for j in range(nshape1):
            if matrix[i][j] > cutoff:
                matrixout[i][j] = 1
            elif matrix[i][j] < -cutoff:
                matrixout[i][j] = -1
            else:
                matrixout[i][j] = 0
    return matrixout


# def add_matrixes(matrixA, matrixB):
#     nshape0 = matrixA.shape[0]
#     nshape1 = matrixA.shape[1]
#     matrixout = np.zeros((nshape0, nshape1))
#     for i in range(nshape0):
#         for j in range(nshape1):
#             matrixout[i][j] = matrixA[i][j] + matrixB[i][j]
#     return matrixout



def get_ratio_matrix_diff(matrix_freq_1, matrix_freq_2):
    matrix_i = matrix_freq_1.shape[0]
    matrix_j = matrix_freq_1.shape[1]
    matrix_diff = np.zeros((matrix_i, matrix_j))
    # calculate distance matrix
    for i in range(matrix_i):
        for j in range(matrix_j):
            mean_dists = (matrix_freq_1[i][j] + matrix_freq_2[i][j])/2
            if abs(matrix_freq_1[i][j] - matrix_freq_2[i][j]) > 0.5:
                matrix_diff[i][j] = (matrix_freq_1[i][j] - matrix_freq_2[i][j])/mean_dists*100
            else:
                matrix_diff[i][j] = 0
    return matrix_diff

def get_filtered_ratio_matrix_diff(matrix_freq_1, matrix_freq_2):
    matrix_i = matrix_freq_1.shape[0]
    matrix_j = matrix_freq_1.shape[1]
    matrix_diff = np.zeros((matrix_i, matrix_j))
    # calculate distance matrix
    for i in range(matrix_i):
        for j in range(matrix_j):
            mean_dists = (matrix_freq_1[i][j] + matrix_freq_2[i][j])/2
            if abs(matrix_freq_1[i][j] - matrix_freq_2[i][j]) > 0.5:
                matrix_diff[i][j] = (matrix_freq_1[i][j] - matrix_freq_2[i][j])/mean_dists*100
                if abs(matrix_diff[i][j]) < 10:
                    matrix_diff[i][j] = 0
            else:
                matrix_diff[i][j] = 0
    return matrix_diff







'''
Correlation matrix calculation
_calc_corr
get_corr_mat
'''


def calc_corr(i, j, indf):
    A=pd.Series(indf[indf.columns[i]])
    B=pd.Series(indf[indf.columns[j]])
    xx = round(A.corr(B),2)
    return xx


def get_corr_mat(dfsel2):

    nn=len(dfsel2.columns)

    _ls = []
    #for i in range(nn): [_ls.append((i,j,dfsel2)) for j in range(i + 1, nn)]
    for i in range(nn): [_ls.append((i,j,dfsel2)) for j in range(nn)]
    print("generating list")
    print(len(_ls))

    #pool = Pool(ncpus)
    pool = Pool(32)
    _results = pool.starmap(calc_corr, _ls)
    pool.close()
    pool.join()
    assert _results.count(-1)==0   # ensure all pairwise rmsd is calculated
    #resample matrix
    print("resample matrix")
    matrix = np.zeros((nn, nn))
    k = 0
    for i in range(nn):
        #for j in range(i + 1, nn):
        for j in range(nn):
            matrix[i][j] = float(_results[k])
            k = k + 1
    return matrix




def write_dict_pickle(in_dict, out_filename):
    with open(out_filename+'.p', 'wb') as f:
        pickle.dump(in_dict, f)

def read_dict_pickle(in_filename):
    with open(in_filename+'.p', 'rb') as f:
        out_dict = pickle.load(f)
    return out_dict

def write_df_hdf(in_df, out_filename):
    in_df.to_hdf(out_filename+'.hdf', key='df', mode='w')



from natsort import natsorted
import numpy as np
from statistics import mean, stdev


def get_abs_chi_diff(chi1_A, chi1_B):
    if len(chi1_A) != len(chi1_B):
        print("Check data length")
        return ""
    ls_diff = []
    for i in range(len(chi1_A)):
        chi1_diff_abs = abs(chi1_A[i]-chi1_B[i])
        if chi1_diff_abs > 180:
            chi1_diff_abs = 360 - chi1_diff_abs
        chi1_diff_abs = round(chi1_diff_abs,1)
        chi1_diff_abs = abs(chi1_diff_abs)
        if chi1_diff_abs > 9000:
            chi1_diff_abs = 9999
        ls_diff.append(chi1_diff_abs)
    return ls_diff


def get_obs_sbp_chi1_mean_stdev(_df, OBSorSBP, obs_sbp="obs"):
    if obs_sbp == "obs":
        dflabel = "piabs_chi1_"
    elif obs_sbp == "sbp":
        dflabel = "piabs_chi1_"
    else:
        print("type is obs or sbp only")

    for iOBSorSBP in OBSorSBP:
        _sel_colname = str(dflabel) + str(iOBSorSBP)
        _sel_colname_ori = _sel_colname + "_ori"
        _df[_sel_colname_ori] = _df[_sel_colname]
        # convert -180~0 to 180~360 for doing the average
        for irow in range(len(_df[_sel_colname])):
            if isinstance(_df[_sel_colname_ori][irow], str):
                print("warning: ", iOBSorSBP, _df[_sel_colname_ori][irow])
                _df[_sel_colname][irow] = 9999
            else:
                if (float(_df[_sel_colname_ori][irow]) < 0):
                    _df[_sel_colname][irow] = float(_df[_sel_colname_ori][irow]) + 360

    ls_mean = []
    ls_stdev = []
    for iOBSorSBP in OBSorSBP:
        _sel_colname = str(dflabel) + str(iOBSorSBP)
        chi1_mean = mean(_df[_sel_colname])
        if chi1_mean != 9999:
            if chi1_mean > 180:
                chi1_mean = chi1_mean - 360
        chi1_mean = round(chi1_mean, 1)
        ls_mean.append(chi1_mean)
        try:
            ls_stdev.append(round(stdev(_df[_sel_colname]), 1))
        except:
            ls_stdev.append(0)
    return ls_mean, ls_stdev


def get_pia_bs_chi1(_df, OBSorSBP):
    dflabel = "piabs_chi1_"

    for iOBSorSBP in OBSorSBP:
        _sel_colname = str(dflabel) + str(iOBSorSBP)
        _sel_colname_ori = _sel_colname + "_ori"
        _df[_sel_colname_ori] = _df[_sel_colname]
        # convert -180~0 to 180~360 for doing the average
        for irow in range(len(_df[_sel_colname])):
            if isinstance(_df[_sel_colname_ori][irow], str):
                print("warning: ", iOBSorSBP, _df[_sel_colname_ori][irow])
                _df[_sel_colname][irow] = 9999
            else:
                if (float(_df[_sel_colname_ori][irow]) < 0):
                    _df[_sel_colname][irow] = float(_df[_sel_colname_ori][irow]) + 360

    ls_mean = []
    for iOBSorSBP in OBSorSBP:
        _sel_colname = str(dflabel) + str(iOBSorSBP)
        chi1_mean = mean(_df[_sel_colname])
        if chi1_mean != 9999:
            if chi1_mean > 180:
                chi1_mean = chi1_mean - 360
        chi1_mean = round(chi1_mean, 1)
        ls_mean.append(chi1_mean)
    return ls_mean


def get_obs_dist_matrix(_df, OBS):
    '''
    give df & OBS or SBP
    return matrix
    '''
    matrix = np.zeros((len(OBS), len(OBS)))
    for i in range(len(OBS)):
        for j in range(i+1,len(OBS)):
            _sel_colname_ca="piabs_cadist_"+str(OBS[i])+"_"+str(OBS[j])
            matrix[i][j] = round(mean(_df[_sel_colname_ca]), 2)
            matrix[j][i] = round(mean(_df[_sel_colname_ca]),2)
    return matrix

def get_sbp_dist_matrix(_df, SBP):
    '''
    give df & SBP
    return matrix
    '''
    matrix = np.zeros((len(SBP), len(SBP)))
    for i in range(len(SBP)):
        for j in range(i+1,len(SBP)):
            _sel_colname_ca="piabs_cadist_"+str(SBP[i])+"_"+str(SBP[j])
            matrix[i][j] = round(mean(_df[_sel_colname_ca]), 2)
            matrix[j][i] = round(mean(_df[_sel_colname_ca]),2)
    return matrix


def get_seg_lists(df):
    # given df
    # return unique elements

    _segs_i = []
    _segs_j = []
    for i in df.columns:
        if "seg_dist_" in i:
            _pair = i.replace("seg_dist_", "").split("_")
            _segs_i.append(_pair[0])
            _segs_j.append(_pair[1])
        elif "seg_angle_" in i:
            _pair = i.replace("seg_angle_", "").split("_")
            _segs_i.append(_pair[0])
            _segs_j.append(_pair[1])

    _segs_i = list(set(_segs_i))
    _segs_j = list(set(_segs_j))
    _segs_i = natsorted(_segs_i, key=lambda y: y.lower())
    _segs_j = natsorted(_segs_j, key=lambda y: y.lower())
    return _segs_i, _segs_j


def get_matrix(df, receptor, ligand, sel, by_receptor=False):
    # average of all
    segi, segj = get_seg_lists(df)
    matrix = np.zeros((len(segi), len(segj)))

    if sel == "distance":
        sel = "seg_dist_"
    elif sel == "angle":
        sel = "seg_angle_"

    df = df[df["Receptor"] == receptor]
    if by_receptor==False:
        df = df[df["Ligand"] == ligand]
    for i in range(len(segi)):
        for j in range(len(segj)):
            matrix[i][j] = mean(df[sel + segi[i] + "_" + segj[j]])

    return matrix


def convert_segment_pia_df_to_matrix(_indf):
    segi, segj = get_seg_lists(_indf)
    matrix = np.zeros((len(segi), len(segj)))
    for i in range(len(segi)):
        for j in range(len(segj)):
            try:
                try:
                    matrix[i][j] = mean(_indf["seg_dist_" + segi[i] + "_" + segj[j]])
                except:
                    matrix[i][j] = 0
            except:
                try:
                    matrix[i][j] = mean(_indf["seg_angle_" + segi[i] + "_" + segj[j]])
                except:
                    matrix[i][j] = 0

    return matrix






