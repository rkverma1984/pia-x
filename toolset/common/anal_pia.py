# from natsort import natsorted
# import numpy as np
# from statistics import mean, stdev
#
#
# def get_abs_chi_diff(chi1_A, chi1_B):
#     if len(chi1_A) != len(chi1_B):
#         print("Check data length")
#         return ""
#     ls_diff = []
#     for i in range(len(chi1_A)):
#         chi1_diff_abs = abs(chi1_A[i]-chi1_B[i])
#         if chi1_diff_abs > 180:
#             chi1_diff_abs = 360 - chi1_diff_abs
#         chi1_diff_abs = round(chi1_diff_abs,1)
#         chi1_diff_abs = abs(chi1_diff_abs)
#         if chi1_diff_abs > 9000:
#             chi1_diff_abs = 9999
#         ls_diff.append(chi1_diff_abs)
#     return ls_diff
#
#
# def get_obs_sbp_chi1_mean_stdev(_df, OBSorSBP, obs_sbp="obs"):
#     if obs_sbp == "obs":
#         dflabel = "piabs_chi1_"
#     elif obs_sbp == "sbp":
#         dflabel = "piabs_chi1_"
#     else:
#         print("type is obs or sbp only")
#
#     for iOBSorSBP in OBSorSBP:
#         _sel_colname = str(dflabel) + str(iOBSorSBP)
#         _sel_colname_ori = _sel_colname + "_ori"
#         _df[_sel_colname_ori] = _df[_sel_colname]
#         # convert -180~0 to 180~360 for doing the average
#         for irow in range(len(_df[_sel_colname])):
#             if isinstance(_df[_sel_colname_ori][irow], str):
#                 print("warning: ", iOBSorSBP, _df[_sel_colname_ori][irow])
#                 _df[_sel_colname][irow] = 9999
#             else:
#                 if (float(_df[_sel_colname_ori][irow]) < 0):
#                     _df[_sel_colname][irow] = float(_df[_sel_colname_ori][irow]) + 360
#
#     ls_mean = []
#     ls_stdev = []
#     for iOBSorSBP in OBSorSBP:
#         _sel_colname = str(dflabel) + str(iOBSorSBP)
#         chi1_mean = mean(_df[_sel_colname])
#         if chi1_mean != 9999:
#             if chi1_mean > 180:
#                 chi1_mean = chi1_mean - 360
#         chi1_mean = round(chi1_mean, 1)
#         ls_mean.append(chi1_mean)
#         try:
#             ls_stdev.append(round(stdev(_df[_sel_colname]), 1))
#         except:
#             ls_stdev.append(0)
#     return ls_mean, ls_stdev
#
#
# def get_pia_bs_chi1(_df, OBSorSBP):
#     dflabel = "piabs_chi1_"
#
#     for iOBSorSBP in OBSorSBP:
#         _sel_colname = str(dflabel) + str(iOBSorSBP)
#         _sel_colname_ori = _sel_colname + "_ori"
#         _df[_sel_colname_ori] = _df[_sel_colname]
#         # convert -180~0 to 180~360 for doing the average
#         for irow in range(len(_df[_sel_colname])):
#             if isinstance(_df[_sel_colname_ori][irow], str):
#                 print("warning: ", iOBSorSBP, _df[_sel_colname_ori][irow])
#                 _df[_sel_colname][irow] = 9999
#             else:
#                 if (float(_df[_sel_colname_ori][irow]) < 0):
#                     _df[_sel_colname][irow] = float(_df[_sel_colname_ori][irow]) + 360
#
#     ls_mean = []
#     for iOBSorSBP in OBSorSBP:
#         _sel_colname = str(dflabel) + str(iOBSorSBP)
#         chi1_mean = mean(_df[_sel_colname])
#         if chi1_mean != 9999:
#             if chi1_mean > 180:
#                 chi1_mean = chi1_mean - 360
#         chi1_mean = round(chi1_mean, 1)
#         ls_mean.append(chi1_mean)
#     return ls_mean
#
#
# def get_obs_dist_matrix(_df, OBS):
#     '''
#     give df & OBS or SBP
#     return matrix
#     '''
#     matrix = np.zeros((len(OBS), len(OBS)))
#     for i in range(len(OBS)):
#         for j in range(i+1,len(OBS)):
#             _sel_colname_ca="piabs_cadist_"+str(OBS[i])+"_"+str(OBS[j])
#             matrix[i][j] = round(mean(_df[_sel_colname_ca]), 2)
#             matrix[j][i] = round(mean(_df[_sel_colname_ca]),2)
#     return matrix
#
# def get_sbp_dist_matrix(_df, SBP):
#     '''
#     give df & SBP
#     return matrix
#     '''
#     matrix = np.zeros((len(SBP), len(SBP)))
#     for i in range(len(SBP)):
#         for j in range(i+1,len(SBP)):
#             _sel_colname_ca="piabs_cadist_"+str(SBP[i])+"_"+str(SBP[j])
#             matrix[i][j] = round(mean(_df[_sel_colname_ca]), 2)
#             matrix[j][i] = round(mean(_df[_sel_colname_ca]),2)
#     return matrix
#
#
# def get_seg_lists(df):
#     # given df
#     # return unique elements
#
#     _segs_i = []
#     _segs_j = []
#     for i in df.columns:
#         if "seg_dist_" in i:
#             _pair = i.replace("seg_dist_", "").split("_")
#             _segs_i.append(_pair[0])
#             _segs_j.append(_pair[1])
#         elif "seg_angle_" in i:
#             _pair = i.replace("seg_angle_", "").split("_")
#             _segs_i.append(_pair[0])
#             _segs_j.append(_pair[1])
#
#     _segs_i = list(set(_segs_i))
#     _segs_j = list(set(_segs_j))
#     _segs_i = natsorted(_segs_i, key=lambda y: y.lower())
#     _segs_j = natsorted(_segs_j, key=lambda y: y.lower())
#     return _segs_i, _segs_j
#
#
# def get_matrix(df, receptor, ligand, sel, by_receptor=False):
#     # average of all
#     segi, segj = get_seg_lists(df)
#     matrix = np.zeros((len(segi), len(segj)))
#
#     if sel == "distance":
#         sel = "seg_dist_"
#     elif sel == "angle":
#         sel = "seg_angle_"
#
#     df = df[df["Receptor"] == receptor]
#     if by_receptor==False:
#         df = df[df["Ligand"] == ligand]
#     for i in range(len(segi)):
#         for j in range(len(segj)):
#             matrix[i][j] = mean(df[sel + segi[i] + "_" + segj[j]])
#
#     return matrix
#
#
# def convert_segment_pia_df_to_matrix(_indf):
#     segi, segj = get_seg_lists(_indf)
#     matrix = np.zeros((len(segi), len(segj)))
#     for i in range(len(segi)):
#         for j in range(len(segj)):
#             matrix[i][j] = mean(_indf["seg_dist_" + segi[i] + "_" + segj[j]])
#     return matrix
#
#
#
#
#
