# import MDAnalysis as mda
# import pandas as pd
# from natsort import natsorted
# import from toolset
import sys
from pathlib import Path
home = str(Path.home())
toolset_dir = home+'/repositories/pia-x/toolset'
toolset_common = toolset_dir+"/common"
toolset_gpcr = toolset_dir+"/gpcr"
sys.path.insert(0, toolset_dir)
sys.path.insert(0, toolset_common)
sys.path.insert(0, toolset_gpcr)
from common_gpcr import *
from heatmapplus import *



d2table = convert_Residue2TMs("DRD2")
d3table = convert_Residue2TMs("DRD3")

d3loop = {
58: 'IL1_58_63',
59: 'IL1_59_64',
60: 'IL1_60_65',
61: 'IL1_61_66',
62: 'IL1_62_67',
135: 'IL2_135_139',
136: 'IL2_136_140',
137: 'IL2_137_141',
138: 'IL2_138_142',
139: 'IL2_139_143',
140: 'IL2_140_x',
141: 'IL2_141_x',
142: 'IL2_142_144',
143: 'IL2_143_145',
144: 'IL2_144_146',
219: 'IL3_219_220',
220: 'IL3_220_221',
221: 'IL3_221_222',
222: 'IL3_222_223',
223: 'IL3_223_224',
224: 'IL3_224_225',
225: 'IL3_225_226',
226: 'IL3_226_227',
227: 'IL3_227_228',
228: 'IL3_228_229',
229: 'IL3_229_230_GLY',
230: 'IL3_230_231_GLY',
231: 'IL3_231_232_GLY',
232: 'IL3_232_233_GLY',
233: 'IL3_233_234_GLY',
234: 'IL3_234_235_GLY',
235: 'IL3_235_236_GLY',
236: 'IL3_236_237_GLY',
237: 'IL3_237_238_GLY',
316: 'IL3_316_360',
317: 'IL3_317_361',
318: 'IL3_318_362',
319: 'IL3_319_363',
320: 'IL3_320_364',
321: 'IL3_321_365'
}


d2loop = {
63:'IL1_58_63',
64:'IL1_59_64',
65:'IL1_60_65',
66:'IL1_61_66',
67:'IL1_62_67',
139:'IL2_135_139',
140:'IL2_136_140',
141:'IL2_137_141',
142:'IL2_138_142',
143:'IL2_139_143',
144:'IL2_142_144',
145:'IL2_143_145',
146:'IL2_144_146',
220:'IL3_219_220',
221:'IL3_220_221',
222:'IL3_221_222',
223:'IL3_222_223',
224:'IL3_223_224',
225:'IL3_224_225',
226:'IL3_225_226',
227:'IL3_226_227',
228:'IL3_227_228',
229:'IL3_228_229',
230:'IL3_229_230_GLY',
231:'IL3_230_231_GLY',
232:'IL3_231_232_GLY',
233:'IL3_232_233_GLY',
234:'IL3_233_234_GLY',
235:'IL3_234_235_GLY',
236:'IL3_235_236_GLY',
237:'IL3_236_237_GLY',
238:'IL3_237_238_GLY',
360:'IL3_316_360',
361:'IL3_317_361',
362:'IL3_318_362',
363:'IL3_319_363',
364:'IL3_320_364',
365:'IL3_321_365'
}


gi_table = {
1:"M_1_M_1",
2:"G_2_G_2",
3:"C_3_C_3",
4:"T_4_T_4",
5:"L_5_L_5",
6:"S_6_S_6",
7:"A_7_A_7",
8:"E_8_E_8",
9:"D_9_E_9",
10:"K_10_R_10",
11:"A_11_A_11",
12:"A_12_A_12",
13:"V_13_L_13",
14:"E_14_E_14",
15:"R_15_R_15",
16:"S_16_S_16",
17:"K_17_K_17",
18:"M_18_A_18",
19:"I_19_I_19",
20:"D_20_E_20",
21:"R_21_K_21",
22:"N_22_N_22",
23:"L_23_L_23",
24:"R_24_K_24",
25:"E_25_E_25",
26:"D_26_D_26",
27:"G_27_G_27",
28:"E_28_I_28",
29:"K_29_S_29",
30:"A_30_A_30",
31:"A_31_A_31",
32:"R_32_K_32",
33:"E_33_D_33",
34:"V_34_V_34",
35:"K_35_K_35",
36:"L_36_L_36",
37:"L_37_L_37",
38:"L_38_L_38",
39:"L_39_L_39",
40:"G_40_G_40",
41:"A_41_A_41",
42:"G_42_G_42",
43:"E_43_E_43",
44:"S_44_S_44",
45:"G_45_G_45",
46:"K_46_K_46",
47:"S_47_S_47",
48:"T_48_T_48",
49:"I_49_I_49",
50:"V_50_V_50",
51:"K_51_K_51",
52:"Q_52_Q_52",
53:"M_53_M_53",
54:"K_54_K_54",
55:"I_55_I_55",
56:"I_56_I_56",
57:"H_57_H_57",
177:"T_177_T_178",
178:"R_178_R_179",
179:"V_179_V_180",
180:"K_180_K_181",
181:"T_181_T_182",
182:"T_182_T_183",
183:"G_183_G_184",
184:"I_184_I_185",
185:"V_185_V_186",
186:"E_186_E_187",
187:"T_187_T_188",
188:"H_188_H_189",
189:"F_189_F_190",
190:"T_190_T_191",
191:"F_191_F_192",
192:"K_192_K_193",
193:"D_193_N_194",
194:"L_194_L_195",
195:"H_195_H_196",
196:"F_196_F_197",
197:"K_197_R_198",
198:"M_198_L_199",
199:"F_199_F_200",
200:"D_200_D_201",
201:"V_201_V_202",
202:"G_202_G_203",
203:"G_203_G_204",
204:"Q_204_Q_205",
205:"R_205_R_206",
206:"S_206_S_207",
207:"E_207_E_208",
208:"R_208_R_209",
209:"K_209_K_210",
210:"K_210_K_211",
211:"W_211_W_212",
212:"I_212_I_213",
213:"H_213_H_214",
214:"C_214_C_215",
215:"F_215_F_216",
216:"E_216_E_217",
217:"G_217_D_218",
218:"V_218_V_219",
219:"T_219_T_220",
220:"A_220_A_221",
221:"I_221_I_222",
222:"I_222_I_223",
223:"F_223_F_224",
224:"C_224_C_225",
225:"V_225_V_226",
226:"A_226_A_227",
227:"L_227_L_228",
228:"S_228_S_229",
229:"D_229_G_230",
230:"Y_230_Y_231",
231:"D_231_D_232",
232:"L_232_Q_233",
233:"V_233_V_234",
234:"L_234_L_235",
235:"A_235_H_236",
236:"E_236_E_237",
237:"D_237_D_238",
238:"E_238_E_239",
239:"E_239_T_240",
240:"M_240_T_241",
241:"N_241_N_242",
242:"R_242_R_243",
243:"M_243_M_244",
244:"H_244_H_245",
245:"E_245_E_246",
246:"S_246_S_247",
247:"M_247_L_248",
248:"K_248_M_249",
249:"L_249_L_250",
250:"F_250_F_251",
251:"D_251_D_252",
252:"S_252_S_253",
253:"I_253_I_254",
254:"C_254_C_255",
255:"N_255_N_256",
256:"N_256_N_257",
257:"K_257_K_258",
258:"W_258_F_259",
259:"F_259_F_260",
260:"T_260_I_261",
261:"D_261_D_262",
262:"T_262_T_263",
263:"S_263_S_264",
264:"I_264_I_265",
265:"I_265_I_266",
266:"L_266_L_267",
267:"F_267_F_268",
268:"L_268_L_269",
269:"N_269_N_270",
270:"K_270_K_271",
271:"K_271_K_272",
272:"D_272_D_273",
273:"L_273_L_274",
274:"F_274_F_275",
275:"E_275_G_276",
276:"E_276_E_277",
277:"K_277_K_278",
278:"I_278_I_279",
279:"K_279_K_280",
280:"K_280_K_281",
281:"S_281_S_282",
282:"P_282_P_283",
283:"L_283_L_284",
284:"T_284_T_285",
285:"I_285_I_286",
286:"C_286_C_287",
287:"Y_287_F_288",
288:"P_288_P_289",
289:"E_289_E_290",
290:"Y_290_Y_291",
291:"A_291_T_292",
292:"G_292_G_293",
293:"S_293_P_294",
294:"N_294_N_295",
295:"T_295_T_296",
296:"Y_296_Y_297",
297:"E_297_E_298",
298:"E_298_D_299",
299:"A_299_A_300",
300:"A_300_A_301",
301:"A_301_A_302",
302:"Y_302_Y_303",
303:"I_303_I_304",
304:"Q_304_Q_305",
305:"C_305_A_306",
306:"Q_306_Q_307",
307:"F_307_F_308",
308:"E_308_E_309",
309:"D_309_S_310",
310:"L_310_K_311",
311:"N_311_N_312",
312:"K_312_x_x",
313:"R_313_R_313",
314:"K_314_S_314",
315:"D_315_P_315",
316:"T_316_N_316",
317:"K_317_K_317",
318:"E_318_E_318",
319:"I_319_I_319",
320:"Y_320_Y_320",
321:"T_321_C_321",
322:"H_322_H_322",
323:"F_323_M_323",
324:"T_324_T_324",
325:"C_325_C_325",
326:"A_326_A_326",
327:"T_327_T_327",
328:"D_328_D_328",
329:"T_329_T_329",
330:"K_330_N_330",
331:"N_331_N_331",
332:"V_332_I_332",
333:"Q_333_Q_333",
334:"F_334_V_334",
335:"V_335_V_335",
336:"F_336_F_336",
337:"D_337_D_337",
338:"A_338_A_338",
339:"V_339_V_339",
340:"T_340_T_340",
341:"D_341_D_341",
342:"V_342_I_342",
343:"I_343_I_343",
344:"I_344_I_344",
345:"K_345_A_345",
346:"N_346_N_346",
347:"N_347_N_347",
348:"L_348_L_348",
349:"K_349_R_349",
350:"D_350_G_350",
351:"C_351_C_351",
352:"G_352_G_352",
353:"L_353_L_353",
354:"F_354_Y_354"
}


go_table={
1:"M_1_M_1",
2:"G_2_G_2",
3:"C_3_C_3",
4:"T_4_T_4",
5:"L_5_L_5",
6:"S_6_S_6",
7:"A_7_A_7",
8:"E_8_E_8",
9:"D_9_E_9",
10:"K_10_R_10",
11:"A_11_A_11",
12:"A_12_A_12",
13:"V_13_L_13",
14:"E_14_E_14",
15:"R_15_R_15",
16:"S_16_S_16",
17:"K_17_K_17",
18:"M_18_A_18",
19:"I_19_I_19",
20:"D_20_E_20",
21:"R_21_K_21",
22:"N_22_N_22",
23:"L_23_L_23",
24:"R_24_K_24",
25:"E_25_E_25",
26:"D_26_D_26",
27:"G_27_G_27",
28:"E_28_I_28",
29:"K_29_S_29",
30:"A_30_A_30",
31:"A_31_A_31",
32:"R_32_K_32",
33:"E_33_D_33",
34:"V_34_V_34",
35:"K_35_K_35",
36:"L_36_L_36",
37:"L_37_L_37",
38:"L_38_L_38",
39:"L_39_L_39",
40:"G_40_G_40",
41:"A_41_A_41",
42:"G_42_G_42",
43:"E_43_E_43",
44:"S_44_S_44",
45:"G_45_G_45",
46:"K_46_K_46",
47:"S_47_S_47",
48:"T_48_T_48",
49:"I_49_I_49",
50:"V_50_V_50",
51:"K_51_K_51",
52:"Q_52_Q_52",
53:"M_53_M_53",
54:"K_54_K_54",
55:"I_55_I_55",
56:"I_56_I_56",
57:"H_57_H_57",
178:"T_177_T_178",
179:"R_178_R_179",
180:"V_179_V_180",
181:"K_180_K_181",
182:"T_181_T_182",
183:"T_182_T_183",
184:"G_183_G_184",
185:"I_184_I_185",
186:"V_185_V_186",
187:"E_186_E_187",
188:"T_187_T_188",
189:"H_188_H_189",
190:"F_189_F_190",
191:"T_190_T_191",
192:"F_191_F_192",
193:"K_192_K_193",
194:"D_193_N_194",
195:"L_194_L_195",
196:"H_195_H_196",
197:"F_196_F_197",
198:"K_197_R_198",
199:"M_198_L_199",
200:"F_199_F_200",
201:"D_200_D_201",
202:"V_201_V_202",
203:"G_202_G_203",
204:"G_203_G_204",
205:"Q_204_Q_205",
206:"R_205_R_206",
207:"S_206_S_207",
208:"E_207_E_208",
209:"R_208_R_209",
210:"K_209_K_210",
211:"K_210_K_211",
212:"W_211_W_212",
213:"I_212_I_213",
214:"H_213_H_214",
215:"C_214_C_215",
216:"F_215_F_216",
217:"E_216_E_217",
218:"G_217_D_218",
219:"V_218_V_219",
220:"T_219_T_220",
221:"A_220_A_221",
222:"I_221_I_222",
223:"I_222_I_223",
224:"F_223_F_224",
225:"C_224_C_225",
226:"V_225_V_226",
227:"A_226_A_227",
228:"L_227_L_228",
229:"S_228_S_229",
230:"D_229_G_230",
231:"Y_230_Y_231",
232:"D_231_D_232",
233:"L_232_Q_233",
234:"V_233_V_234",
235:"L_234_L_235",
236:"A_235_H_236",
237:"E_236_E_237",
238:"D_237_D_238",
239:"E_238_E_239",
240:"E_239_T_240",
241:"M_240_T_241",
242:"N_241_N_242",
243:"R_242_R_243",
244:"M_243_M_244",
245:"H_244_H_245",
246:"E_245_E_246",
247:"S_246_S_247",
248:"M_247_L_248",
249:"K_248_M_249",
250:"L_249_L_250",
251:"F_250_F_251",
252:"D_251_D_252",
253:"S_252_S_253",
254:"I_253_I_254",
255:"C_254_C_255",
256:"N_255_N_256",
257:"N_256_N_257",
258:"K_257_K_258",
259:"W_258_F_259",
260:"F_259_F_260",
261:"T_260_I_261",
262:"D_261_D_262",
263:"T_262_T_263",
264:"S_263_S_264",
265:"I_264_I_265",
266:"I_265_I_266",
267:"L_266_L_267",
268:"F_267_F_268",
269:"L_268_L_269",
270:"N_269_N_270",
271:"K_270_K_271",
272:"K_271_K_272",
273:"D_272_D_273",
274:"L_273_L_274",
275:"F_274_F_275",
276:"E_275_G_276",
277:"E_276_E_277",
278:"K_277_K_278",
279:"I_278_I_279",
280:"K_279_K_280",
281:"K_280_K_281",
282:"S_281_S_282",
283:"P_282_P_283",
284:"L_283_L_284",
285:"T_284_T_285",
286:"I_285_I_286",
287:"C_286_C_287",
288:"Y_287_F_288",
289:"P_288_P_289",
290:"E_289_E_290",
291:"Y_290_Y_291",
292:"A_291_T_292",
293:"G_292_G_293",
294:"S_293_P_294",
295:"N_294_N_295",
296:"T_295_T_296",
297:"Y_296_Y_297",
298:"E_297_E_298",
299:"E_298_D_299",
300:"A_299_A_300",
301:"A_300_A_301",
302:"A_301_A_302",
303:"Y_302_Y_303",
304:"I_303_I_304",
305:"Q_304_Q_305",
306:"C_305_A_306",
307:"Q_306_Q_307",
308:"F_307_F_308",
309:"E_308_E_309",
310:"D_309_S_310",
311:"L_310_K_311",
312:"N_311_N_312",
None:"K_312_x_x",
313:"R_313_R_313",
314:"K_314_S_314",
315:"D_315_P_315",
316:"T_316_N_316",
317:"K_317_K_317",
318:"E_318_E_318",
319:"I_319_I_319",
320:"Y_320_Y_320",
321:"T_321_C_321",
322:"H_322_H_322",
323:"F_323_M_323",
324:"T_324_T_324",
325:"C_325_C_325",
326:"A_326_A_326",
327:"T_327_T_327",
328:"D_328_D_328",
329:"T_329_T_329",
330:"K_330_N_330",
331:"N_331_N_331",
332:"V_332_I_332",
333:"Q_333_Q_333",
334:"F_334_V_334",
335:"V_335_V_335",
336:"F_336_F_336",
337:"D_337_D_337",
338:"A_338_A_338",
339:"V_339_V_339",
340:"T_340_T_340",
341:"D_341_D_341",
342:"V_342_I_342",
343:"I_343_I_343",
344:"I_344_I_344",
345:"K_345_A_345",
346:"N_346_N_346",
347:"N_347_N_347",
348:"L_348_L_348",
349:"K_349_R_349",
350:"D_350_G_350",
351:"C_351_C_351",
352:"G_352_G_352",
353:"L_353_L_353",
354:"F_354_Y_354"
}
# def get_newlabel_bwidx_resn(in_resids,dict_resid_resn,dict_bwidx):
#     labelout=[]
#     _resns=resid_to_resname(dict_resid_resn,in_resids)
#     _bwidx=resid_to_bwidx(dict_bwidx, in_resids)
#     assert(len(_resns)==len(_bwidx))
#     for i in range(len(_resns)):
#         labelout.append(str(_bwidx[i])+"_"+str(_resns[i]))
#     return labelout

# def get_newlabel_resid_resn(in_resids,dict_resid_resn):
#     labelout=[]
#     _resns=resid_to_resname(dict_resid_resn,in_resids)
#     for i in range(len(_resns)):
#         labelout.append(str(in_resids[i])+"_"+str(_resns[i]))
#     return labelout




# def get_matrix_contacts(indf, receptors, resids_A, resids_B, dist_cutoff=8, freq_cutoff=0.6):
#     # select subset of dataframe
#     df_sel = indf[indf['Receptor'] == receptors]
#     df_sel = df_sel.dropna(axis=1, how='all')
#     df_sel = df_sel.reset_index()


# def get_matrix_nrow_ncol(indf, receptors, dist_cutoff=8, freq_cutoff=0.6):
#
#     # select subset of dataframe
#     df_sel = indf[indf['Receptor'] == receptors]
#     df_sel = df_sel.dropna(axis=1, how='all')
#     df_sel = df_sel.reset_index()
#
#     # only need the df.columns start with "cadist_"
#     labels = []
#     for i in list(df_sel.columns):
#         if "cadist_" in i:
#             labels.append(i)
#
#     # replace "cadist_"
#     labels = [item.replace("cadist_", "") for item in labels]
#
#     # parse chain ID
#     sel0_chainid = list(set([item.split("_")[0] for item in labels]))[0]
#     sel1_chainid = list(set([item.split("_")[2] for item in labels]))[0]
#     assert len(sel0_chainid) == len(sel1_chainid)
#
#     # parse residue IDs
#     resids_seli = [item.split("_")[1] for item in labels]
#     resids_seli = natsorted(list(set(resids_seli)))
#     resids_selj = [item.split("_")[3] for item in labels]
#     resids_selj = natsorted(list(set(resids_selj)))
#
#     # convert distance dataframe into contact frequency matrix based on distance cutoff (dist_cutoff)
#     matrix = np.zeros((len(resids_seli), len(resids_selj)))
#     for i in range(len(resids_seli)):
#         for j in range(len(resids_selj)):
#             obs_sel = 'cadist_' + sel0_chainid + '_' + resids_seli[i] + '_' + sel1_chainid + '_' + resids_selj[j]
#             try:
#                 matrix[i][j] = round(sum(df_sel[obs_sel] < dist_cutoff) / len(df_sel[obs_sel]), 2)
#             except:
#                 matrix[i][j] = 0
#
#     # selected index
#     idx_seli = sorted(list(set(np.argwhere(matrix > freq_cutoff)[:, 0])))
#     idx_selj = sorted(list(set(np.argwhere(matrix > freq_cutoff)[:, 1])))
#
#     return matrix, idx_seli, idx_selj, resids_seli, resids_selj















# prep; preparing bwindex table
# d3_bwtable=convert_Residue2TMs("DRD3")
# d2_bwtable=convert_Residue2TMs("DRD2")


# def get_dict_resid_resname(inpdb, segid="PROR"):
#     u = mda.Universe(inpdb)
#     _resids=list(u.select_atoms("not resname NMA and name CA and segid "+segid).residues.resids)
#     _resids=list(set(_resids))
#     _resnames=list(u.select_atoms("not resname NMA and name CA and segid "+segid).residues.resnames)
#     assert len(_resids)==len(_resnames)
#     _dict_resid_resname = dict(zip(_resids, _resnames))
#     return _resids, _resnames, _dict_resid_resname




# d3
# inpdb="/home/khlee/desmond/output/d3gi/pd/d3gi_pd.f1/mda_protein.pdb"
# d3_resids, d3_resnames, d3_dict_resid_resname=get_dict_resid_resname(inpdb, segid="PROR")

# d2
# inpdb="/home/khlee/desmond/output/d2gi/bro/d2gi_bro.f1/mda_protein.pdb"
# d2_resids, d2_resnames, d2_dict_resid_resname=get_dict_resid_resname(inpdb, segid="PROR")

# gi
# inpdb="/home/khlee/desmond/output/d3gi/pd/d3gi_pd.f1/mda_protein.pdb"
# gi_resids, gi_resnames, gi_dict_resid_resname=get_dict_resid_resname(inpdb, segid="PROA")

# go
# inpdb="/home/khlee/desmond/output/d3go/pd/d3go_pd.f1/mda_protein.pdb"
# go_resids, go_resnames, go_dict_resid_resname=get_dict_resid_resname(inpdb, segid="PROA")


# dataframe
# df_d3=pd.DataFrame({'resid': d3_resids,'resn':  d3_resnames,'bwidx': resid_to_bwidx(d3_bwtable, d3_resids)})
# df_d2=pd.DataFrame({'resid': d2_resids,'resn':  d2_resnames,'bwidx': resid_to_bwidx(d2_bwtable, d2_resids)})
# df_gi=pd.DataFrame({'resid': gi_resids,'resn':  gi_resnames})
# df_go=pd.DataFrame({'resid': go_resids,'resn':  go_resnames})



# def get_union_d3(d3gi_idxR, d3go_idxR, d3gi_residsR, d3go_residsR):
#     # union for d3
#     union_d3_idx = list(set(d3gi_idxR + d3go_idxR))  # union of d3gi & d3go index for d3 receptor
#     union_d3_idx.sort()
#     union_d3_resids = []
#     union_d3_resnames = []
#     for i in range(len(union_d3_idx)):
#         assert d3gi_residsR[union_d3_idx[i]] == d3go_residsR[union_d3_idx[i]]
#         _resid = d3gi_residsR[union_d3_idx[i]]
#         _resname = d3_dict_resid_resname[int(_resid)]
#         union_d3_resids.append(_resid)
#         union_d3_resnames.append(_resname)
#     union_d3_bwidx = resid_to_bwidx(d3_bwtable, union_d3_resids)
#     return union_d3_resids, union_d3_resnames, union_d3_bwidx


# def get_union_d2(d2gi_idxR, d2go_idxR, d2gi_residsR, d2go_residsR):
#     # union for d2
#     union_d2_idx = list(set(d2gi_idxR + d2go_idxR))  # union of d3gi & d3go index for d3 receptor
#     union_d2_idx.sort()
#     union_d2_resids = []
#     union_d2_resnames = []
#     for i in range(len(union_d2_idx)):
#         assert d2gi_residsR[union_d2_idx[i]] == d2go_residsR[union_d2_idx[i]]
#         _resid = d2gi_residsR[union_d2_idx[i]]
#         _resname = d2_dict_resid_resname[int(_resid)]
#         union_d2_resids.append(_resid)
#         union_d2_resnames.append(_resname)
#     union_d2_bwidx = resid_to_bwidx(d2_bwtable, union_d2_resids)
#     return union_d2_resids, union_d2_resnames, union_d2_bwidx


# def get_union_gi(d2gi_idxA, d3gi_idxA, d2gi_residsA, d3gi_residsA):
#     # union for gi
#     union_gi_idx = list(set(d2gi_idxA + d3gi_idxA))  # union of d2gi & d3gi index for gi
#     union_gi_idx.sort()
#     union_gi_resids = []
#     union_gi_resnames = []
#     for i in range(len(union_gi_idx)):
#         assert d2gi_residsA[union_gi_idx[i]] == d3gi_residsA[union_gi_idx[i]]
#         _resid = d2gi_residsA[union_gi_idx[i]]
#         _resname = gi_dict_resid_resname[int(_resid)]
#         union_gi_resids.append(_resid)
#         union_gi_resnames.append(_resname)
#     return union_gi_resids, union_gi_resnames


# def get_union_go(d2go_idxA, d3go_idxA, d2go_residsA, d3go_residsA):
#     # union for gi
#     union_go_idx = list(set(d2go_idxA + d3go_idxA))  # union of d2gi & d3gi index for gi
#     union_go_idx.sort()
#     union_go_resids = []
#     union_go_resnames = []
#     for i in range(len(union_go_idx)):
#         assert d2go_residsA[union_go_idx[i]] == d3go_residsA[union_go_idx[i]]
#         _resid = d2go_residsA[union_go_idx[i]]
#         _resname = go_dict_resid_resname[int(_resid)]
#         union_go_resids.append(_resid)
#         union_go_resnames.append(_resname)
#     return union_go_resids, union_go_resnames


# def get_union_gigo(union_gi_resids, union_go_resids):
#     # union_gi_resids
#     # union_go_resids
#
#     union_gi_resids = [int(i) for i in union_gi_resids]
#     union_go_resids = [int(i) for i in union_go_resids]
#
#     # covert residue index from gi to go
#     union_gi_resids_in_go = []
#     for i in union_gi_resids:
#         if (int(i) >= 177 and int(i) <= 311):
#             union_gi_resids_in_go.append(i + 1)
#         else:
#             union_gi_resids_in_go.append(i)
#     union_gigo_resids_in_go = list(set(union_gi_resids_in_go + union_go_resids))
#     union_gigo_resids_in_go.sort()
#
#     # covert residue index go to gi
#     union_gigo_resids_in_gi = []
#     for i in union_gigo_resids_in_go:
#         if (int(i) >= 178 and int(i) <= 312):
#             union_gigo_resids_in_gi.append(i - 1)
#         else:
#             union_gigo_resids_in_gi.append(i)
#
#     union_gigo_resids_in_gi = [str(i) for i in union_gigo_resids_in_gi]
#     union_gigo_resids_in_go = [str(i) for i in union_gigo_resids_in_go]
#
#     return union_gigo_resids_in_gi, union_gigo_resids_in_go


# def get_union_d3_resids_idxR(union_d3_resids, d3_residsR):
#     '''
#     convert union_d3_resids to union_d3_resids_idxR
#     '''
#     union_d3_resids_idxR = []
#     for _resid in union_d3_resids:
#         _idx = d3_residsR.index(_resid)
#         union_d3_resids_idxR.append(_idx)
#     return union_d3_resids_idxR


# def get_union_d2_resids_idxR(union_d2_resids, d2_residsR):
#     '''
#     covert union_d2_resids to union_d2_resids_idxR
#     '''
#     union_d2_resids_idxR = []
#     for _resid in union_d2_resids:
#         _idx = d2_residsR.index(_resid)
#         union_d2_resids_idxR.append(_idx)
#     return union_d2_resids_idxR

#
# def get_union_gigo_resids_in_gi_idxA(union_gigo_resids_in_gi, gi_residsA):
#     union_gigo_resids_in_gi_idxA = []
#     for _resid in union_gigo_resids_in_gi:
#         _idx = gi_residsA.index(_resid)
#         union_gigo_resids_in_gi_idxA.append(_idx)
#     return union_gigo_resids_in_gi_idxA
#
#
# def get_union_gigo_resids_in_go_idxA(union_gigo_resids_in_go, go_residsA):
#     union_gigo_resids_in_go_idxA = []
#     for _resid in union_gigo_resids_in_go:
#         _idx = go_residsA.index(_resid)
#         union_gigo_resids_in_go_idxA.append(_idx)
#     return union_gigo_resids_in_go_idxA
#

# from IPython.core.interactiveshell import InteractiveShell
# InteractiveShell.ast_node_interactivity = 'all'



def convert_rec_resids_to_indexlabels(indata, indict_tm, indict_loop):
    '''
    '''
    _indexlabels = []
    for i in indata:
        try:
            _indexlabels.append(indict_tm[int(i)])
        except:
            _indexlabels.append(indict_loop[int(i)])
    return _indexlabels


def get_key(my_dict, val):
    for key, value in my_dict.items():
        if val == value:
            return key


def convert_rec_indexlabels_to_resids(in_union, in_tmtable, in_lptable):
    ls_resids = []
    for istr in in_union:
        if 'IL' in istr:
            ls_resids.append(get_key(in_lptable, istr))
        elif 'TM' in istr:
            ls_resids.append(get_key(in_tmtable, istr))
    return ls_resids


def convert_ga_resids_to_indexlabels(indata, indict):
    '''given residue ids as the indata, corresponding dictronary and return index labels'''
    _indexlabels = []
    for i in indata:
        _indexlabels.append(indict[int(i)])

    return _indexlabels


def convert_ga_indexlabels_to_resids(in_union, in_dict):
    '''convert index labels back to residue ids'''
    _resids = []
    for istr in in_union:
        for key, value in in_dict.items():
            if value == istr:
                _resids.append(key)
    # _resids.sort()
    return _resids


def get_union(inlistA, inlistB):
    '''given two lists, return the union'''
    ls_union = list(set(inlistA + inlistB))
    ls_union.sort()
    return ls_union


def prep_resids(indf, receptor, method='hhdist'):

    # method
    if method=='hhdist':
        inmethod = "hhdist_"
    elif method=='cadist':
        inmethod = "cadist_"
    else:
        print("method needed. either using cadist or hhhist.")
        return 0

    # select subset of dataframe
    df_sel = indf[indf['Receptor'] == receptor]
    df_sel = df_sel.dropna(axis=1, how='all')
    df_sel = df_sel.reset_index()

    # only need the df.columns start with "cadist_"
    labels = []
    for i in list(df_sel.columns):
        if inmethod in i:
            labels.append(i)

    labels = [item.replace(inmethod, "") for item in labels]

    # parse residue IDs
    resids_seli = [item.split("_")[1] for item in labels]
    resids_seli = natsorted(list(set(resids_seli)))
    resids_selj = [item.split("_")[3] for item in labels]
    resids_selj = natsorted(list(set(resids_selj)))

    return resids_seli, resids_selj

def get_contacts_frequency_matrix(indf, receptor, resids_A, resids_B, dist_cutoff=5, method='hhdist'):

    # method
    if method=='hhdist':
        inmethod = "hhdist_"
    elif method=='cadist':
        inmethod = "cadist_"
    else:
        print("method needed. either using cadist or hhhist.")
        return 0

    # select subset of dataframe
    df_sel = indf[indf['Receptor'] == receptor]
    df_sel = df_sel.reset_index()

    matrix = np.zeros((len(resids_A), len(resids_B)))
    for i in range(len(resids_A)):
        for j in range(len(resids_B)):
            obs_sel = inmethod + 'R' + '_' + str(resids_A[i]) + '_' + 'A' + '_' + str(resids_B[j])
            try:
                matrix[i][j] = round(sum(df_sel[obs_sel] < dist_cutoff) / len(df_sel[obs_sel]), 2)
            except:
                matrix[i][j] = 0

    return matrix


# def reduce_matrix(matrixdiff, xlabels, ylabels, freq_cutoff=0.45):
#     idx_selia = sorted(list(set(np.argwhere(matrixdiff > freq_cutoff)[:, 0])))
#     idx_selja = sorted(list(set(np.argwhere(matrixdiff > freq_cutoff)[:, 1])))
#     idx_selib = sorted(list(set(np.argwhere(matrixdiff < -freq_cutoff)[:, 0])))
#     idx_seljb = sorted(list(set(np.argwhere(matrixdiff < -freq_cutoff)[:, 1])))
#     idx_seli = sorted(list(set(idx_selia + idx_selib)))
#     idx_selj = sorted(list(set(idx_selja + idx_seljb)))
#
#     _xlabels = []
#     [_xlabels.append(xlabels[i]) for i in idx_seli]
#     _ylabels = []
#     [_ylabels.append(ylabels[i]) for i in idx_selj]
#
#     _matrix = matrixdiff[idx_seli, :]
#     _matrix = _matrix[:, idx_selj]
#
#     matrixout = _matrix
#     ylabels = _ylabels
#     xlabels = _xlabels
#
#     plt.figure(figsize=np.shape(matrixout))
#     plot_index = 1
#     ax = plt.subplot(1, 1, plot_index)
#     # im,cbar = heatmap(matrixout, _xlabels, _ylabels, ax=ax, cmap="bwr", cbarlabel="Frequency")
#     im, cbar = heatmap(matrixout, xlabels, ylabels, ax=ax, cmap="bwr", cbarlabel="Frequency", cbar_min=-1, cbar_max=1)
#     texts = annotate_heatmap(im, valfmt="{x:.1f}", size=8)
#     # ax.set_title(title_diff)
#     plt.xlabel(x_axis_label)
#     plt.ylabel(y_axis_label)




# def reduce_matrix(matrixdiff, xlabels, ylabels, x_axis_label="", y_axis_label="", freq_cutoff=0.45, out_png_name="png"):
#     idx_selia = sorted(list(set(np.argwhere(matrixdiff > freq_cutoff)[:, 0])))
#     idx_selja = sorted(list(set(np.argwhere(matrixdiff > freq_cutoff)[:, 1])))
#     idx_selib = sorted(list(set(np.argwhere(matrixdiff < -freq_cutoff)[:, 0])))
#     idx_seljb = sorted(list(set(np.argwhere(matrixdiff < -freq_cutoff)[:, 1])))
#     idx_seli = sorted(list(set(idx_selia + idx_selib)))
#     idx_selj = sorted(list(set(idx_selja + idx_seljb)))
#
#     _xlabels = []
#     [_xlabels.append(xlabels[i]) for i in idx_seli]
#     _ylabels = []
#     [_ylabels.append(ylabels[i]) for i in idx_selj]
#
#     _matrix = matrixdiff[idx_seli, :]
#     _matrix = _matrix[:, idx_selj]
#
#     matrixout = _matrix
#     ylabels = _ylabels
#     xlabels = _xlabels
#
#     plt.figure(figsize=np.shape(matrixout))
#     plot_index = 1
#     ax = plt.subplot(1, 1, plot_index)
#     # im,cbar = heatmap(matrixout, _xlabels, _ylabels, ax=ax, cmap="bwr", cbarlabel="Frequency")
#     im, cbar = heatmap(matrixout, xlabels, ylabels, ax=ax, cmap="bwr", cbarlabel="Frequency", cbar_min=-1, cbar_max=1)
#     texts = annotate_heatmap(im, valfmt="{x:.1f}", size=8)
#     # ax.set_title(title_diff)
#     plt.xlabel(x_axis_label)
#     plt.ylabel(y_axis_label)
#
#     plt.savefig(out_png_name, dpi=300, bbox_inches='tight')
#
#







def reduce_matrix(matrixdiff, xlabels, ylabels, freq_cutoff=0.45, x_axis_label="", y_axis_label="", out_png_name=None):
    idx_selia = sorted(list(set(np.argwhere(matrixdiff > freq_cutoff)[:,0])))
    idx_selja = sorted(list(set(np.argwhere(matrixdiff > freq_cutoff)[:,1])))
    idx_selib = sorted(list(set(np.argwhere(matrixdiff < -freq_cutoff)[:,0])))
    idx_seljb = sorted(list(set(np.argwhere(matrixdiff < -freq_cutoff)[:,1])))
    idx_seli = sorted(list(set(idx_selia + idx_selib)))
    idx_selj = sorted(list(set(idx_selja + idx_seljb)))

    _xlabels=[]
    [_xlabels.append(xlabels[i]) for i in idx_seli]
    _ylabels=[]
    [_ylabels.append(ylabels[i]) for i in idx_selj]


    _matrix = matrixdiff[idx_seli, :]
    _matrix = _matrix[:, idx_selj]

    matrixout=_matrix
    ylabels=_ylabels
    xlabels=_xlabels

    #plt.figure(figsize=np.shape(matrixout))
    # plt.figure(figsize=(2*len(_xlabels),2*len(_ylabels)))
    # plot_index = 1
    # ax = plt.subplot(1, 1, plot_index)
    # #im,cbar = heatmap(matrixout, _xlabels, _ylabels, ax=ax, cmap="bwr", cbarlabel="Frequency")
    # im,cbar = heatmap(matrixout,xlabels, ylabels,ax=ax,cmap="bwr", cbarlabel="Frequency",cbar_min=-1, cbar_max=1)
    # texts = annotate_heatmap(im, valfmt="{x:.1f}", size=8)
    # #ax.set_title(title_diff)
    # plt.xlabel(x_axis_label)
    # plt.ylabel(y_axis_label)
    # if out_png_name!=None:
    #     print("save to png")
    #     plt.savefig(out_png_name, dpi=300, bbox_inches='tight')
    return matrixout, xlabels, ylabels

def reduce_matrix_delpolyG(matrixdiff, xlabels, ylabels, freq_cutoff=0.45, x_axis_label="", y_axis_label=""):
    idx_selia = sorted(list(set(np.argwhere(matrixdiff > freq_cutoff)[:,0])))
    idx_selja = sorted(list(set(np.argwhere(matrixdiff > freq_cutoff)[:,1])))
    idx_selib = sorted(list(set(np.argwhere(matrixdiff < -freq_cutoff)[:,0])))
    idx_seljb = sorted(list(set(np.argwhere(matrixdiff < -freq_cutoff)[:,1])))
    #idx_seli = sorted(list(set(idx_selia + idx_selib)))
    idx_selj = sorted(list(set(idx_selja + idx_seljb)))

    ls_delployG = []
    for i in range(len(xlabels)):
        if xlabels[i].split("_")[-1] == "GLY":
            ls_delployG.append(i)
    idx_seli = sorted(list(set(idx_selia+idx_selib) - set(ls_delployG)))

    _xlabels=[]
    [_xlabels.append(xlabels[i]) for i in idx_seli]
    _ylabels=[]
    [_ylabels.append(ylabels[i]) for i in idx_selj]


    _matrix = matrixdiff[idx_seli, :]
    _matrix = _matrix[:, idx_selj]

    matrixout=_matrix
    ylabels=_ylabels
    xlabels=_xlabels

#     #plt.figure(figsize=np.shape(matrixout))
#     plt.figure(figsize=(2*len(_xlabels),2*len(_ylabels)))
#     plot_index = 1
#     ax = plt.subplot(1, 1, plot_index)
#     #im,cbar = heatmap(matrixout, _xlabels, _ylabels, ax=ax, cmap="bwr", cbarlabel="Frequency")
#     im,cbar = heatmap(matrixout,xlabels, ylabels,ax=ax,cmap="bwr", cbarlabel="Frequency",cbar_min=-1, cbar_max=1)
#     texts = annotate_heatmap(im, valfmt="{x:.1f}", size=8)
#     #ax.set_title(title_diff)
#     plt.xlabel(x_axis_label)
#     plt.ylabel(y_axis_label)
    return matrixout, xlabels, ylabels



def get_rec_union_index(matrixes, freq_cutoff, union_rec, noployG=True):
    '''
    given matrixes, freq_cutoff, union_rec, noployG flag
    return union indexes for i & j
    '''
    _seli = []
    _selj = []
    for imatrix in matrixes:
        idx_selia = sorted(list(set(np.argwhere(imatrix > freq_cutoff)[:,0])))
        idx_selja = sorted(list(set(np.argwhere(imatrix > freq_cutoff)[:,1])))
        idx_selib = sorted(list(set(np.argwhere(imatrix < -freq_cutoff)[:,0])))
        idx_seljb = sorted(list(set(np.argwhere(imatrix < -freq_cutoff)[:,1])))
        if noployG:
            ls_delployG = []
            for i in range(len(union_rec)):
                if union_rec[i].split("_")[-1] == "GLY":
                    ls_delployG.append(i)
            idx_seli = sorted(list(set(idx_selia+idx_selib) - set(ls_delployG)))
        else:
            idx_seli = sorted(list(set(idx_selia + idx_selib)))
        idx_selj = sorted(list(set(idx_selja + idx_seljb)))
        _seli = _seli + idx_seli
        _selj = _selj + idx_selj
    xxseli = sorted(list(set(_seli)))
    xxselj = sorted(list(set(_selj)))

    return xxseli, xxselj



def get_union_matrix(matrixin, xlabels, ylabels, seli, selj):
    '''get union for all 4 cases'''
    _xlabels=[]
    [_xlabels.append(xlabels[i]) for i in seli]
    _ylabels=[]
    [_ylabels.append(ylabels[i]) for i in selj]

    _matrix = matrixin[seli, :]
    _matrix = _matrix[:, selj]

    matrixout=_matrix
    ylabels=_ylabels
    xlabels=_xlabels
    return matrixout, xlabels, ylabels


def get_counts(matrixin):
    ''' get_counts
    given a matrix (matrixin) and counting the elements with frequency > 0.9, 0.8, 0.7, 0.6 and 0.5'''
    ls_counts = []
    for i in range(5):
        try:
            ls_counts.append(sum(matrixin.flatten(order='C') >= 1-0.1*i))
        except:
            ls_counts.append(0)
    return ls_counts

def get_counts2(matrixin):
    ''' get_counts
    given a matrix (matrixin) and counting the elements with frequency > 0.9, 0.8, 0.7, 0.6 and 0.5
                                                        with frequency < -0.9, -0.8, -0.7, -0.6 and -0.5
    '''
    ls_counts = []
    ls_counts.append(sum(matrixin.flatten(order='C') >= 0.9))
    ls_counts.append(sum(matrixin.flatten(order='C') >= 0.8))
    ls_counts.append(sum(matrixin.flatten(order='C') >= 0.7))
    ls_counts.append(sum(matrixin.flatten(order='C') >= 0.6))
    ls_counts.append(sum(matrixin.flatten(order='C') >= 0.5))
    ls_counts.append(sum(matrixin.flatten(order='C') <= -0.5))
    ls_counts.append(sum(matrixin.flatten(order='C') <= -0.6))
    ls_counts.append(sum(matrixin.flatten(order='C') <= -0.7))
    ls_counts.append(sum(matrixin.flatten(order='C') <= -0.8))
    ls_counts.append(sum(matrixin.flatten(order='C') <= -0.9))
    return ls_counts



def combine_labels(ylabels):
    ylabels2=[]
    for ilabel in ylabels:
        ylabels2.append(ilabel.split("_")[0]+ilabel.split("_")[1]+"_"+ilabel.split("_")[2]+ilabel.split("_")[3])
    return ylabels2

def combine_labels_rec(reclabels):
    ylabels2=[]
    for ilabel in reclabels:
        ylabels2.append(ilabel.split("_")[0]+"_"+ilabel.split("_")[1]+ilabel.split("_")[2]+"_"+ilabel.split("_")[3]+ilabel.split("_")[4])
    return ylabels2


# def color_xlables(ylabels):
#     for ii in range(len(ylabels)):
#         ilabel = ylabels[ii]
#         #print(ilabel.split("_")[0],ilabel.split("_")[2])
#         if ilabel.split("_")[0] != ilabel.split("_")[2]:
#             ax.get_xticklabels()[ii].set_color("red")