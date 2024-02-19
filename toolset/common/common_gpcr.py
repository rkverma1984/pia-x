from vmd import atomsel
from vmd import molecule
import pandas as pd
import os
from pathlib import Path

import sys, os
from pathlib import Path
home = str(Path.home())
toolset_dir = home+'/repositories/pia-x/toolset'
toolset_common = toolset_dir+"/common"
toolset_gpcr = toolset_dir+"/gpcr"
sys.path.insert(0, toolset_dir)
sys.path.insert(0, toolset_common)
sys.path.insert(0, toolset_gpcr)
from common_vmd import *




def prep_vmd_atomsel(vmd_prot_R, vmd_prot_A, vmd_prot_B, vmd_prot_C,
                     vmd_ignore_prot_A_sel,
                     vmd_add_backbone, vmd_add_heavyatoms, vmd_add_CA,
                     vmd_add_rec_tm, OBS_resid, vmd_ligand, nframes):


    # parsing atomselection string
    _prot_R_tm = vmd_prot_R+vmd_add_rec_tm
    _prot_A_sel = vmd_prot_A+vmd_ignore_prot_A_sel
    _prot_B = vmd_prot_B
    _prot_C = vmd_prot_C
    _prot_sel = "("+_prot_R_tm+" or "+_prot_A_sel+" or "+_prot_B+" or "+_prot_C+")"
    _all = vmd_all
    _ligand = vmd_ligand
    _OBS = vmd_prot_R+" and resid "+str(OBS_resid).replace(",","").replace("]","").replace("[","")
    _backbone = vmd_add_backbone
    _heavyatoms = vmd_add_heavyatoms
    _CA = vmd_add_CA

    sel_allatom = atomsel(_all)
    sel_backbone = atomsel(_all+_backbone)
    ref_allatom = atomsel(_all, frame=0)
    ref_backbone = atomsel(_all+_backbone, frame=0)

    sel_prot_sel = atomsel(_prot_sel)
    sel_prot_sel_backbone = atomsel(_prot_sel+_backbone)
    ref_prot_sel=atomsel(_prot_sel, frame=0)
    ref_prot_sel_backbone=atomsel(_prot_sel+vmd_add_backbone, frame=0)

    sel_prot_R_tm = atomsel(_prot_R_tm)
    sel_prot_R_tm_backbone = atomsel(_prot_R_tm+_backbone)
    ref_prot_R_tm = atomsel(_prot_R_tm, frame=0)
    ref_prot_R_tm_backbone = atomsel(_prot_R_tm+vmd_add_backbone, frame=0)

    sel_prot_A_sel = atomsel(_prot_A_sel)
    sel_prot_A_sel_backbone = atomsel(_prot_A_sel+_backbone)
    ref_prot_A_sel =atomsel(_prot_A_sel, frame=0)
    ref_prot_A_sel_backbone=atomsel(_prot_A_sel+vmd_add_backbone, frame=0)

    sel_OBS = atomsel(_OBS)
    sel_OBS_backbone = atomsel(_OBS+_backbone)
    ref_OBS = atomsel(_OBS, frame=0)
    ref_OBS_backbone = atomsel(_OBS+_backbone, frame=0)

    sel_ligand = atomsel(_ligand)
    ref_ligand = atomsel(_ligand, frame=0)
    ref_ligand_last = atomsel(_ligand, frame=(nframes-1))

    return sel_allatom, sel_backbone, ref_allatom, ref_backbone, \
           sel_prot_sel, sel_prot_sel_backbone, ref_prot_sel, ref_prot_sel_backbone,\
           sel_prot_R_tm, sel_prot_R_tm_backbone, ref_prot_R_tm, ref_prot_R_tm_backbone, \
           sel_prot_A_sel, sel_prot_A_sel_backbone, ref_prot_A_sel, ref_prot_A_sel_backbone,\
           sel_OBS, sel_OBS_backbone, ref_OBS, ref_OBS_backbone,\
           sel_ligand, ref_ligand, ref_ligand_last



def get_convert(resid_str=31,resid_end=62,BW_end=60,TMx='TM1'):
    GPCR_BW_dic = {}
    for i in range(resid_end, resid_str-1, -1): #TM1
        GPCR_BW_dic[i]=TMx+"."+str(BW_end)
        BW_end=BW_end-1
    return GPCR_BW_dic


def convert_Residue2TMs(gpcrin):
    BW_dic = {}
    if gpcrin=="DRD3":
        BW_dic.update(get_convert(17,57,60,'TM1'))
        BW_dic.update(get_convert(63,91,66,'TM2'))
        BW_dic.update(get_convert(100,134,56,'TM3'))
        BW_dic.update(get_convert(145,170,62,'TM4'))
        BW_dic.update(get_convert(186,218,68,'TM5'))
        BW_dic.update(get_convert(322,355,61,'TM6'))
        BW_dic.update(get_convert(361,400,70,'TM7'))
        BW_dic.update({96: 'EL1.50'})
        BW_dic.update({181: 'EL2.50'})
        BW_dic.update({182: 'EL2.51'})
        BW_dic.update({183: 'EL2.52'})
        BW_dic.update({184: 'EL2.53'})
        BW_dic.update({185: 'EL2.54'})
        # beginnings
        assert BW_dic[17] == 'TM1.20'
        assert BW_dic[63] == 'TM2.38'
        assert BW_dic[100] == 'TM3.22'
        assert BW_dic[145] == 'TM4.37'
        assert BW_dic[186] == 'TM5.36'
        assert BW_dic[322] == 'TM6.28'
        assert BW_dic[361] == 'TM7.31'
        # endings
        assert BW_dic[57] == 'TM1.60'
        assert BW_dic[91] == 'TM2.66'
        assert BW_dic[134] == 'TM3.56'
        assert BW_dic[170] == 'TM4.62'
        assert BW_dic[218] == 'TM5.68'
        assert BW_dic[354] == 'TM6.60'
        assert BW_dic[362] == 'TM7.32'
    elif gpcrin == "DRD2":
        BW_dic.update(get_convert(21, 62, 60, 'TM1'))
        BW_dic.update(get_convert(68, 96, 66, 'TM2'))
        BW_dic.update(get_convert(104, 138, 56, 'TM3'))
        BW_dic.update(get_convert(147, 172, 62, 'TM4'))
        BW_dic.update(get_convert(187, 219, 68, 'TM5'))
        BW_dic.update(get_convert(366, 399, 61, 'TM6'))
        BW_dic.update(get_convert(404, 443, 70, 'TM7'))
        BW_dic.update({100: 'EL1.50'})
        BW_dic.update({182: 'EL2.50'})
        BW_dic.update({183: 'EL2.51'})
        BW_dic.update({184: 'EL2.52'})
        BW_dic.update({185: 'EL2.53'})
        BW_dic.update({186: 'EL2.54'})
        # beginnings
        assert BW_dic[21] == 'TM1.19'
        assert BW_dic[68] == 'TM2.38'
        assert BW_dic[104] == 'TM3.22'
        assert BW_dic[147] == 'TM4.37'
        assert BW_dic[187] == 'TM5.36'
        assert BW_dic[366] == 'TM6.28'
        assert BW_dic[404] == 'TM7.31'
        # endings
        assert BW_dic[62] == 'TM1.60'
        assert BW_dic[96] == 'TM2.66'
        assert BW_dic[138] == 'TM3.56'
        assert BW_dic[172] == 'TM4.62'
        assert BW_dic[219] == 'TM5.68'
        assert BW_dic[398] == 'TM6.60'
        assert BW_dic[405] == 'TM7.32'
    elif gpcrin == "DRD2_short":
        BW_dic.update(get_convert(21, 62, 60, 'TM1'))
        BW_dic.update(get_convert(68, 96, 66, 'TM2'))
        BW_dic.update(get_convert(104, 138, 56, 'TM3'))
        BW_dic.update(get_convert(147, 172, 62, 'TM4'))
        BW_dic.update(get_convert(187, 219, 68, 'TM5'))
        BW_dic.update(get_convert(337, 370, 61, 'TM6'))
        BW_dic.update(get_convert(375, 414, 70, 'TM7'))
        BW_dic.update({100: 'EL1.50'})
        BW_dic.update({182: 'EL2.50'})
        BW_dic.update({183: 'EL2.51'})
        BW_dic.update({184: 'EL2.52'})
        BW_dic.update({185: 'EL2.53'})
        BW_dic.update({186: 'EL2.54'})
        # beginnings
        assert BW_dic[21] == 'TM1.19'
        assert BW_dic[68] == 'TM2.38'
        assert BW_dic[104] == 'TM3.22'
        assert BW_dic[147] == 'TM4.37'
        assert BW_dic[187] == 'TM5.36'
        assert BW_dic[337] == 'TM6.28'
        assert BW_dic[375] == 'TM7.31'
        # endings
        assert BW_dic[62] == 'TM1.60'
        assert BW_dic[96] == 'TM2.66'
        assert BW_dic[138] == 'TM3.56'
        assert BW_dic[172] == 'TM4.62'
        assert BW_dic[219] == 'TM5.68'
        assert BW_dic[369] == 'TM6.60'
        assert BW_dic[376] == 'TM7.32'
    elif gpcrin == "agrl3_mouse":
        # skip this one
        BW_dic.update(get_convert(940, 970, 64, 'TM1'))
        BW_dic.update(get_convert(975, 1001, 70, 'TM2'))
        BW_dic.update(get_convert(1005, 1038, 59, 'TM3'))
        BW_dic.update(get_convert(1046, 1069, 63, 'TM4'))
        BW_dic.update(get_convert(1089, 1120, 67, 'TM5'))
        BW_dic.update(get_convert(1139, 1163, 58, 'TM6'))
        BW_dic.update(get_convert(1174, 1197, 63, 'TM7'))
        assert BW_dic[956] == 'TM1.50'
        assert BW_dic[981] == 'TM2.50'
        assert BW_dic[1029] == 'TM3.50'
        assert BW_dic[1056] == 'TM4.50'
        assert BW_dic[1103] == 'TM5.50'
        assert BW_dic[1155] == 'TM6.50'
        assert BW_dic[1184] == 'TM7.50'
    elif gpcrin == "agrg3_human":
        # skip this one
        BW_dic.update(get_convert(263, 293, 64, 'TM1'))
        BW_dic.update(get_convert(303, 329, 70, 'TM2'))
        BW_dic.update(get_convert(335, 368, 59, 'TM3'))
        BW_dic.update(get_convert(377, 400, 63, 'TM4'))
        BW_dic.update(get_convert(428, 459, 67, 'TM5'))
        BW_dic.update(get_convert(471, 495, 58, 'TM6'))
        BW_dic.update(get_convert(504, 527, 63, 'TM7'))
        assert BW_dic[279] == 'TM1.50'
        assert BW_dic[309] == 'TM2.50'
        assert BW_dic[359] == 'TM3.50'
        assert BW_dic[387] == 'TM4.50'
        assert BW_dic[442] == 'TM5.50'
        assert BW_dic[487] == 'TM6.50'
        assert BW_dic[514] == 'TM7.50'
    elif gpcrin == "DRD1":
        BW_dic.update(get_convert(11,51,60,'TM1'))
        BW_dic.update(get_convert(58,86,66,'TM2'))
        BW_dic.update(get_convert(93,127,56,'TM3'))
        BW_dic.update(get_convert(136,161,63,'TM4'))
        BW_dic.update(get_convert(192,224,68,'TM5'))
        BW_dic.update(get_convert(265,298,61,'TM6'))
        # BW_dic.update(get_convert(310,332,54,'TM7'))
        BW_dic.update(get_convert(310,347,69,'TM7'))
        BW_dic.update({90: 'EL1.50'})
        BW_dic.update({186: 'EL2.50'})
        BW_dic.update({187: 'EL2.51'})
        BW_dic.update({188: 'EL2.52'})
        BW_dic.update({189: 'EL2.53'})
        BW_dic.update({190: 'EL2.54'})
        BW_dic.update({191: 'EL2.55'})
        # beginnings
        assert BW_dic[11] == 'TM1.20'
        assert BW_dic[58] == 'TM2.38'
        assert BW_dic[93] == 'TM3.22'
        assert BW_dic[136] == 'TM4.38'
        assert BW_dic[192] == 'TM5.36'
        assert BW_dic[265] == 'TM6.28'
        assert BW_dic[310] == 'TM7.32'
        # endings
        assert BW_dic[51] == 'TM1.60'
        assert BW_dic[86] == 'TM2.66'
        assert BW_dic[127] == 'TM3.56'
        assert BW_dic[161] == 'TM4.63'
        assert BW_dic[224] == 'TM5.68'
        assert BW_dic[297] == 'TM6.60'
        assert BW_dic[347] == 'TM7.69'
    elif gpcrin == "DRD4":
        BW_dic.update(get_convert(22,62,60,'TM1'))
        BW_dic.update(get_convert(68,96,66,'TM2'))
        BW_dic.update(get_convert(105,139,56,'TM3'))
        BW_dic.update(get_convert(147,172,62,'TM4'))
        BW_dic.update(get_convert(190,222,68,'TM5'))
        BW_dic.update(get_convert(387,420,61,'TM6'))
        # BW_dic.update(get_convert(427,449,54,'TM7'))
        BW_dic.update(get_convert(427,464,69,'TM7'))
        BW_dic.update({101: 'EL1.50'})
        BW_dic.update({185: 'EL2.50'})
        BW_dic.update({186: 'EL2.51'})
        BW_dic.update({187: 'EL2.52'})
        BW_dic.update({188: 'EL2.53'})
        BW_dic.update({189: 'EL2.54'})
        # beginnings
        assert BW_dic[22] == 'TM1.20'
        assert BW_dic[68] == 'TM2.38'
        assert BW_dic[105] == 'TM3.22'
        assert BW_dic[147] == 'TM4.37'
        assert BW_dic[190] == 'TM5.36'
        assert BW_dic[387] == 'TM6.28'
        assert BW_dic[427] == 'TM7.32'
        # endings
        assert BW_dic[62] == 'TM1.60'
        assert BW_dic[96] == 'TM2.66'
        assert BW_dic[139] == 'TM3.56'
        assert BW_dic[172] == 'TM4.62'
        assert BW_dic[222] == 'TM5.68'
        assert BW_dic[419] == 'TM6.60'
        assert BW_dic[464] == 'TM7.69'
    elif gpcrin == "DRD4_mouse":
        BW_dic.update(get_convert(19, 59, 60, 'TM1'))
        BW_dic.update(get_convert(65, 93, 66, 'TM2'))
        BW_dic.update(get_convert(102, 136, 56, 'TM3'))
        BW_dic.update(get_convert(144, 167, 62, 'TM4'))
        BW_dic.update(get_convert(185, 217, 68, 'TM5'))
        BW_dic.update(get_convert(307, 339, 60, 'TM6'))
        BW_dic.update(get_convert(347, 384, 69, 'TM7'))
        BW_dic.update({98: 'EL1.50'})
        BW_dic.update({180: 'EL2.50'})
        BW_dic.update({181: 'EL2.51'})
        BW_dic.update({182: 'EL2.52'})
        BW_dic.update({183: 'EL2.53'})
        BW_dic.update({184: 'EL2.54'})
        # beginnings
        assert BW_dic[19] == 'TM1.20'
        assert BW_dic[65] == 'TM2.38'
        assert BW_dic[102] == 'TM3.22'
        assert BW_dic[144] == 'TM4.39'
        assert BW_dic[185] == 'TM5.36'
        assert BW_dic[307] == 'TM6.28'
        assert BW_dic[347] == 'TM7.32'
        # endings
        assert BW_dic[59] == 'TM1.60'
        assert BW_dic[93] == 'TM2.66'
        assert BW_dic[136] == 'TM3.56'
        assert BW_dic[167] == 'TM4.62'
        assert BW_dic[217] == 'TM5.68'
        assert BW_dic[339] == 'TM6.60'
        assert BW_dic[384] == 'TM7.69'
    elif gpcrin == "B2" or gpcrin == "adrb2" or gpcrin == "ADRB2":
        BW_dic.update(get_convert(21, 61, 60, 'TM1'))
        BW_dic.update(get_convert(67, 95, 66, 'TM2'))
        BW_dic.update(get_convert(103, 137, 56, 'TM3'))
        BW_dic.update(get_convert(146, 171, 63, 'TM4'))
        BW_dic.update(get_convert(197, 229, 68, 'TM5'))
        BW_dic.update(get_convert(266, 299, 61, 'TM6'))
        # BW_dic.update(get_convert(305, 327, 54, 'TM7'))
        BW_dic.update(get_convert(305, 342, 69, 'TM7'))
        BW_dic.update({99: 'EL1.50'})
        BW_dic.update({191: 'EL2.50'})
        BW_dic.update({192: 'EL2.51'})
        BW_dic.update({193: 'EL2.52'})
        BW_dic.update({194: 'EL2.53'})
        BW_dic.update({195: 'EL2.54'})
        BW_dic.update({196: 'EL2.55'})
        # beginnings
        assert BW_dic[21] == 'TM1.20'
        assert BW_dic[67] == 'TM2.38'
        assert BW_dic[103] == 'TM3.22'
        assert BW_dic[146] == 'TM4.38'
        assert BW_dic[197] == 'TM5.36'
        assert BW_dic[266] == 'TM6.28'
        assert BW_dic[305] == 'TM7.32'
        # endings
        assert BW_dic[61] == 'TM1.60'
        assert BW_dic[95] == 'TM2.66'
        assert BW_dic[137] == 'TM3.56'
        assert BW_dic[171] == 'TM4.63'
        assert BW_dic[229] == 'TM5.68'
        assert BW_dic[298] == 'TM6.60'
        assert BW_dic[342] == 'TM7.69'
    elif gpcrin == "B1" or gpcrin == "adrb1" or gpcrin == "adrb1_human" or gpcrin == "ADRB1":
        BW_dic.update(get_convert(46, 86, 60, 'TM1'))  # TM1_start -10
        BW_dic.update(get_convert(92, 120, 66, 'TM2'))
        BW_dic.update(get_convert(128, 162, 56, 'TM3'))
        BW_dic.update(get_convert(171, 196, 63, 'TM4'))
        BW_dic.update(get_convert(222, 254, 68, 'TM5'))
        BW_dic.update(get_convert(317, 350, 61, 'TM6'))  # TM6_end +1
        # BW_dic.update(get_convert(356, 378, 54, 'TM7')) # TM7 & extend to H8
        BW_dic.update(get_convert(356, 393, 69, 'TM7'))
        BW_dic.update({124: 'EL1.50'})
        BW_dic.update({216: 'EL2.50'})
        BW_dic.update({217: 'EL2.51'})
        BW_dic.update({218: 'EL2.52'})
        BW_dic.update({219: 'EL2.53'})
        BW_dic.update({220: 'EL2.54'})
        BW_dic.update({221: 'EL2.55'})
        # beginnings
        assert BW_dic[46] == 'TM1.20'
        assert BW_dic[92] == 'TM2.38'
        assert BW_dic[128] == 'TM3.22'
        assert BW_dic[171] == 'TM4.38'
        assert BW_dic[222] == 'TM5.36'
        assert BW_dic[317] == 'TM6.28'
        assert BW_dic[356] == 'TM7.32'
        # endings
        assert BW_dic[86] == 'TM1.60'
        assert BW_dic[120] == 'TM2.66'
        assert BW_dic[162] == 'TM3.56'
        assert BW_dic[195] == 'TM4.62'
        assert BW_dic[254] == 'TM5.68'
        assert BW_dic[349] == 'TM6.60'
        assert BW_dic[393] == 'TM7.69'
    elif gpcrin == "adrb1_melga" or gpcrin == "mADRB1":
        BW_dic.update(get_convert(29, 69, 60, 'TM1'))  # TM1_start -10
        BW_dic.update(get_convert(75, 103, 66, 'TM2'))
        BW_dic.update(get_convert(111, 145, 56, 'TM3'))
        BW_dic.update(get_convert(154, 179, 63, 'TM4'))
        BW_dic.update(get_convert(205, 237, 68, 'TM5'))
        BW_dic.update(get_convert(283, 316, 61, 'TM6'))  # TM6_end +1
        # BW_dic.update(get_convert(322, 344, 54, 'TM7')) # TM7 & extend to H8
        BW_dic.update(get_convert(322, 359, 69, 'TM7'))
        BW_dic.update({107: 'EL1.50'})
        BW_dic.update({199: 'EL2.50'})
        BW_dic.update({200: 'EL2.51'})
        BW_dic.update({201: 'EL2.52'})
        BW_dic.update({202: 'EL2.53'})
        BW_dic.update({203: 'EL2.54'})
        BW_dic.update({204: 'EL2.55'})
        # beginnings
        assert BW_dic[29] == 'TM1.20'
        assert BW_dic[75] == 'TM2.38'
        assert BW_dic[111] == 'TM3.22'
        assert BW_dic[154] == 'TM4.38'
        assert BW_dic[205] == 'TM5.36'
        assert BW_dic[283] == 'TM6.28'
        assert BW_dic[322] == 'TM7.32'
        # endings
        assert BW_dic[69] == 'TM1.60'
        assert BW_dic[103] == 'TM2.66'
        assert BW_dic[145] == 'TM3.56'
        assert BW_dic[178] == 'TM4.62'
        assert BW_dic[237] == 'TM5.68'
        assert BW_dic[315] == 'TM6.60'
        assert BW_dic[359] == 'TM7.69'
    elif gpcrin == "M2" or gpcrin == "ACM2":
        BW_dic.update(get_convert(11, 51, 60, 'TM1'))
        BW_dic.update(get_convert(57, 85, 66, 'TM2'))
        BW_dic.update(get_convert(93, 127, 56, 'TM3'))
        BW_dic.update(get_convert(136, 161, 63, 'TM4'))
        BW_dic.update(get_convert(184, 216, 68, 'TM5'))
        BW_dic.update(get_convert(380, 413, 61, 'TM6'))
        # BW_dic.update(get_convert(419, 441, 54, 'TM7'))
        BW_dic.update(get_convert(419, 457, 70, 'TM7'))
        BW_dic.update({89: 'EL1.50'})
        BW_dic.update({179: 'EL2.50'}) # this is hack EL2.53 to EL2.50
        BW_dic.update({180: 'EL2.51'}) # this is hack EL2.54 to EL2.51
        BW_dic.update({181: 'EL2.52'}) # this is hack EL2.55 to EL2.52
        BW_dic.update({182: 'EL2.53'}) # this is hack EL2.56 to EL2.53
        BW_dic.update({183: 'EL2.54'}) # this is hack EL2.57 to EL2.54
        # beginnings
        assert BW_dic[11] == 'TM1.20'
        assert BW_dic[57] == 'TM2.38'
        assert BW_dic[93] == 'TM3.22'
        assert BW_dic[136] == 'TM4.38'
        assert BW_dic[184] == 'TM5.36'
        assert BW_dic[380] == 'TM6.28'
        assert BW_dic[419] == 'TM7.32'
        # endings
        assert BW_dic[51] == 'TM1.60'
        assert BW_dic[85] == 'TM2.66'
        assert BW_dic[127] == 'TM3.56'
        assert BW_dic[161] == 'TM4.63'
        assert BW_dic[216] == 'TM5.68'
        assert BW_dic[412] == 'TM6.60'
        assert BW_dic[457] == 'TM7.70'
    elif gpcrin == "ACM1":
        BW_dic.update(get_convert(23, 53, 60, 'TM1'))
        BW_dic.update(get_convert(59, 87, 66, 'TM2'))
        BW_dic.update(get_convert(95, 129, 56, 'TM3'))
        BW_dic.update(get_convert(138, 163, 63, 'TM4'))
        BW_dic.update(get_convert(186, 218, 68, 'TM5'))
        BW_dic.update(get_convert(358, 390, 60, 'TM6'))
        BW_dic.update(get_convert(397, 436, 71, 'TM7'))
        BW_dic.update({91: 'EL1.50'})
        BW_dic.update({178: 'EL2.50'})
        BW_dic.update({179: 'EL2.51'})
        BW_dic.update({180: 'EL2.52'})
        BW_dic.update({181: 'EL2.53'})
        BW_dic.update({182: 'EL2.54'})
        BW_dic.update({183: 'EL2.55'})
        BW_dic.update({184: 'EL2.56'})
        BW_dic.update({185: 'EL2.57'})
        # beginnings
        assert BW_dic[23] == 'TM1.30'
        assert BW_dic[59] == 'TM2.38'
        assert BW_dic[95] == 'TM3.22'
        assert BW_dic[138] == 'TM4.38'
        assert BW_dic[186] == 'TM5.36'
        assert BW_dic[358] == 'TM6.28'
        assert BW_dic[397] == 'TM7.32'
        # endings
        assert BW_dic[53] == 'TM1.60'
        assert BW_dic[87] == 'TM2.66'
        assert BW_dic[129] == 'TM3.56'
        assert BW_dic[163] == 'TM4.63'
        assert BW_dic[218] == 'TM5.68'
        assert BW_dic[390] == 'TM6.60'
        assert BW_dic[436] == 'TM7.71'
    elif gpcrin == "5HT2A":
        BW_dic.update(get_convert(62, 102, 60, 'TM1'))
        BW_dic.update(get_convert(108, 136, 66, 'TM2'))
        BW_dic.update(get_convert(145, 179, 56, 'TM3'))
        BW_dic.update(get_convert(188, 212, 63, 'TM4'))
        BW_dic.update(get_convert(232, 264, 68, 'TM5'))
        BW_dic.update(get_convert(316, 348, 60, 'TM6'))
        BW_dic.update(get_convert(359, 398, 71, 'TM7'))
        BW_dic.update({141: 'EL1.50'})
        BW_dic.update({227: 'EL2.50'})
        BW_dic.update({228: 'EL2.51'})
        BW_dic.update({229: 'EL2.52'})
        BW_dic.update({230: 'EL2.53'})
        BW_dic.update({231: 'EL2.54'})
        # beginnings
        assert BW_dic[62]  == 'TM1.20'
        assert BW_dic[108] == 'TM2.38'
        assert BW_dic[145] == 'TM3.22'
        # ?  188 not 4.38?
        assert BW_dic[188] == 'TM4.39'
        assert BW_dic[232] == 'TM5.36'
        assert BW_dic[316] == 'TM6.28'
        assert BW_dic[359] == 'TM7.32'
        # endings
        assert BW_dic[102] == 'TM1.60'
        assert BW_dic[136] == 'TM2.66'
        assert BW_dic[179] == 'TM3.56'
        assert BW_dic[212] == 'TM4.63'
        assert BW_dic[264] == 'TM5.68'
        assert BW_dic[348] == 'TM6.60'
        assert BW_dic[398] == 'TM7.71'
    elif gpcrin == "HRH1":
        BW_dic.update(get_convert(15, 55, 60, 'TM1'))
        BW_dic.update(get_convert(59, 87, 66, 'TM2'))
        BW_dic.update(get_convert(97, 131, 56, 'TM3'))
        BW_dic.update(get_convert(140, 164, 62, 'TM4'))
        BW_dic.update(get_convert(188, 220, 68, 'TM5'))
        BW_dic.update(get_convert(408, 440, 60, 'TM6'))
        BW_dic.update(get_convert(447, 485, 70, 'TM7'))
        BW_dic.update({93: 'EL1.50'})
        BW_dic.update({180: 'EL2.50'})
        BW_dic.update({181: 'EL2.51'})
        BW_dic.update({182: 'EL2.52'})
        BW_dic.update({183: 'EL2.53'})
        BW_dic.update({184: 'EL2.54'})
        BW_dic.update({185: 'EL2.55'})
        BW_dic.update({186: 'EL2.56'})
        BW_dic.update({187: 'EL2.57'})
        # beginnings
        assert BW_dic[15] == 'TM1.20'
        assert BW_dic[59] == 'TM2.38'
        assert BW_dic[97] == 'TM3.22'
        assert BW_dic[140] == 'TM4.38'
        assert BW_dic[188] == 'TM5.36'
        assert BW_dic[408] == 'TM6.28'
        assert BW_dic[447] == 'TM7.32'
        # endings
        assert BW_dic[55] == 'TM1.60'
        assert BW_dic[87] == 'TM2.66'
        assert BW_dic[131] == 'TM3.56'
        assert BW_dic[164] == 'TM4.62'
        assert BW_dic[220] == 'TM5.68'
        assert BW_dic[440] == 'TM6.60'
        assert BW_dic[485] == 'TM7.70'
    elif gpcrin == "bOPSD":
        BW_dic.update(get_convert(35, 69, 64, 'TM1'))
        BW_dic.update(get_convert(71, 100, 67, 'TM2'))
        BW_dic.update(get_convert(107, 141, 56, 'TM3'))
        BW_dic.update(get_convert(149, 173, 62, 'TM4'))
        BW_dic.update(get_convert(200, 237, 72, 'TM5'))
        BW_dic.update(get_convert(240, 277, 60, 'TM6'))
        BW_dic.update(get_convert(285, 323, 70, 'TM7'))
        BW_dic.update({187: 'EL2.50'})
        BW_dic.update({188: 'EL2.51'})
        BW_dic.update({189: 'EL2.52'})
        BW_dic.update({190: 'EL2.53'})
        BW_dic.update({191: 'EL2.54'})
        # beginnings
        assert BW_dic[35] == 'TM1.30'
        assert BW_dic[71] == 'TM2.38'
        assert BW_dic[107] == 'TM3.22'
        assert BW_dic[149] == 'TM4.38'
        assert BW_dic[200] == 'TM5.35'
        assert BW_dic[240] == 'TM6.23'
        assert BW_dic[286] == 'TM7.33'
        # endings
        assert BW_dic[69] == 'TM1.64'
        assert BW_dic[100] == 'TM2.67'
        assert BW_dic[141] == 'TM3.56'
        assert BW_dic[173] == 'TM4.62'
        assert BW_dic[237] == 'TM5.72'
        assert BW_dic[277] == 'TM6.60'
        assert BW_dic[323] == 'TM7.70'
    # elif gpcrin == "5HT1B":
    #     BW_dic.update(get_convert(37, 77, 60, 'TM1'))
    #     BW_dic.update(get_convert(83, 111, 66, 'TM2'))
    #     BW_dic.update(get_convert(119, 153, 56, 'TM3'))
    #     BW_dic.update(get_convert(162, 187, 62, 'TM4'))
    #     BW_dic.update(get_convert(186, 238, 68, 'TM5'))
    #     BW_dic.update(get_convert(322, 340, 61, 'TM6'))
    #     BW_dic.update(get_convert(361, 385, 70, 'TM7'))
    #     BW_dic.update({115: 'EL1.50'})
    #     BW_dic.update({199: 'EL2.50'})
    #     BW_dic.update({200: 'EL2.51'})
    #     BW_dic.update({201: 'EL2.52'})
    #     BW_dic.update({202: 'EL2.53'})
    #     BW_dic.update({203: 'EL2.54'})
    #     # beginnings
    #     assert BW_dic[37] == 'TM1.20'
    #     assert BW_dic[83] == 'TM2.38'
    #     assert BW_dic[119] == 'TM3.22'
    #     assert BW_dic[162] == 'TM4.38'
    #     assert BW_dic[206] == 'TM5.36'
    #     assert BW_dic[307] == 'TM6.28'
    #     assert BW_dic[347] == 'TM7.31'
    #     # endings
    #     assert BW_dic[77] == 'TM1.60'
    #     assert BW_dic[111] == 'TM2.66'
    #     assert BW_dic[153] == 'TM3.56'
    #     assert BW_dic[187] == 'TM4.62'
    #     assert BW_dic[238] == 'TM5.68'
    #     assert BW_dic[339] == 'TM6.60'
    #     assert BW_dic[348] == 'TM7.32'
    elif gpcrin == "ADA2A":
        BW_dic.update(get_convert(36, 76, 60, 'TM1'))
        BW_dic.update(get_convert(82, 110, 66, 'TM2'))
        BW_dic.update(get_convert(118, 152, 56, 'TM3'))
        BW_dic.update(get_convert(161, 185, 62, 'TM4'))
        BW_dic.update(get_convert(209, 241, 68, 'TM5'))
        BW_dic.update(get_convert(382, 415, 61, 'TM6'))
        BW_dic.update(get_convert(418, 459, 71, 'TM7'))
        BW_dic.update({114: 'EL1.50'})
        BW_dic.update({203: 'EL2.50'})
        BW_dic.update({204: 'EL2.51'})
        BW_dic.update({205: 'EL2.52'})
        BW_dic.update({206: 'EL2.53'})
        BW_dic.update({207: 'EL2.54'})
        # beginnings
        assert BW_dic[36] == 'TM1.20'
        assert BW_dic[82] == 'TM2.38'
        assert BW_dic[118] == 'TM3.22'
        assert BW_dic[161] == 'TM4.38'
        assert BW_dic[209] == 'TM5.36'
        assert BW_dic[382] == 'TM6.28'
        assert BW_dic[419] == 'TM7.31'
        # endings
        assert BW_dic[76] == 'TM1.60'
        assert BW_dic[110] == 'TM2.66'
        assert BW_dic[152] == 'TM3.56'
        assert BW_dic[185] == 'TM4.62'
        assert BW_dic[241] == 'TM5.68'
        assert BW_dic[414] == 'TM6.60'
        assert BW_dic[420] == 'TM7.32'
    elif gpcrin == "5HT5A":
        BW_dic.update(get_convert(28, 68, 60, 'TM1'))
        BW_dic.update(get_convert(74, 102, 66, 'TM2'))
        BW_dic.update(get_convert(111, 145, 56, 'TM3'))
        BW_dic.update(get_convert(154, 179, 62, 'TM4'))
        BW_dic.update(get_convert(198, 230, 68, 'TM5'))
        BW_dic.update(get_convert(278, 311, 61, 'TM6'))
        BW_dic.update(get_convert(316, 354, 69, 'TM7'))
        BW_dic.update({107: 'EL1.50'})
        BW_dic.update({192: 'EL2.50'})
        BW_dic.update({193: 'EL2.51'})
        BW_dic.update({194: 'EL2.52'})
        BW_dic.update({195: 'EL2.53'})
        BW_dic.update({196: 'EL2.54'})
        # beginnings
        assert BW_dic[28] == 'TM1.20'
        assert BW_dic[74] == 'TM2.38'
        assert BW_dic[111] == 'TM3.22'
        assert BW_dic[154] == 'TM4.37'
        assert BW_dic[198] == 'TM5.36'
        assert BW_dic[278] == 'TM6.28'
        assert BW_dic[316] == 'TM7.31'
        # endings
        assert BW_dic[68] == 'TM1.60'
        assert BW_dic[102] == 'TM2.66'
        assert BW_dic[145] == 'TM3.56'
        assert BW_dic[179] == 'TM4.62'
        assert BW_dic[230] == 'TM5.68'
        assert BW_dic[310] == 'TM6.60'
        assert BW_dic[317] == 'TM7.32'
    else:
        print("Need to provide GPCR_input. Using DRD3, DRD2, agrl3_mouse or agrl3_human.")
    return BW_dic



def get_key(my_dict, val):
    '''
    :param my_dict:
    :param val:
    :return:
    '''
    for key, value in my_dict.items():
        if val == value:
            return key


# covert bw index (receptor only) into residue id
def convert_bwidx_to_resid(list_bwidx, receptor_df):
    list_resid=list(prep_bindingsite_df(list_bwidx, receptor_df)['resid'])
    return list_resid

# covert residue id (receptor only) into bw index
def convert_resid_to_bwidx(list_resid, receptor_df):
    list_bwidx=[]
    for i in range(len(list_resid)):
        list_bwidx.append(list(receptor_df[receptor_df['resid']==int(list_resid[i])]['bwid'])[0])
    list_bwidx
    return list_bwidx

def prep_protein_df(bwtable,_sel,_sel_CA):

    receptor_chain=_sel_CA.chain
    receptor_resname=_sel_CA.resname
    receptor_index=_sel_CA.resid
    receptor_bwindex=[]
    for i in range(len(_sel_CA.resid)):
        try:
            receptor_bwindex.append(bwtable[_sel_CA.resid[i]])
        except:
            receptor_bwindex.append(_sel_CA.resid[i])

    # receptor dataframe
    receptor_df = pd.DataFrame({"chain": receptor_chain,
                                "resname": receptor_resname,
                                "resid": receptor_index,
                                "bwid": receptor_bwindex})
    return receptor_df



def prep_receptor_df(molid, bwtable):
    # atom selection for receptor
    receptor = atomsel(vmd_prot_R)
    receptor_CA = atomsel(vmd_prot_R+vmd_add_CA)
    receptor_backbone = atomsel(vmd_prot_R+vmd_add_backbone)

    # commom df
    receptor_df = prep_protein_df(bwtable, receptor, receptor_CA)
    return receptor_df


def prep_bindingsite_df(_bindingsite_bwindex, receptor_df):
    # print(_bindingsite_bwindex)
    # print("receptor_df")
    #_bindingsite=OBS
    _list_chain=[]
    _list_resid=[]
    _list_bwid=[]
    _list_resname=[]
    for i in range(len(_bindingsite_bwindex)):
        if "TM" in _bindingsite_bwindex[i]:
            _bwid=_bindingsite_bwindex[i]
        if "EL" in _bindingsite_bwindex[i]:
            _bwid = _bindingsite_bwindex[i]
        else:
            _bwid="TM" + _bindingsite_bwindex[i]
        try:
            _resid=list(receptor_df[receptor_df['bwid']==_bwid]['resid'])[0]
            _resname=list(receptor_df[receptor_df['resid']==_resid]['resname'])[0]
            _chain=list(receptor_df[receptor_df['resid']==_resid]['chain'])[0]
        except:
            _resid=list(receptor_df[receptor_df['bwid']==_bwid]['resid'])
            _resname=list(receptor_df[receptor_df['resid']==_resid]['resname'])
            _chain=list(receptor_df[receptor_df['resid']==_resid]['chain'])
        _list_chain.append(_chain)
        _list_resid.append(_resid)
        _list_bwid.append(_bwid)
        _list_resname.append(_resname)
        #print(i,_bindingsite[i],_bwid,_resid,_resname,_chain)

    # print("df")
    _df = pd.DataFrame({"chain": _list_chain,
                        "resname": _list_resname,
                        "resid": _list_resid,
                        "bwid": _list_bwid})
    return _df


# move this one to argv
# bw table
#bwtable=convert_Residue2TMs(gpcrin="DRD2")






def resid_to_bwidx(dict_resid_bwidx, ls_resid):
    """

    Args:
        dict_resid_bwidx: preparing dictionary for paticular gpcr system
        ls_resid:  list of residue ids that will be convert into BW index

    Returns: BW index will be returns

    """
    # dict_resid_bwidx: dictionary of resid verse bw_index
    ls_bwidx=[]
    for resid in ls_resid:
        try:
            bwidx=dict_resid_bwidx[resid]
        except:
            bwidx=resid
        ls_bwidx.append(bwidx)
    return ls_bwidx



def resid_to_resname(dict_resid_resname, ls_resid):
    ls_resname=[]
    for resid in ls_resid:
        try:
            resname=dict_resid_resname[int(resid)]
        except:
            resname=resid
        ls_resname.append(resname)
    return ls_resname


def combind_idx_resn(ls_idx,ls_resn):
    ls_out=[]
    for i in range(len(ls_idx)):
        ls_out.append(str(ls_idx[i])+"_"+str(ls_resn[i]))
    return ls_out




## general atom selection
prot = 'protein'
protbb = 'protein and backbone'
protca = 'protein and name CA'
lig = 'not protein and not resname SOD CLA'
lig_noh = 'not type H and not protein and not resname SOD CLA'
ions3 = 'resname SOD CLA'
lig_neighbor = 'protein and not type H and around 7 (not protein and not resname SOD CLA)'


# il2 str & end
il2_str=139
il2_end=146




# binding sites
OBS = ["3.32","3.33","3.36","EL2.52","5.42","5.46","6.48","6.51","6.52","6.55","7.39","7.43"] # from review paper
# OBS=['3.32','3.33','3.36','3.37','5.42','5.43','5.46','6.48','6.51','6.52','6.55','7.35','7.38','7.39','7.42','7.43']
# SBP=["2.61","2.64","2.65","EL1.50","3.28","3.29","3.37","3.40","4.57","5.38","5.39","5.43","5.47","6.44","6.56","6.58","6.59","7.32","7.35","7.36","7.40"]
SBP = ["2.61","2.64","2.65","3.28","3.29","3.37","3.40","4.57","5.38","5.39","5.43","5.47","5.50","6.44","6.56","6.58","6.59","7.32","7.35","7.36","7.40"]
# add 5.50 (2022-01-22)

# more binding sites from ligand contacts
#OBSplus=["3.32","3.33","3.36","3.37","3.40","183","5.38","5.42","5.43","5.46","6.48","6.51","6.52","6.55","7.35","7.38","7.39","7.40","7.42","7.43"]
OBSplus=["3.32","3.33","3.36","3.37","3.40","5.38","5.42","5.43","5.46","6.48","6.51","6.52","6.55","7.35","7.38","7.39","7.40","7.42","7.43"]

# binding site of interested
#BSI = OBS + SBP
# BSI = ['2.61', '2.64', '2.65', '3.28', '3.29', '3.32', '3.33', '3.36', '3.37', '3.40', '4.57', 'EL2.52', '5.38', '5.39', '5.42', '5.43', '5.46', '5.47', '5.50', '6.44', '6.48', '6.51', '6.52', '6.55', '6.56', '6.58', '6.59', '7.32', '7.35', '7.36', '7.39', '7.40', '7.43']
# add EL1.50 & 7.42
BSI = ['2.61', '2.64', '2.65', 'EL1.50', '3.28', '3.29', '3.32', '3.33', '3.36', '3.37', '3.40', '4.57', 'EL2.52', '5.38', '5.39', '5.42', '5.43', '5.46', '5.47', '5.50', '6.44', '6.48', '6.51', '6.52', '6.55', '6.56', '6.58', '6.59', '7.32', '7.35', '7.36', '7.39', '7.40', '7.42', '7.43']

def get_rec_info(bwtable):
    rec_segments_pia=[str(get_key(bwtable, "TM1.30"))+":"+str(get_key(bwtable, "TM1.36")),
                      str(get_key(bwtable, "TM1.37"))+":"+str(get_key(bwtable, "TM1.45")),
                      str(get_key(bwtable, "TM1.46"))+":"+str(get_key(bwtable, "TM1.59")),
                      str(get_key(bwtable, "TM2.38"))+":"+str(get_key(bwtable, "TM2.51")),
                      str(get_key(bwtable, "TM2.52"))+":"+str(get_key(bwtable, "TM2.60")),
                      str(get_key(bwtable, "TM2.61"))+":"+str(get_key(bwtable, "TM2.66")),
                      str(get_key(bwtable, "TM3.22"))+":"+str(get_key(bwtable, "TM3.35")),
                      str(get_key(bwtable, "TM3.36"))+":"+str(get_key(bwtable, "TM3.40")),
                      str(get_key(bwtable, "TM3.41"))+":"+str(get_key(bwtable, "TM3.55")),
                      str(get_key(bwtable, "TM4.39"))+":"+str(get_key(bwtable, "TM4.49")),
                      str(get_key(bwtable, "TM4.50"))+":"+str(get_key(bwtable, "TM4.55")),
                      str(get_key(bwtable, "TM4.56"))+":"+str(get_key(bwtable, "TM4.62")),
                      str(get_key(bwtable, "TM5.36"))+":"+str(get_key(bwtable, "TM5.45")),
                      str(get_key(bwtable, "TM5.46"))+":"+str(get_key(bwtable, "TM5.50")),
                      str(get_key(bwtable, "TM5.51"))+":"+str(get_key(bwtable, "TM5.63")),
                      str(get_key(bwtable, "TM6.30"))+":"+str(get_key(bwtable, "TM6.43")),
                      str(get_key(bwtable, "TM6.44"))+":"+str(get_key(bwtable, "TM6.48")),
                      str(get_key(bwtable, "TM6.49"))+":"+str(get_key(bwtable, "TM6.60")),
                      str(get_key(bwtable, "TM7.32"))+":"+str(get_key(bwtable, "TM7.43")),
                      str(get_key(bwtable, "TM7.44"))+":"+str(get_key(bwtable, "TM7.54")),
                      str(get_key(bwtable, "TM7.58"))+":"+str(get_key(bwtable, "TM7.67"))]
    rec_labels_pia=["1e", "1m", "1i",
                    "2i", "2m", "2e",
                    "3e", "3m", "3i",
                    "4i", "4m", "4e",
                    "5e", "5m", "5i",
                    "6i", "6m", "6e",
                    "7e", "7i", "8h"]
    rec_chainid="R"
    return rec_segments_pia, rec_labels_pia, rec_chainid


ga_labels=['H2m','H3m','HGm','H4','H5','H5n','H5c','HN2t']
gi_segments=['211:214','242:255','271:280','294:309','329:353','329:333','347:351','24:32']
go_segments=['212:215','243:256','272:281','295:310','329:353','329:333','347:351','24:32']
ga_chainid="A"



def check_output_dir(output_dir="output"):
    '''
    :param output_dir: make sure output_dir exists
    '''
    if not os.path.exists(output_dir):
        Path(output_dir).mkdir(parents=True, exist_ok=True)

