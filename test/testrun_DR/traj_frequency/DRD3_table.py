#!/usr/bin/env python
"""
DRD3 TMs table
"""
def get_convert(resid_str=31,resid_end=62,BW_end=60,TMx='TM1'):
    GPCR_BW_dic = {}
    for i in range(resid_end, resid_str-1, -1): #TM1
        GPCR_BW_dic[i]=TMx+"."+str(BW_end)
        BW_end=BW_end-1
    return GPCR_BW_dic

def convert_Residue2TMs():
    DRD3_BW_dic = {}
    DRD3_BW_dic.update(get_convert(27,57,60,'TM1'))
    DRD3_BW_dic.update(get_convert(63,91,66,'TM2'))
    DRD3_BW_dic.update(get_convert(100,134,56,'TM3'))
    DRD3_BW_dic.update(get_convert(145,170,62,'TM4'))
    DRD3_BW_dic.update(get_convert(186,218,68,'TM5'))
    DRD3_BW_dic.update(get_convert(322,354,60,'TM6'))
    DRD3_BW_dic.update(get_convert(362,400,70,'TM7'))
    assert DRD3_BW_dic[57] == 'TM1.60'
    assert DRD3_BW_dic[91] == 'TM2.66'
    assert DRD3_BW_dic[134] == 'TM3.56'
    assert DRD3_BW_dic[170] == 'TM4.62'
    assert DRD3_BW_dic[218] == 'TM5.68'
    assert DRD3_BW_dic[354] == 'TM6.60'
    assert DRD3_BW_dic[362] == 'TM7.32'
    return DRD3_BW_dic
