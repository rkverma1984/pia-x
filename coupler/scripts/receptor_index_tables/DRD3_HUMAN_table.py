#!/usr/bin/env python
"""
DRD3_HUMAN TMs table
"""
def convert_Residue2TMs():
    DRD3_HUMAN_BW_dic = {}
    for i in range(27-5, 56+2): #TM1
        DRD3_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(27-5))

    for i in range(63-5, 92+2): #TM2
        DRD3_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(63-5))

    for i in range(100-5, 132+2): #TM3
        DRD3_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(100-5))

    for i in range(148-5, 170+2): #TM4
        DRD3_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(148-5))

    for i in range(185-5, 210+2): #TM5
        DRD3_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(185-5))

    for i in range(324-5, 354+2): #TM6
        DRD3_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(324-5))

    for i in range(363-5, 398+2): #TM7
        DRD3_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(363-5))

    assert DRD3_HUMAN_BW_dic[47] == 'TM1.50'
    assert DRD3_HUMAN_BW_dic[75] == 'TM2.50'
    assert DRD3_HUMAN_BW_dic[128] == 'TM3.50'
    assert DRD3_HUMAN_BW_dic[158] == 'TM4.50'
    assert DRD3_HUMAN_BW_dic[200] == 'TM5.50'
    assert DRD3_HUMAN_BW_dic[344] == 'TM6.50'
    assert DRD3_HUMAN_BW_dic[380] == 'TM7.50'

    return DRD3_HUMAN_BW_dic
