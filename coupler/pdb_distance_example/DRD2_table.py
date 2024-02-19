#!/usr/bin/env python
"""
DRD2 TMs table
"""
def convert_Residue2TMs():
    DRD2_BW_dic = {}
    for i in range(32-5, 61+2): #TM1
        DRD2_BW_dic[i]='TM1.'+str(30-5+i-(32-5))

    for i in range(68-5, 97+2): #TM2
        DRD2_BW_dic[i]='TM2.'+str(38-5+i-(68-5))

    for i in range(104-5, 136+2): #TM3
        DRD2_BW_dic[i]='TM3.'+str(22-5+i-(104-5))

    for i in range(150-5, 172+2): #TM4
        DRD2_BW_dic[i]='TM4.'+str(40-5+i-(150-5))

    for i in range(186-5, 211+2): #TM5
        DRD2_BW_dic[i]='TM5.'+str(35-5+i-(186-5))

    for i in range(368-5, 398+2): #TM6
        DRD2_BW_dic[i]='TM6.'+str(30-5+i-(368-5))

    for i in range(406-5, 441+2): #TM7
        DRD2_BW_dic[i]='TM7.'+str(33-5+i-(406-5))

    assert DRD2_BW_dic[52] == 'TM1.50'
    assert DRD2_BW_dic[80] == 'TM2.50'
    assert DRD2_BW_dic[132] == 'TM3.50'
    assert DRD2_BW_dic[160] == 'TM4.50'
    assert DRD2_BW_dic[201] == 'TM5.50'
    assert DRD2_BW_dic[388] == 'TM6.50'
    assert DRD2_BW_dic[423] == 'TM7.50'

    return DRD2_BW_dic
