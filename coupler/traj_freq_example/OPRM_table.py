#!/usr/bin/env python
"""
OPRM TMs table
"""
def convert_Residue2TMs():
    OPRM_BW_dic = {}
    for i in range(66-5+2, 95+2+2): #TM1
        OPRM_BW_dic[i]='TM1.'+str(30-5+i-(66-5+2))

    for i in range(102-5+2, 131+2+2): #TM2
        OPRM_BW_dic[i]='TM2.'+str(38-5+i-(102-5+2))

    for i in range(137-5+2, 169+2+2): #TM3
        OPRM_BW_dic[i]='TM3.'+str(22-5+i-(137-5+2))

    for i in range(182-5+2, 204+2+2): #TM4
        OPRM_BW_dic[i]='TM4.'+str(40-5+i-(182-5+2))

    for i in range(229-5+2, 254+2+2): #TM5
        OPRM_BW_dic[i]='TM5.'+str(35-5+i-(229-5+2))

    for i in range(275-5+2, 305+2+2): #TM6
        OPRM_BW_dic[i]='TM6.'+str(30-5+i-(275-5+2))

    for i in range(316-5+2, 351+2+2): #TM7
        OPRM_BW_dic[i]='TM7.'+str(33-5+i-(316-5+2))

    assert OPRM_BW_dic[86+2] == 'TM1.50'
    assert OPRM_BW_dic[114+2] == 'TM2.50'
    assert OPRM_BW_dic[165+2] == 'TM3.50'
    assert OPRM_BW_dic[192+2] == 'TM4.50'
    assert OPRM_BW_dic[244+2] == 'TM5.50'
    assert OPRM_BW_dic[295+2] == 'TM6.50'
    assert OPRM_BW_dic[333+2] == 'TM7.50'

    return OPRM_BW_dic
