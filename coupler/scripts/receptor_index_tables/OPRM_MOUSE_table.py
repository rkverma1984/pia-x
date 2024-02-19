#!/usr/bin/env python
"""
OPRM_MOUSE TMs table
"""
def convert_Residue2TMs():
    OPRM_MOUSE_BW_dic = {}
    for i in range(66-5, 95+2): #TM1
        OPRM_MOUSE_BW_dic[i]='TM1.'+str(30-5+i-(66-5))

    for i in range(102-5, 131+2): #TM2
        OPRM_MOUSE_BW_dic[i]='TM2.'+str(38-5+i-(102-5))

    for i in range(137-5, 169+2): #TM3
        OPRM_MOUSE_BW_dic[i]='TM3.'+str(22-5+i-(137-5))

    for i in range(182-5, 204+2): #TM4
        OPRM_MOUSE_BW_dic[i]='TM4.'+str(40-5+i-(182-5))

    for i in range(229-5, 254+2): #TM5
        OPRM_MOUSE_BW_dic[i]='TM5.'+str(35-5+i-(229-5))

    for i in range(275-5, 305+2): #TM6
        OPRM_MOUSE_BW_dic[i]='TM6.'+str(30-5+i-(275-5))

    for i in range(316-5, 351+2): #TM7
        OPRM_MOUSE_BW_dic[i]='TM7.'+str(33-5+i-(316-5))

    assert OPRM_MOUSE_BW_dic[86] == 'TM1.50'
    assert OPRM_MOUSE_BW_dic[114] == 'TM2.50'
    assert OPRM_MOUSE_BW_dic[165] == 'TM3.50'
    assert OPRM_MOUSE_BW_dic[192] == 'TM4.50'
    assert OPRM_MOUSE_BW_dic[244] == 'TM5.50'
    assert OPRM_MOUSE_BW_dic[295] == 'TM6.50'
    assert OPRM_MOUSE_BW_dic[333] == 'TM7.50'

    return OPRM_MOUSE_BW_dic
