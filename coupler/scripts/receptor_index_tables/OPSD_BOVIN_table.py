#!/usr/bin/env python
"""
OPSD_BOVIN TMs table
"""
def convert_Residue2TMs():
    OPSD_BOVIN_BW_dic = {}
    for i in range(35-5, 64+2): #TM1
        OPSD_BOVIN_BW_dic[i]='TM1.'+str(30-5+i-(35-5))

    for i in range(71-5, 100+2): #TM2
        OPSD_BOVIN_BW_dic[i]='TM2.'+str(38-5+i-(71-5))

    for i in range(107-5, 139+2): #TM3
        OPSD_BOVIN_BW_dic[i]='TM3.'+str(22-5+i-(107-5))

    for i in range(151-5, 173+2): #TM4
        OPSD_BOVIN_BW_dic[i]='TM4.'+str(40-5+i-(151-5))

    for i in range(200-5, 225+2): #TM5
        OPSD_BOVIN_BW_dic[i]='TM5.'+str(35-5+i-(200-5))

    for i in range(247-5, 277+2): #TM6
        OPSD_BOVIN_BW_dic[i]='TM6.'+str(30-5+i-(247-5))

    for i in range(286-5, 321+2): #TM7
        OPSD_BOVIN_BW_dic[i]='TM7.'+str(33-5+i-(286-5))

    assert OPSD_BOVIN_BW_dic[55] == 'TM1.50'
    assert OPSD_BOVIN_BW_dic[83] == 'TM2.50'
    assert OPSD_BOVIN_BW_dic[135] == 'TM3.50'
    assert OPSD_BOVIN_BW_dic[161] == 'TM4.50'
    assert OPSD_BOVIN_BW_dic[215] == 'TM5.50'
    assert OPSD_BOVIN_BW_dic[267] == 'TM6.50'
    assert OPSD_BOVIN_BW_dic[303] == 'TM7.50'

    return OPSD_BOVIN_BW_dic
