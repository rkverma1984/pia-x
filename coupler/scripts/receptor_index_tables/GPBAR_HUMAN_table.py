#!/usr/bin/env python
"""
GPBAR_HUMAN TMs table
"""
def convert_Residue2TMs():
    GPBAR_HUMAN_BW_dic = {}
    for i in range(12-5, 41+2): #TM1
        GPBAR_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(12-5))

    for i in range(49-5, 78+2): #TM2
        GPBAR_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(49-5))

    for i in range(82-5, 114+2): #TM3
        GPBAR_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(82-5))

    for i in range(122-5, 144+2): #TM4
        GPBAR_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(122-5))

    for i in range(161-5, 186+2): #TM5
        GPBAR_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(161-5))

    for i in range(219-5, 249+2): #TM6
        GPBAR_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(219-5))

    for i in range(260-5, 295+2): #TM7
        GPBAR_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(260-5))

    assert GPBAR_HUMAN_BW_dic[32] == 'TM1.50'
    assert GPBAR_HUMAN_BW_dic[61] == 'TM2.50'
    assert GPBAR_HUMAN_BW_dic[110] == 'TM3.50'
    assert GPBAR_HUMAN_BW_dic[132] == 'TM4.50'
    assert GPBAR_HUMAN_BW_dic[176] == 'TM5.50'
    assert GPBAR_HUMAN_BW_dic[239] == 'TM6.50'
    assert GPBAR_HUMAN_BW_dic[277] == 'TM7.50'

    return GPBAR_HUMAN_BW_dic
