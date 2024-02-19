#!/usr/bin/env python
"""
CXCR2_HUMAN TMs table
"""
def convert_Residue2TMs():
    CXCR2_HUMAN_BW_dic = {}
    for i in range(46-5, 75+2): #TM1
        CXCR2_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(46-5))

    for i in range(82-5, 111+2): #TM2
        CXCR2_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(82-5))

    for i in range(116-5, 148+2): #TM3
        CXCR2_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(116-5))

    for i in range(160-5, 182+2): #TM4
        CXCR2_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(160-5))

    for i in range(208-5, 233+2): #TM5
        CXCR2_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(208-5))

    for i in range(246-5, 276+2): #TM6
        CXCR2_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(246-5))

    for i in range(294-5, 329+2): #TM7
        CXCR2_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(294-5))

    assert CXCR2_HUMAN_BW_dic[66] == 'TM1.50'
    assert CXCR2_HUMAN_BW_dic[94] == 'TM2.50'
    assert CXCR2_HUMAN_BW_dic[144] == 'TM3.50'
    assert CXCR2_HUMAN_BW_dic[170] == 'TM4.50'
    assert CXCR2_HUMAN_BW_dic[223] == 'TM5.50'
    assert CXCR2_HUMAN_BW_dic[266] == 'TM6.50'
    assert CXCR2_HUMAN_BW_dic[311] == 'TM7.50'

    return CXCR2_HUMAN_BW_dic
