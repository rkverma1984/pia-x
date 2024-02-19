#!/usr/bin/env python
"""
FPR2_HUMAN TMs table
"""
def convert_Residue2TMs():
    FPR2_HUMAN_BW_dic = {}
    for i in range(24-5, 53+2): #TM1
        FPR2_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(24-5))

    for i in range(59-5, 88+2): #TM2
        FPR2_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(59-5))

    for i in range(95-5, 127+2): #TM3
        FPR2_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(95-5))

    for i in range(140-5, 162+2): #TM4
        FPR2_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(140-5))

    for i in range(198-5, 223+2): #TM5
        FPR2_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(198-5))

    for i in range(236-5, 266+2): #TM6
        FPR2_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(236-5))

    for i in range(282-5, 317+2): #TM7
        FPR2_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(282-5))

    assert FPR2_HUMAN_BW_dic[44] == 'TM1.50'
    assert FPR2_HUMAN_BW_dic[71] == 'TM2.50'
    assert FPR2_HUMAN_BW_dic[123] == 'TM3.50'
    assert FPR2_HUMAN_BW_dic[150] == 'TM4.50'
    assert FPR2_HUMAN_BW_dic[213] == 'TM5.50'
    assert FPR2_HUMAN_BW_dic[256] == 'TM6.50'
    assert FPR2_HUMAN_BW_dic[299] == 'TM7.50'

    return FPR2_HUMAN_BW_dic
