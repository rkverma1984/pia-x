#!/usr/bin/env python
"""
AA1R_HUMAN TMs table
"""
def convert_Residue2TMs():
    AA1R_HUMAN_BW_dic = {}
    for i in range(7-5, 36+2): #TM1
        AA1R_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(7-5))

    for i in range(43-5, 72+2): #TM2
        AA1R_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(43-5))

    for i in range(77-5, 109+2): #TM3
        AA1R_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(77-5))

    for i in range(122-5, 144+2): #TM4
        AA1R_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(122-5))

    for i in range(177-5, 202+2): #TM5
        AA1R_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(177-5))

    for i in range(229-5, 259+2): #TM6
        AA1R_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(229-5))

    for i in range(268-5, 303+2): #TM7
        AA1R_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(268-5))

    assert AA1R_HUMAN_BW_dic[27] == 'TM1.50'
    assert AA1R_HUMAN_BW_dic[55] == 'TM2.50'
    assert AA1R_HUMAN_BW_dic[105] == 'TM3.50'
    assert AA1R_HUMAN_BW_dic[132] == 'TM4.50'
    assert AA1R_HUMAN_BW_dic[192] == 'TM5.50'
    assert AA1R_HUMAN_BW_dic[249] == 'TM6.50'
    assert AA1R_HUMAN_BW_dic[285] == 'TM7.50'

    return AA1R_HUMAN_BW_dic
