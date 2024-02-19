#!/usr/bin/env python
"""
CNR2_HUMAN TMs table
"""
def convert_Residue2TMs():
    CNR2_HUMAN_BW_dic = {}
    for i in range(31-5, 60+2): #TM1
        CNR2_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(31-5))

    for i in range(68-5, 97+2): #TM2
        CNR2_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(68-5))

    for i in range(103-5, 135+2): #TM3
        CNR2_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(103-5))

    for i in range(148-5, 170+2): #TM4
        CNR2_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(148-5))

    for i in range(186-5, 211+2): #TM5
        CNR2_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(186-5))

    for i in range(240-5, 270+2): #TM6
        CNR2_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(240-5))

    for i in range(279-5, 314+2): #TM7
        CNR2_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(279-5))

    assert CNR2_HUMAN_BW_dic[51] == 'TM1.50'
    assert CNR2_HUMAN_BW_dic[80] == 'TM2.50'
    assert CNR2_HUMAN_BW_dic[131] == 'TM3.50'
    assert CNR2_HUMAN_BW_dic[158] == 'TM4.50'
    assert CNR2_HUMAN_BW_dic[201] == 'TM5.50'
    assert CNR2_HUMAN_BW_dic[260] == 'TM6.50'
    assert CNR2_HUMAN_BW_dic[296] == 'TM7.50'

    return CNR2_HUMAN_BW_dic
