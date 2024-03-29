#!/usr/bin/env python
"""
CNR1_HUMAN TMs table
"""
def convert_Residue2TMs():
    CNR1_HUMAN_BW_dic = {}
    for i in range(114-5, 143+2+1): #TM1
        CNR1_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(114-5))

    for i in range(151-5, 180+2): #TM2
        CNR1_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(151-5))

    for i in range(186-5, 228): #TM3  #Manually
        CNR1_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(186-5))

    for i in range(231-3, 253+2): #TM4 #Manually
        CNR1_HUMAN_BW_dic[i]='TM4.'+str(40-3+i-(231-3))

    for i in range(271-5, 324): #TM5  #Manually
        CNR1_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(271-5))

    for i in range(338-14, 368+2): #TM6 #Manually
        CNR1_HUMAN_BW_dic[i]='TM6.'+str(30-14+i-(338-14))

    for i in range(377-5, 412+2): #TM7
        CNR1_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(377-5))

    assert CNR1_HUMAN_BW_dic[134] == 'TM1.50'
    assert CNR1_HUMAN_BW_dic[163] == 'TM2.50'
    assert CNR1_HUMAN_BW_dic[214] == 'TM3.50'
    assert CNR1_HUMAN_BW_dic[241] == 'TM4.50'
    assert CNR1_HUMAN_BW_dic[286] == 'TM5.50'
    assert CNR1_HUMAN_BW_dic[358] == 'TM6.50'
    assert CNR1_HUMAN_BW_dic[394] == 'TM7.50'

    return CNR1_HUMAN_BW_dic
