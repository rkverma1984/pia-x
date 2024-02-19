#!/usr/bin/env python
"""
ACM1_HUMAN TMs table
"""
def convert_Residue2TMs():
    ACM1_HUMAN_BW_dic = {}
    for i in range(23-5, 52+2): #TM1
        ACM1_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(23-5))

    for i in range(59-5, 88+2): #TM2
        ACM1_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(59-5))

    for i in range(95-5, 127+2): #TM3
        ACM1_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(95-5))

    for i in range(140-5, 162+2): #TM4
        ACM1_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(140-5))

    for i in range(185-5, 210+2): #TM5
        ACM1_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(185-5))

    for i in range(360-5, 390+2): #TM6
        ACM1_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(360-5))

    for i in range(398-5, 433+2): #TM7
        ACM1_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(398-5))

    assert ACM1_HUMAN_BW_dic[43] == 'TM1.50'
    assert ACM1_HUMAN_BW_dic[71] == 'TM2.50'
    assert ACM1_HUMAN_BW_dic[123] == 'TM3.50'
    assert ACM1_HUMAN_BW_dic[150] == 'TM4.50'
    assert ACM1_HUMAN_BW_dic[200] == 'TM5.50'
    assert ACM1_HUMAN_BW_dic[380] == 'TM6.50'
    assert ACM1_HUMAN_BW_dic[415] == 'TM7.50'

    return ACM1_HUMAN_BW_dic
