#!/usr/bin/env python
"""
CCR6_HUMAN TMs table
"""
def convert_Residue2TMs():
    CCR6_HUMAN_BW_dic = {}
    for i in range(44-5, 73+2): #TM1
        CCR6_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(44-5))

    for i in range(80-5, 109+2): #TM2
        CCR6_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(80-5))

    for i in range(115-5, 147+2): #TM3
        CCR6_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(115-5))

    for i in range(162-5, 184+2): #TM4
        CCR6_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(162-5))

    for i in range(211-5, 236+2): #TM5
        CCR6_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(211-5))

    for i in range(249-5, 279+2): #TM6
        CCR6_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(249-5))

    for i in range(296-5, 331+2): #TM7
        CCR6_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(296-5))

    assert CCR6_HUMAN_BW_dic[64] == 'TM1.50'
    assert CCR6_HUMAN_BW_dic[92] == 'TM2.50'
    assert CCR6_HUMAN_BW_dic[143] == 'TM3.50'
    assert CCR6_HUMAN_BW_dic[172] == 'TM4.50'
    assert CCR6_HUMAN_BW_dic[226] == 'TM5.50'
    assert CCR6_HUMAN_BW_dic[269] == 'TM6.50'
    assert CCR6_HUMAN_BW_dic[313] == 'TM7.50'

    return CCR6_HUMAN_BW_dic
