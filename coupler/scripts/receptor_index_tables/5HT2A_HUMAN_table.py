#!/usr/bin/env python
"""
5HT2A_HUMAN TMs table
"""
def convert_Residue2TMs():
    5HT2A_HUMAN_BW_dic = {}
    for i in range(72-5, 101+2): #TM1
        5HT2A_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(72-5))

    for i in range(108-5, 137+2): #TM2
        5HT2A_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(108-5))

    for i in range(145-5, 177+2): #TM3
        5HT2A_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(145-5))

    for i in range(190-5, 212+2): #TM4
        5HT2A_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(190-5))

    for i in range(231-5, 256+2): #TM5
        5HT2A_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(231-5))

    for i in range(318-5, 348+2): #TM6
        5HT2A_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(318-5))

    for i in range(360-5, 395+2): #TM7
        5HT2A_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(360-5))

    assert 5HT2A_HUMAN_BW_dic[92] == 'TM1.50'
    assert 5HT2A_HUMAN_BW_dic[120] == 'TM2.50'
    assert 5HT2A_HUMAN_BW_dic[173] == 'TM3.50'
    assert 5HT2A_HUMAN_BW_dic[200] == 'TM4.50'
    assert 5HT2A_HUMAN_BW_dic[246] == 'TM5.50'
    assert 5HT2A_HUMAN_BW_dic[338] == 'TM6.50'
    assert 5HT2A_HUMAN_BW_dic[377] == 'TM7.50'

    return 5HT2A_HUMAN_BW_dic
