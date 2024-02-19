#!/usr/bin/env python
"""
5HT1B_HUMAN TMs table
"""
def convert_Residue2TMs():
    5HT1B_HUMAN_BW_dic = {}
    for i in range(47-5, 76+2): #TM1
        5HT1B_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(47-5))

    for i in range(83-5, 112+2): #TM2
        5HT1B_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(83-5))

    for i in range(119-5, 151+2): #TM3
        5HT1B_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(119-5))

    for i in range(164-5, 186+2): #TM4
        5HT1B_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(164-5))

    for i in range(205-5, 230+2): #TM5
        5HT1B_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(205-5))

    for i in range(309-5, 339+2): #TM6
        5HT1B_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(309-5))

    for i in range(349-5, 384+2): #TM7
        5HT1B_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(349-5))

    assert 5HT1B_HUMAN_BW_dic[67] == 'TM1.50'
    assert 5HT1B_HUMAN_BW_dic[95] == 'TM2.50'
    assert 5HT1B_HUMAN_BW_dic[147] == 'TM3.50'
    assert 5HT1B_HUMAN_BW_dic[174] == 'TM4.50'
    assert 5HT1B_HUMAN_BW_dic[220] == 'TM5.50'
    assert 5HT1B_HUMAN_BW_dic[329] == 'TM6.50'
    assert 5HT1B_HUMAN_BW_dic[366] == 'TM7.50'

    return 5HT1B_HUMAN_BW_dic
