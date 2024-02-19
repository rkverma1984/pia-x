#!/usr/bin/env python
"""
ADA2B_HUMAN TMs table
"""
def convert_Residue2TMs():
    ADA2B_HUMAN_BW_dic = {}
    for i in range(10-5, 39+2): #TM1
        ADA2B_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(10-5))

    for i in range(46-5, 75+2): #TM2
        ADA2B_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(46-5))

    for i in range(82-5, 114+2): #TM3
        ADA2B_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(82-5))

    for i in range(127-5, 149+2): #TM4
        ADA2B_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(127-5))

    for i in range(169-5, 194+2): #TM5
        ADA2B_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(169-5))

    for i in range(366-5, 396+2): #TM6
        ADA2B_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(366-5))

    for i in range(406-5, 441+2): #TM7
        ADA2B_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(406-5))

    assert ADA2B_HUMAN_BW_dic[30] == 'TM1.50'
    assert ADA2B_HUMAN_BW_dic[58] == 'TM2.50'
    assert ADA2B_HUMAN_BW_dic[110] == 'TM3.50'
    assert ADA2B_HUMAN_BW_dic[137] == 'TM4.50'
    assert ADA2B_HUMAN_BW_dic[184] == 'TM5.50'
    assert ADA2B_HUMAN_BW_dic[386] == 'TM6.50'
    assert ADA2B_HUMAN_BW_dic[423] == 'TM7.50'

    return ADA2B_HUMAN_BW_dic
