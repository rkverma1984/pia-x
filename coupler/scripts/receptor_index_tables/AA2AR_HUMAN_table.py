#!/usr/bin/env python
"""
AA2AR_HUMAN TMs table
"""
def convert_Residue2TMs():
    AA2AR_HUMAN_BW_dic = {}
    for i in range(4-5, 33+2): #TM1
        AA2AR_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(4-5))

    for i in range(40-5, 69+2): #TM2
        AA2AR_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(40-5))

    for i in range(74-5, 106+2): #TM3
        AA2AR_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(74-5))

    for i in range(119-5, 141+2): #TM4
        AA2AR_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(119-5))

    for i in range(174-5, 199+2): #TM5
        AA2AR_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(174-5))

    for i in range(228-5, 258+2): #TM6
        AA2AR_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(228-5))

    for i in range(268-5, 303+2): #TM7
        AA2AR_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(268-5))

    assert AA2AR_HUMAN_BW_dic[24] == 'TM1.50'
    assert AA2AR_HUMAN_BW_dic[52] == 'TM2.50'
    assert AA2AR_HUMAN_BW_dic[102] == 'TM3.50'
    assert AA2AR_HUMAN_BW_dic[129] == 'TM4.50'
    assert AA2AR_HUMAN_BW_dic[189] == 'TM5.50'
    assert AA2AR_HUMAN_BW_dic[248] == 'TM6.50'
    assert AA2AR_HUMAN_BW_dic[285] == 'TM7.50'

    return AA2AR_HUMAN_BW_dic
