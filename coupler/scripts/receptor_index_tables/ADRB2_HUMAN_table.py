#!/usr/bin/env python
"""
ADRB2_HUMAN TMs table
"""
def convert_Residue2TMs():
    ADRB2_HUMAN_BW_dic = {}
    for i in range(31-5, 60+2): #TM1
        ADRB2_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(31-5))

    for i in range(67-5, 96+2): #TM2
        ADRB2_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(67-5))

    for i in range(103-5, 135+2): #TM3
        ADRB2_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(103-5))

    for i in range(148-5, 170+2): #TM4
        ADRB2_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(148-5))

    for i in range(196-5, 221+2): #TM5
        ADRB2_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(196-5))

    for i in range(268-5, 298+2): #TM6
        ADRB2_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(268-5))

    for i in range(306-5, 341+2): #TM7
        ADRB2_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(306-5))

    assert ADRB2_HUMAN_BW_dic[51] == 'TM1.50'
    assert ADRB2_HUMAN_BW_dic[79] == 'TM2.50'
    assert ADRB2_HUMAN_BW_dic[131] == 'TM3.50'
    assert ADRB2_HUMAN_BW_dic[158] == 'TM4.50'
    assert ADRB2_HUMAN_BW_dic[211] == 'TM5.50'
    assert ADRB2_HUMAN_BW_dic[288] == 'TM6.50'
    assert ADRB2_HUMAN_BW_dic[323] == 'TM7.50'

    return ADRB2_HUMAN_BW_dic
