#!/usr/bin/env python
"""
NTR1_HUMAN TMs table
"""
def convert_Residue2TMs():
    NTR1_HUMAN_BW_dic = {}
    for i in range(61-5, 90+2): #TM1
        NTR1_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(61-5))

    for i in range(100-5, 129+2): #TM2
        NTR1_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(100-5))

    for i in range(138-5, 170+2): #TM3
        NTR1_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(138-5))

    for i in range(183-5, 205+2): #TM4
        NTR1_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(183-5))

    for i in range(233-5, 258+2): #TM5
        NTR1_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(233-5))

    for i in range(298-5, 328+2): #TM6
        NTR1_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(298-5))

    for i in range(344-5, 379+2): #TM7
        NTR1_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(344-5))

    assert NTR1_HUMAN_BW_dic[81] == 'TM1.50'
    assert NTR1_HUMAN_BW_dic[112] == 'TM2.50'
    assert NTR1_HUMAN_BW_dic[166] == 'TM3.50'
    assert NTR1_HUMAN_BW_dic[193] == 'TM4.50'
    assert NTR1_HUMAN_BW_dic[248] == 'TM5.50'
    assert NTR1_HUMAN_BW_dic[318] == 'TM6.50'
    assert NTR1_HUMAN_BW_dic[361] == 'TM7.50'

    return NTR1_HUMAN_BW_dic
