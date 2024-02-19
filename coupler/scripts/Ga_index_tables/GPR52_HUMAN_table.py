#!/usr/bin/env python
"""
GPR52_HUMAN TMs table
"""
def convert_Residue2TMs():
    GPR52_HUMAN_BW_dic = {}
    for i in range(38-5, 67+2): #TM1
        GPR52_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(38-5))

    for i in range(75-5, 104+2): #TM2
        GPR52_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(75-5))

    for i in range(111-5, 143+2): #TM3
        GPR52_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(111-5))

    for i in range(156-5, 178+2): #TM4
        GPR52_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(156-5))

    for i in range(199-5, 224+2): #TM5
        GPR52_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(199-5))

    for i in range(260-5, 290+2): #TM6
        GPR52_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(260-5))

    for i in range(297-5, 332+2): #TM7
        GPR52_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(297-5))

    assert GPR52_HUMAN_BW_dic[58] == 'TM1.50'
    assert GPR52_HUMAN_BW_dic[87] == 'TM2.50'
    assert GPR52_HUMAN_BW_dic[139] == 'TM3.50'
    assert GPR52_HUMAN_BW_dic[166] == 'TM4.50'
    assert GPR52_HUMAN_BW_dic[214] == 'TM5.50'
    assert GPR52_HUMAN_BW_dic[280] == 'TM6.50'
    assert GPR52_HUMAN_BW_dic[314] == 'TM7.50'

    return GPR52_HUMAN_BW_dic
