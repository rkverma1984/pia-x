#!/usr/bin/env python
"""
ACM2_HUMAN TMs table
"""
def convert_Residue2TMs():
    ACM2_HUMAN_BW_dic = {}
    for i in range(21-5, 50+2): #TM1
        ACM2_HUMAN_BW_dic[i]='TM1.'+str(30-5+i-(21-5))

    for i in range(57-5, 86+2): #TM2
        ACM2_HUMAN_BW_dic[i]='TM2.'+str(38-5+i-(57-5))

    for i in range(93-5, 125+2): #TM3
        ACM2_HUMAN_BW_dic[i]='TM3.'+str(22-5+i-(93-5))

    for i in range(138-5, 160+2): #TM4
        ACM2_HUMAN_BW_dic[i]='TM4.'+str(40-5+i-(138-5))

    for i in range(183-5, 208+2): #TM5
        ACM2_HUMAN_BW_dic[i]='TM5.'+str(35-5+i-(183-5))

    for i in range(382-5, 412+2): #TM6
        ACM2_HUMAN_BW_dic[i]='TM6.'+str(30-5+i-(382-5))

    for i in range(420-5, 455+2): #TM7
        ACM2_HUMAN_BW_dic[i]='TM7.'+str(33-5+i-(420-5))

    assert ACM2_HUMAN_BW_dic[41] == 'TM1.50'
    assert ACM2_HUMAN_BW_dic[69] == 'TM2.50'
    assert ACM2_HUMAN_BW_dic[121] == 'TM3.50'
    assert ACM2_HUMAN_BW_dic[148] == 'TM4.50'
    assert ACM2_HUMAN_BW_dic[198] == 'TM5.50'
    assert ACM2_HUMAN_BW_dic[402] == 'TM6.50'
    assert ACM2_HUMAN_BW_dic[437] == 'TM7.50'

    return ACM2_HUMAN_BW_dic
