#!/usr/bin/env python
"""
ADRB1_MELGA TMs table
"""
def convert_Residue2TMs():
    ADRB1_MELGA_BW_dic = {}
    for i in range(39-5, 68+2): #TM1
        ADRB1_MELGA_BW_dic[i]='TM1.'+str(30-5+i-(39-5))

    for i in range(75-5, 104+2): #TM2
        ADRB1_MELGA_BW_dic[i]='TM2.'+str(38-5+i-(75-5))

    for i in range(111-5, 143+2): #TM3
        ADRB1_MELGA_BW_dic[i]='TM3.'+str(22-5+i-(111-5))

    for i in range(156-5, 178+2): #TM4
        ADRB1_MELGA_BW_dic[i]='TM4.'+str(40-5+i-(156-5))

    for i in range(204-5, 229+2): #TM5
        ADRB1_MELGA_BW_dic[i]='TM5.'+str(35-5+i-(204-5))

    for i in range(285-5, 315+2): #TM6
        ADRB1_MELGA_BW_dic[i]='TM6.'+str(30-5+i-(285-5))

    for i in range(323-5, 358+2): #TM7
        ADRB1_MELGA_BW_dic[i]='TM7.'+str(33-5+i-(323-5))

    assert ADRB1_MELGA_BW_dic[59] == 'TM1.50'
    assert ADRB1_MELGA_BW_dic[87] == 'TM2.50'
    assert ADRB1_MELGA_BW_dic[139] == 'TM3.50'
    assert ADRB1_MELGA_BW_dic[166] == 'TM4.50'
    assert ADRB1_MELGA_BW_dic[219] == 'TM5.50'
    assert ADRB1_MELGA_BW_dic[305] == 'TM6.50'
    assert ADRB1_MELGA_BW_dic[340] == 'TM7.50'

    return ADRB1_MELGA_BW_dic
