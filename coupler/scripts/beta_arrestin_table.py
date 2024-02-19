#!/usr/bin/env python
"""
beta arrestin table
Based on https://gpcrdb.org/signprot/family/200_000_001/

"""
def convert_Residue2TMs():
    arrestin_BW_dic = {}
    for i in range(9, 14):
        arrestin_BW_dic[i] = 'N.S1.'+"%02d"%(i-8)
    for i in range(14,19):
        arrestin_BW_dic[i] = 'N.s1s2.'+"%02d"%(i-13)
    
    #arrestin_BW_dic[16] = 'N.S1S2.03'
    #arrestin_BW_dic[17] = 'N.S1S2.04'
    #arrestin_BW_dic[18] = 'N.S1S2.05'
    for i in range(19, 24): 
        arrestin_BW_dic[i] = 'N.S2.'+"%02d"%(i-18)

    #arrestin_BW_dic[19] = 'N.S2.01'
    #arrestin_BW_dic[20] = 'N.S2.02'
    #arrestin_BW_dic[21] = 'N.S2.03'
    #arrestin_BW_dic[22] = 'N.S2.04'
    #arrestin_BW_dic[23] = 'N.S2.05'
    for i in range(24, 27):
        arrestin_BW_dic[i] = 'N.s2s3.'+"%02d"%(i-23)

    #arrestin_BW_dic[24] = 'N.S2S3.01'
    #arrestin_BW_dic[25] = 'N.S2S3.02'
    #arrestin_BW_dic[26] = 'N.S2S3.03'
    for i in range(27, 31):
        arrestin_BW_dic[i] = 'N.S3.'+"%02d"%(i-26)

    #arrestin_BW_dic[27] = 'N.S3.01'
    #arrestin_BW_dic[28] = 'N.S3.02'
    #arrestin_BW_dic[29] = 'N.S3.03'
    #arrestin_BW_dic[30] = 'N.S3.04'
    for i in range(31, 38):
        arrestin_BW_dic[i] = 'N.s3s4.'+"%02d"%(i-30)

    #arrestin_BW_dic[31] = 'N.S3S4.01'
    #arrestin_BW_dic[32] = 'N.S3S4.02'
    #arrestin_BW_dic[33] = 'N.S3S4.03'
    #arrestin_BW_dic[34] = 'N.S3S4.04'
    #arrestin_BW_dic[35] = 'N.S3S4.05'
    #arrestin_BW_dic[36] = 'N.S3S4.06'
    #arrestin_BW_dic[37] = 'N.S3S4.07'
    for i in range(38,45):
        arrestin_BW_dic[i] = 'N.S4.'+"%02d"%(i-37)

    #arrestin_BW_dic[38] = 'N.S4.01'
    #arrestin_BW_dic[39] = 'N.S4.02'
    #arrestin_BW_dic[40] = 'N.S4.03'
    #arrestin_BW_dic[41] = 'N.S4.04'
    #arrestin_BW_dic[42] = 'N.S4.05'
    #arrestin_BW_dic[43] = 'N.S4.06'
    #arrestin_BW_dic[44] = 'N.S4.07'
    for i in range(45, 53):
        arrestin_BW_dic[i] = 'N.s4s5.'+"%02d"%(i-44)

    #arrestin_BW_dic[45] = 'N.S4S5.01'
    #arrestin_BW_dic[46] = 'N.S4S5.02'
    #arrestin_BW_dic[47] = 'N.S4S5.03'
    #arrestin_BW_dic[48] = 'N.S4S5.04'
    #arrestin_BW_dic[49] = 'N.S4S5.05'
    #arrestin_BW_dic[50] = 'N.S4S5.06'
    #arrestin_BW_dic[51] = 'N.S4S5.07'
    #arrestin_BW_dic[52] = 'N.S4S5.08'
    for i in range(53, 65):
        arrestin_BW_dic[i] = 'N.S5.'+"%02d"%(i-52)

    #arrestin_BW_dic[53] = 'N.S5.01'
    #arrestin_BW_dic[54] = 'N.S5.02'
    #arrestin_BW_dic[55] = 'N.S5.03'
    #arrestin_BW_dic[56] = 'N.S5.04'
    #arrestin_BW_dic[57] = 'N.S5.05'
    #arrestin_BW_dic[58] = 'N.S5.06'
    #arrestin_BW_dic[59] = 'N.S5.07'
    #arrestin_BW_dic[60] = 'N.S5.08'
    #arrestin_BW_dic[61] = 'N.S5.09'
    #arrestin_BW_dic[62] = 'N.S5.10'
    #arrestin_BW_dic[63] = 'N.S5.11'
    #arrestin_BW_dic[64] = 'N.S5.12'
    for in range(65, 76):
        arrestin_BW_dic[i] = 'N.s5s6.'+"%02d"%(i-64)
 
    #arrestin_BW_dic[65] = 'N.S5S6.01'
    #arrestin_BW_dic[66] = 'N.S5S6.02'
    #arrestin_BW_dic[67] = 'N.S5S6.03'
    #arrestin_BW_dic[68] = 'N.S5S6.04'
    #arrestin_BW_dic[69] = 'N.S5S6.05'
    #arrestin_BW_dic[70] = 'N.S5S6.06'
    #arrestin_BW_dic[71] = 'N.S5S6.07'
    #arrestin_BW_dic[72] = 'N.S5S6.08'
    #arrestin_BW_dic[73] = 'N.S5S6.09'
    #arrestin_BW_dic[74] = 'N.S5S6.10'
    #arrestin_BW_dic[75] = 'N.S5S6.11'
    for i in range(76, 89):
        arrestin_BW_dic[i] = 'N.S6.'+"%02d"%(i-75)
 
    #arrestin_BW_dic[76] = 'N.S6.01'
    #arrestin_BW_dic[88] = 'N.S6.13'
    for i in range(89, 100):
        arrestin_BW_dic[i] = 'N.s6h1.'+"%02d"%(i-88)
    for i in range(100, 110):
        arrestin_BW_dic[i] = 'N.H1.'+"%02d"%(i-99)
    for i in range(110, 113):
        arrestin_BW_dic[i] = 'N.h1s7.'+"%02d"%(i-109)
    for i in range(113, 119):
        arrestin_BW_dic[i] = 'N.S7.'+"%02d"%(i-112)
    for i in range(119, 128):
        arrestin_BW_dic[i] = 'N.s7s8.'+"%02d"%(i-118)
    for i in range(128, 131):
        arrestin_BW_dic[i] = 'N.S8.'+"%02d"%(i-127)
    for i in range(131, 141):
        arrestin_BW_dic[i] = 'N.s8s9.'+"%02d"%(i-130)
    for i in range(141, 153):
        arrestin_BW_dic[i] = 'N.S9.'+"%02d"%(i-140)
    for i in range(153, 164):
        arrestin_BW_dic[i] = 'N.s9s10.'+"%02d"%(i-152)
    for i in range(164, 175):
        arrestin_BW_dic[i] = 'N.S10.'+"%02d"%(i-163)
    for i in range(175, 185):
        arrestin_BW_dic[i] = 'C.s10s11.'+"%02d"%(i-174)
    for i in range(185, 191):
        arrestin_BW_dic[i] = 'C.S11.'+"%02d"%(i-184)
    for i in range(191, 198):
        arrestin_BW_dic[i] = 'C.s11s12.'+"%02d"%(i-190)
    for i in range(198, 205):
        arrestin_BW_dic[i] = 'C.S12.'+"%02d"%(i-197)
    for i in range(205, 208):
        arrestin_BW_dic[i] = 'C.s12s13.'+"%02d"%(i-204)
    for i in range(208, 211):
        arrestin_BW_dic[i] = 'C.S13.'+"%02d"%(i-207)
    for i in range(211, 215):
        arrestin_BW_dic[i] = 'C.s13s14.'+"%02d"%(i-210)
    for i in range(215, 224):
        arrestin_BW_dic[i] = 'C.S14.'+"%02d"%(i-214)
    for i in range(224, 229):
        arrestin_BW_dic[i] = 'C.s14s15.'+"%02d"%(i-223)
    for i in range(229, 243):
        arrestin_BW_dic[i] = 'C.S15.'+"%02d"%(i-228)
    for i in range(243, 248):
        arrestin_BW_dic[i] = 'C.s15s16.'+"%02d"%(i-242)
    for i in range(248, 260):
        arrestin_BW_dic[i] = 'C.S16.'+"%02d"%(i-247)
    for i in range(260, 267):
        arrestin_BW_dic[i] = 'C.s16s17.'+"%02d"%(i-259)
    for i in range(267, 276):
        arrestin_BW_dic[i] = 'C.S17.'+"%02d"%(i-266)
    for i in range(276, 318):
        arrestin_BW_dic[i] = 'C.s17s18.'+"%02d"%(i-275)
    for i in range(318, 331):
        arrestin_BW_dic[i] = 'C.S18.'+"%02d"%(i-317)
    for i in range(331, 342):
        arrestin_BW_dic[i] = 'C.s18s19.'+"%02d"%(i-330)
    for i in range(342, 355):
        arrestin_BW_dic[i] = 'C.S19.'+"%02d"%(i-341)


    return arrestin_BW_dic
