
# those are for D3R...
# receptor
# rec_chainid='R'
rec_segments=['27:57', '63:91', '100:134', '145:170', '186:218', '322:354', '362:400']
# rec_segments2=['63:91', '100:134', '145:170', '186:218']
rec_labels=['TM1','TM2','TM3','TM4','TM5','TM6','TM7']

#il1
il1_vmd="58 to 62"
# il1_mda=['58:62']
#il2
il2_vmd="135 to 144"
# il2_mda=['135:144']
# il3 (QPRGVPLR)
il3_vmd="316 to 323"
# il3_mda=['316:323']
# poly-GLY region
polyg_vmd="229 to 237"
# polyg_mda=['229:237']


# for D2
# # rec_chainid='R'
# rec_segments=['31:62', '68:96', '104:138', '147:172', '187:219', '366:398', '405:443']
# # rec_segments2=['68:96', '104:138', '147:172', '187:219']
# rec_labels=['TM1','TM2','TM3','TM4','TM5','TM6','TM7']
#
#
# # il1
# il1_vmd = "63 to 67"
# il1_mda = ['63:67']
# #il2
# il2_vmd="139 to 146"
# il2_mda=['139:146']
# # il3 (RRKLSQQK)
# il3_vmd="360 to 367"
# il3_mda=['360:367']
# # poly-GLY region
# polyg_vmd="230 to 238"
# polyg_mda=['230:238']




# common vmd atom selection key string
vmd_prot_R = "(protein and chain R)"
vmd_prot_A = "(protein and chain A)"
vmd_ignore_prot_A_sel = " and (not (resid 1 to 8 31 to 34 308 to 320 191 to 195 215 to 220))"
vmd_prot_B = "(protein and chain B)"
vmd_prot_C = "(protein and chain C)"
vmd_all = "all"
#
vmd_add_backbone = ' and backbone'
vmd_add_heavyatoms = ' and not hydrogen'
vmd_add_CA = ' and name CA'
# tm region vmd string
vmd_add_rec_tm = " and resid "
for i in range(len(rec_segments)):
    vmd_add_rec_tm = vmd_add_rec_tm+rec_segments[i].split(":")[0]+" to "+rec_segments[i].split(":")[1]+" "
#print(vmd_rec_tm)
# ligand
#vmd_ligand = "not protein and (not resname ACE NMA) and not hydrogen and (resname PRM P12 BRO qnp QNP)"
vmd_ligand = "not hydrogen and (resname PRM P12 BRO qnp QNP)"
