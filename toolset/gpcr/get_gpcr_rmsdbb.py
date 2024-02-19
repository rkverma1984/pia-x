from sys import argv
from vmd import atomsel
from vmd import molecule
import pandas as pd
import numpy as np

# define enviroment
import sys, os
from pathlib import Path
home = str(Path.home())
toolset_dir = home+'/repositories/pia-x/toolset'
toolset_common = toolset_dir+"/common"
toolset_gpcr = toolset_dir+"/gpcr"
sys.path.insert(0, toolset_dir)
sys.path.insert(0, toolset_common)
sys.path.insert(0, toolset_gpcr)

from common_gpcr import *
from common_vmd import *

def get_gpcr_RMSD(alignment_str):

    _list_prot=[]
    _list_prot_R=[]
    _list_prot_A=[]
    _list_lignad=[]
    # align to the last frame of ligand
    _list_lignad_last=[]

    for t in range(1,nframes):

        sel_allatom.frame=t
        sel_allatom.update()

        sel_prot_sel_backbone.frame=t
        sel_prot_sel_backbone.update()
        
        sel_prot_R_tm_backbone.frame=t
        sel_prot_R_tm_backbone.update()
        
        sel_prot_A_sel_backbone.frame=t
        sel_prot_A_sel_backbone.update()

        sel_OBS_backbone.frame=t
        sel_OBS_backbone.update()

        sel_ligand.frame=t
        sel_ligand.update()
        
        sel_backbone.frame=t
        sel_backbone.update()

        if alignment_str==arr_align[0]:
            M=sel_prot_sel_backbone.fit(selection=ref_prot_sel_backbone)
        elif alignment_str==arr_align[1]:
            M=sel_prot_A_sel_backbone.fit(selection=ref_prot_A_sel_backbone)
        elif alignment_str==arr_align[2]:
            M=sel_prot_R_tm_backbone.fit(selection=ref_prot_R_tm_backbone)
        elif alignment_str==arr_align[3]:
            M=sel_OBS_backbone.fit(selection=ref_OBS_backbone)
        else:
            print("error")

        sel_allatom.move(M)

        _list_prot.append(sel_prot_sel_backbone.rmsd(ref_prot_sel_backbone))
        _list_prot_R.append(sel_prot_R_tm_backbone.rmsd(ref_prot_R_tm_backbone))
        _list_prot_A.append(sel_prot_A_sel_backbone.rmsd(ref_prot_A_sel_backbone))
        _list_lignad.append(sel_ligand.rmsd(ref_ligand))
        _list_lignad_last.append(sel_ligand.rmsd(ref_ligand_last))
        _df = pd.DataFrame({"all":    _list_prot,
                            "prot_R": _list_prot_R,
                            "prot_A": _list_prot_A, 
                            "lig":    _list_lignad,
                            "lig_last": _list_lignad_last})
    return _df

if __name__ == '__main__':

    # argvs
    pdbin = argv[1]
    dcdin = argv[2]
    gpcrin = argv[3]

    # check output
    output_dir="output"
    check_output_dir(output_dir)

    if gpcrin == "d3": bwtable = convert_Residue2TMs(gpcrin="DRD3")
    elif gpcrin == "d2": bwtable = convert_Residue2TMs(gpcrin="DRD2")
    elif gpcrin == "d1": bwtable = convert_Residue2TMs(gpcrin="D1R")
    elif gpcrin == "d4": bwtable = convert_Residue2TMs(gpcrin="D4R")
    elif gpcrin == "b2": bwtable = convert_Residue2TMs(gpcrin="B2")
    elif gpcrin == "m2": bwtable = convert_Residue2TMs(gpcrin="M2")
    else: print("grcp_input needs")

    # load pdb/dcd
    molid = molecule.load("pdb", pdbin)
    molecule.read(molid, "dcd", dcdin, stride=1, waitfor=-1)
    nframes = molecule.numframes(molid)  # number of frames (PDB+DCD)


    rec_df = prep_receptor_df(molid, bwtable)
    # OBS_df = prep_bindingsite_df(OBS, rec_df)
    # SBP_df = prep_bindingsite_df(SBP, rec_df)
    # OBS_resid = list(OBS_df['resid'])
    # OBS_df = prep_bindingsite_df(OBS, rec_df)
    # SBP_df = prep_bindingsite_df(SBP, rec_df)
    # OBS_resid = list(OBS_df['resid'])

    #print(OBS)
    # get OBS resid
    OBS_resid = []
    for i in OBS:
        if "EL2.52" not in i:
            #print(get_key(bwtable, "TM"+str(i)))
            OBS_resid.append(get_key(bwtable, "TM"+str(i)))

    sel_allatom, sel_backbone, ref_allatom, ref_backbone, \
    sel_prot_sel, sel_prot_sel_backbone, ref_prot_sel, ref_prot_sel_backbone, \
    sel_prot_R_tm, sel_prot_R_tm_backbone, ref_prot_R_tm, ref_prot_R_tm_backbone, \
    sel_prot_A_sel, sel_prot_A_sel_backbone, ref_prot_A_sel, ref_prot_A_sel_backbone, \
    sel_OBS, sel_OBS_backbone, ref_OBS, ref_OBS_backbone, \
    sel_ligand, ref_ligand, ref_ligand_last = \
        prep_vmd_atomsel(vmd_prot_R, vmd_prot_A, vmd_prot_B, vmd_prot_C,
                         vmd_ignore_prot_A_sel,
                         vmd_add_backbone, vmd_add_heavyatoms, vmd_add_CA,
                         vmd_add_rec_tm, OBS_resid, vmd_ligand, nframes)

    arr_align=["all_backbone","prot_A_backbone","prot_R_backbone","OBS_backbone"]

    for i in range(len(arr_align)):
        #print("RMSDBB analyses: alignmnet by "+arr_align[i])
        outputname="rmsdbb_alin_"+arr_align[i]+".hdf"
        #print(outputname)
        _df=get_gpcr_RMSD(arr_align[i])
        _df.to_hdf(output_dir+'/'+outputname,'df', mode='w')
