import MDAnalysis as mda
from MDAnalysis.analysis import align
# from MDAnalysis.analysis.rms import rmsd
# from MDAnalysis.analysis import dihedrals
from sys import argv
import pickle
import warnings

warnings.filterwarnings('ignore')
from statistics import mean
import os
# define enviroment
import sys
from pathlib import Path

home = str(Path.home())
toolset_dir = home + '/repositories/pia-x/toolset'
toolset_common = toolset_dir + "/common"
toolset_gpcr = toolset_dir + "/gpcr"
sys.path.insert(0, toolset_dir)
sys.path.insert(0, toolset_common)
sys.path.insert(0, toolset_gpcr)
from common_gpcr import *

if __name__ == '__main__':

    # check output
    output_dir = "output"
    check_output_dir(output_dir)

    pdbin = argv[1]
    dcdin = argv[2]
    # pdbin="/home/khlee/desmond/output/d2gi/qnp/d2gi_qnp.f8/mda_protein.pdb"
    # dcdin="/home/khlee/desmond/output/d2gi/qnp/d2gi_qnp.f8/d2gi_qnp.f8_e0-1800.protein.wrapped.dcd"
    u = mda.Universe(pdbin, dcdin)  # trajectory

    refpdbin = argv[3]
    # refpdbin="/home/khlee/repositories/dopamine/pj_d2d3gigo/ref_pdb_quinpirole/d2gi_qnp_f6_ref_mda.pdb"
    ref = mda.Universe(refpdbin)  # reference

    gpcrin = argv[4]
    if gpcrin == "d3":
        bwtable = convert_Residue2TMs(gpcrin="DRD3")
    elif gpcrin == "d2":
        bwtable = convert_Residue2TMs(gpcrin="DRD2")
    elif gpcrin == "d1":
        bwtable = convert_Residue2TMs(gpcrin="D1R")
    elif gpcrin == "d4":
        bwtable = convert_Residue2TMs(gpcrin="D4R")
    elif gpcrin == "b2":
        bwtable = convert_Residue2TMs(gpcrin="B2")
    elif gpcrin == "m2":
        bwtable = convert_Residue2TMs(gpcrin="M2")
    else:
        print("gpcrin needs")

    # d2tm_sel_ori = "name CA and resid 38:61,68:96,104:137,149:172,187:214,368:398,405:427 and segid PROR"
    # d3tm_sel_ori = "name CA and resid 33:56,63:91,100:133,147:170,186:213,324:354,362:384 and segid PROR"
    # d2tm_sel = "name CA and resid 70:94,109:133,190:214,388:399,404:423 and segid PROR"
    # d3tm_sel = "name CA and resid 65:89,105:129,189:213,344:355,361:380 and segid PROR"
    tm_backbone = "backbone and segid PROR and resid "
    bwidx_str = "1.36:1.59,2.38:2.66,3.22:3.55,4.39:4.62,5.36:5.63,6.30:6.60,7.32:7.54"
    strin = bwidx_str.replace(":", ",").split(",")
    strin_resi = []
    for istr in strin:
        strin_resi.append(get_key(bwtable, "TM" + istr))
    for i in range(len(strin_resi)):
        bwidx_str = bwidx_str.replace(str(strin[i]), str(strin_resi[i]))
    sel_tm_backbone = tm_backbone + bwidx_str

    # alignment
    align.AlignTraj(u,  # trajectory to align
                    ref,  # reference
                    select=sel_tm_backbone,  # selection of atoms to align
                    filename='_aligned.dcd',  # file to write the trajectory to
                    match_atoms=True,  # whether to match atoms based on mass
                    ).run()
    u_aligned = mda.Universe(pdbin, "_aligned.dcd")

    _ls_pos_z = []
    for ts in u_aligned.trajectory:  # iterate through all frames
        protein = u_aligned.select_atoms('resname qnp and not name pseu and not type H')
        _pos_z = []
        for i in range(len(protein)):
            _pos_z.append(protein[i].position[2])
        _ls_pos_z.append(mean(_pos_z))

    with open(output_dir + '/ligand_z.p', 'wb') as fp:
        pickle.dump(_ls_pos_z, fp)

    os.remove("_aligned.dcd")
