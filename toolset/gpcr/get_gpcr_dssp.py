## define enviroment
import sys
from pathlib import Path
home = str(Path.home())
toolset_dir = home+'/repositories/pia-x/toolset'
toolset_common = toolset_dir+"/common"
toolset_gpcr = toolset_dir+"/gpcr"
sys.path.insert(0, toolset_dir)
sys.path.insert(0, toolset_common)
sys.path.insert(0, toolset_gpcr)
from check import *
from common_gpcr import *

from sys import argv
import warnings
warnings.filterwarnings('ignore')

import mdtraj as md

def get_dssp(traj, table, resid_str, resid_end):

    newtable = table[table['name'] == "CA"]
    newtable = newtable.reset_index()
    newtable = newtable[newtable['segmentID'] == "ROR"]
    del newtable['index']
    newtable = newtable[newtable['resSeq'] >= resid_str]
    newtable = newtable[newtable['resSeq'] <= resid_end]
    sel_index = list(newtable.index)
    del newtable['serial']
    del newtable['element']
    del newtable['name']
    del newtable['chainID']
    del newtable['segmentID']
    newtable = newtable.reset_index()
    del newtable['index']

    # generate labels for final output
    ls_labels = []
    for i in range(len(newtable)):
        ls_labels.append("dssp_"+str(newtable['resSeq'][i]) + "_" + str(newtable['resName'][i]))
    # ls_labels

    # dssp calculation
    dssp = md.compute_dssp(traj)
    #dssp = md.compute_dssp(traj,simplified=False)

    for i in range(dssp.shape[0]):
        _colname = str(i + 1)
        newtable[_colname] = list(dssp[i][sel_index])

    del newtable['resSeq']
    del newtable['resName']
    # transpose newtable
    newtable = newtable.T
    # rename colnames
    newtable.columns = ls_labels

    return newtable

if __name__ == '__main__':
    # argv
    # load pdb/dcd
    PDBin = argv[1]
    DCDin = argv[2]
    gpcrin = argv[3]
    traj = md.load(DCDin, top=PDBin)
    topology = traj.topology
    table, bonds = topology.to_dataframe()
    # check output
    output_dir = "output"
    check_output_dir(output_dir)

    # check gpcrin
    if gpcrin == "d3": bwtable = convert_Residue2TMs(gpcrin="DRD3")
    elif gpcrin == "d2": bwtable = convert_Residue2TMs(gpcrin="DRD2")
    elif gpcrin == "d2short": bwtable = convert_Residue2TMs(gpcrin="D2R_short")
    elif gpcrin == "d1": bwtable = convert_Residue2TMs(gpcrin="D1R")
    elif gpcrin == "d4": bwtable = convert_Residue2TMs(gpcrin="D4R")
    elif gpcrin == "b2": bwtable = convert_Residue2TMs(gpcrin="B2")
    elif gpcrin == "m2": bwtable = convert_Residue2TMs(gpcrin="M2")
    else: print("gpcrin needs")


    resid_str=il2_str
    resid_end=il2_end
    _df=get_dssp(traj,table,resid_str,resid_end)
    _df.to_hdf(output_dir+'/df_dssp_il2.hdf', key='df', mode='w')

    _df = _df.replace(to_replace="H", value="1")
    _df = _df.replace(to_replace="C", value="0")
    _df.to_hdf(output_dir+'/df_dssp_il2_binary.hdf', key='df', mode='w')
