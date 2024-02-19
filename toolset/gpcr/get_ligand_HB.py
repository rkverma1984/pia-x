# define enviroment
import sys
from sys import argv
from pathlib import Path
home = str(Path.home())
toolset_dir = home + '/repositories/pia-x/toolset'
toolset_common = toolset_dir + "/common"
toolset_gpcr = toolset_dir + "/gpcr"
sys.path.insert(0, toolset_dir)
sys.path.insert(0, toolset_common)
sys.path.insert(0, toolset_gpcr)
from common_gpcr import *
import MDAnalysis as mda
import pickle
import warnings
warnings.filterwarnings('ignore')

if __name__ == '__main__':

    # check output
    output_dir = "output"
    check_output_dir(output_dir)

    pdbin = argv[1]
    dcdin = argv[2]
    u = mda.Universe(pdbin, dcdin)

    sel = u.select_atoms('not protein and not (resname NMA or resname CYT) and not type C')  # select ligand (may need to change for other cases)
    # print(sel)
    # print(len(sel))
    ligand_resid = (list(set(sel.residues.resids)))
    print(ligand_resid)
    assert len(ligand_resid) == 1
    ligand_resname = (list(set(sel.residues.resnames)))
    lig_resname = ligand_resname[0]
    print(lig_resname)

    _ls_hhcontacts = []
    for ts in u.trajectory:                # iterate through all frames
        # within 3.5 A, ligand heavy atoms have contact with residue
        sel = u.select_atoms('(around 3.5 resname '+str(lig_resname)+') and not type H and not type C')
        _hhcontact = sel.residues.resids.tolist()
        _ls_hhcontacts.append(_hhcontact)
        print(_hhcontact)

    # save to files
    with open(output_dir + '/ligand_HB.p', 'wb') as fp:
        pickle.dump(_ls_hhcontacts, fp)
