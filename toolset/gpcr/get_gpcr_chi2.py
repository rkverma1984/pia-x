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
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
from sys import argv
import pickle
import warnings
warnings.filterwarnings('ignore')


def get_gpcr_rec_Janin(u, resid, chi="chi2", Verbose=False):
    if Verbose: print(str(resid))
    r = u.select_atoms('protein and segid PROR and resid '+str(resid))
    R = mda.analysis.dihedrals.Janin(r).run()
    arr_angles = R.angles.flatten()
    _n = len(arr_angles)
    assert _n%2 == 0
    n = int(_n/2)
    mda_chi = []

    if chi=="chi2":
        for i in range(n):
            _A = arr_angles[i*2+1]
            if _A > 180:
                _A = _A-360
            mda_chi.append(_A)
        if Verbose:
            print(mda_chi)
    if chi=="chi1":
        for i in range(n):
            _A = arr_angles[i*2]
            if _A > 180:
                _A = _A-360
            mda_chi.append(_A)
        if Verbose:
            print(mda_chi)
    return mda_chi




if __name__ == '__main__':

    # argv
    pdbin = argv[1]
    dcdin = argv[2]
    bwidxin = argv[3]
    gpcrin = argv[4]

    # check
    output_dir = "output"
    check_output_dir(output_dir)

    # check gpcrin
    bwtable = convert_Residue2TMs(gpcrin)



    u = mda.Universe(pdbin, dcdin)
    resid = get_key(bwtable, "TM"+bwidxin)
    mda_chi2 = get_gpcr_rec_Janin(u,resid,chi="chi2")

    # save to files
    with open(output_dir + '/gpcr_chi2_' + str(bwidxin) + '.p', 'wb') as fp:
        pickle.dump(mda_chi2, fp)