from vmd import atomsel
from vmd import molecule
from sys import argv
import pickle
import multiprocessing
from multiprocessing import Process, Pool

ncpus = multiprocessing.cpu_count()

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


def _calc_interface(t):
    # sel_AR = atomsel('protein and chain A B C R')
    # sel_A = atomsel('protein and chain A B C')
    # sel_R = atomsel('protein and chain R')
    '''
    if using above seletion, there will be memeroy leaking
    I did some cutoff test. cutoff > 6 show difference in error  <0.1 Angstrom^2
    '''
    cutoff = 8
    sel_AR = atomsel('(protein and (chain A B) and (within ' + str(
        cutoff) + ' of chain R)) or (protein and (chain R) and (within ' + str(cutoff) + ' of chain A B))')
    sel_A = atomsel('protein and (chain A B) and (within ' + str(cutoff) + ' of chain R)')
    sel_R = atomsel('protein and (chain R) and (within ' + str(cutoff) + ' of chain A B)')
    sel_AR.frame = t
    sel_AR.update()
    sel_A.frame = t
    sel_A.update()
    sel_R.frame = t
    sel_R.update()
    # sasa_AR, points = sel_AR.sasa(srad=1.4, points=True)
    # sasa_A, points = sel_A.sasa(srad=1.4, points=True)
    # sasa_R, points = sel_R.sasa(srad=1.4, points=True)
    sasa_AR = sel_AR.sasa(srad=1.4, points=False)
    sasa_A = sel_A.sasa(srad=1.4, points=False)
    sasa_R = sel_R.sasa(srad=1.4, points=False)
    _interface = round((sasa_A + sasa_R - sasa_AR) / 2, 0)
    return _interface


def get_interface(pdbin, dcdin, verbose=False, mp=True):
    molid = molecule.load("pdb", pdbin)
    sasa = []

    if not mp:
        # single core version
        molecule.read(molid, "dcd", dcdin, stride=1, waitfor=-1)
        # number of frames (PDB+DCD)
        nframes = molecule.numframes(molid)
        sel_AR = atomsel('protein and chain A B C R')
        sel_A = atomsel('protein and chain A B C')
        sel_R = atomsel('protein and chain R')
        # sasa_AR, points = sel_AR.sasa(srad=1.4, points=True)
        # sasa_A, points = sel_A.sasa(srad=1.4, points=True)
        # sasa_R, points = sel_R.sasa(srad=1.4, points=True)

        for t in range(1, nframes):
            sel_AR.frame = t
            sel_AR.update()
            sel_A.frame = t
            sel_A.update()
            sel_R.frame = t
            sel_R.update()
            sasa_AR, points = sel_AR.sasa(srad=1.4, points=True)
            sasa_A, points = sel_A.sasa(srad=1.4, points=True)
            sasa_R, points = sel_R.sasa(srad=1.4, points=True)
            area = round((sasa_A + sasa_R - sasa_AR) / 2, 2)
            sasa.append(area)
            if verbose:
                print(area)
    if mp:
        # multiple core version
        molecule.read(molid, "dcd", dcdin, stride=1, waitfor=-1)
        # number of frames (PDB+DCD)
        nframes = molecule.numframes(molid)
        # nframes = 100
        _ls = []
        [_ls.append(i) for i in range(1, nframes)]
        # pool = Pool(ncpus)
        pool = Pool(16)
        sasa = pool.map(_calc_interface, _ls)
        pool.close()
        pool.join()

    return sasa


if __name__ == '__main__':
    print("interface SASA R/AB")
    pdbin = argv[1]
    dcdin = argv[2]
    output_dir = "output"
    check_output_dir(output_dir)

    _interface = get_interface(pdbin, dcdin)
    with open(output_dir + '/gpcr_interface.p', 'wb') as fp:
        pickle.dump(_interface, fp)
