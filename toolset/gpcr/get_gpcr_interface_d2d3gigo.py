from vmd import atomsel
from vmd import molecule
from sys import argv
import MDAnalysis as mda
import pickle
import multiprocessing
from multiprocessing import Process, Pool

ncpus = multiprocessing.cpu_count()

# define enviroment
import sys, os
from pathlib import Path

home = str(Path.home())
toolset_dir = home + '/repositories/pia-x/toolset'
toolset_common = toolset_dir + "/common"
toolset_gpcr = toolset_dir + "/gpcr"
sys.path.insert(0, toolset_dir)
sys.path.insert(0, toolset_common)
sys.path.insert(0, toolset_gpcr)
# end define enviroment

from common_gpcr import *


def _calc_interface_nopolyg_noil3_noil2_noil1(t):
    # sel_AR = atomsel('protein and chain A B C R')
    # sel_A = atomsel('protein and chain A B C')
    # sel_R = atomsel('protein and chain R')
    '''
    if using above seletion, there will be memeroy leaking
    I did some cutoff test. cutoff > 6 show difference in error  <0.1 Angstrom^2
    '''
    cutoff = 8
    sel_AR = atomsel('(protein and (chain A B) and (within ' + str(cutoff) + ' of chain R)) or (protein and ((not resid '+str(polyg_vmd)+') and (not resid '+str(il3_vmd)+') and (not resid '+str(il2_vmd)+') and (not resid '+str(il1_vmd)+')) and (chain R) and (within ' + str(cutoff) + ' of chain A B))')
    sel_A =  atomsel('protein and (chain A B) and (within ' + str(cutoff) + ' of chain R)')
    sel_R =  atomsel('protein and ((not resid '+str(polyg_vmd)+') and (not resid '+str(il3_vmd)+') and (not resid '+str(il2_vmd)+') and (not resid '+str(il1_vmd)+')) and (chain R) and (within ' + str(cutoff) + ' of chain A B)')
    sel_AR.frame = t
    sel_AR.update()
    sel_A.frame = t
    sel_A.update()
    sel_R.frame = t
    sel_R.update()
    sasa_AR = sel_AR.sasa(srad=1.4, points=False)
    sasa_A = sel_A.sasa(srad=1.4, points=False)
    sasa_R = sel_R.sasa(srad=1.4, points=False)
    _interface = round((sasa_A + sasa_R - sasa_AR) / 2, 0)
    return _interface

def get_interface_nopolyg_noil3_noil2_noil1(pdbin, dcdin, Verbose=False):
    molid = molecule.load("pdb", pdbin)
    if (True):
        # multiple core version
        molecule.read(molid, "dcd", dcdin, stride=1, waitfor=-1)
        # number of frames (PDB+DCD)
        nframes = molecule.numframes(molid)
        # nframes = 100
        _ls = []
        [_ls.append(i) for i in range(1, nframes)]
        # pool = Pool(ncpus)
        pool = Pool(16)
        sasa = pool.map(_calc_interface_nopolyg_noil3_noil2_noil1, _ls)
        pool.close()
        pool.join()

    return sasa


def _calc_interface_nopolyg_noil3_noil2(t):
    # sel_AR = atomsel('protein and chain A B C R')
    # sel_A = atomsel('protein and chain A B C')
    # sel_R = atomsel('protein and chain R')
    '''
    if using above seletion, there will be memeroy leaking
    I did some cutoff test. cutoff > 6 show difference in error  <0.1 Angstrom^2
    '''
    cutoff = 8
    sel_AR = atomsel('(protein and (chain A B) and (within ' + str(cutoff) + ' of chain R)) or (protein and ((not resid '+str(polyg_vmd)+') and (not resid '+str(il3_vmd)+') and (not resid '+str(il2_vmd)+')) and (chain R) and (within ' + str(cutoff) + ' of chain A B))')
    sel_A =  atomsel('protein and (chain A B) and (within ' + str(cutoff) + ' of chain R)')
    sel_R =  atomsel('protein and ((not resid '+str(polyg_vmd)+') and (not resid '+str(il3_vmd)+') and (not resid '+str(il2_vmd)+')) and (chain R) and (within ' + str(cutoff) + ' of chain A B)')
    sel_AR.frame = t
    sel_AR.update()
    sel_A.frame = t
    sel_A.update()
    sel_R.frame = t
    sel_R.update()
    sasa_AR = sel_AR.sasa(srad=1.4, points=False)
    sasa_A = sel_A.sasa(srad=1.4, points=False)
    sasa_R = sel_R.sasa(srad=1.4, points=False)
    _interface = round((sasa_A + sasa_R - sasa_AR) / 2, 0)
    return _interface

def get_interface_nopolyg_noil3_noil2(pdbin, dcdin, Verbose=False):
    molid = molecule.load("pdb", pdbin)
    if (True):
        # multiple core version
        molecule.read(molid, "dcd", dcdin, stride=1, waitfor=-1)
        # number of frames (PDB+DCD)
        nframes = molecule.numframes(molid)
        # nframes = 100
        _ls = []
        [_ls.append(i) for i in range(1, nframes)]
        # pool = Pool(ncpus)
        pool = Pool(16)
        sasa = pool.map(_calc_interface_nopolyg_noil3_noil2, _ls)
        pool.close()
        pool.join()

    return sasa




def _calc_interface_nopolyg_noil3(t):
    # sel_AR = atomsel('protein and chain A B C R')
    # sel_A = atomsel('protein and chain A B C')
    # sel_R = atomsel('protein and chain R')
    '''
    if using above seletion, there will be memeroy leaking
    I did some cutoff test. cutoff > 6 show difference in error  <0.1 Angstrom^2
    '''
    cutoff = 8
    sel_AR = atomsel('(protein and (chain A B) and (within ' + str(cutoff) + ' of chain R)) or (protein and ((not resid '+str(polyg_vmd)+') and (not resid '+str(il3_vmd)+')) and (chain R) and (within ' + str(cutoff) + ' of chain A B))')
    sel_A = atomsel('protein and (chain A B) and (within ' + str(cutoff) + ' of chain R)')
    sel_R = atomsel('protein and ((not resid '+str(polyg_vmd)+') and (not resid '+str(il3_vmd)+')) and (chain R) and (within ' + str(cutoff) + ' of chain A B)')
    sel_AR.frame = t
    sel_AR.update()
    sel_A.frame = t
    sel_A.update()
    sel_R.frame = t
    sel_R.update()
    sasa_AR = sel_AR.sasa(srad=1.4, points=False)
    sasa_A = sel_A.sasa(srad=1.4, points=False)
    sasa_R = sel_R.sasa(srad=1.4, points=False)
    _interface = round((sasa_A + sasa_R - sasa_AR) / 2, 0)
    return _interface

def get_interface_nopolyg_noil3(pdbin, dcdin, Verbose=False):
    molid = molecule.load("pdb", pdbin)
    if (True):
        # multiple core version
        molecule.read(molid, "dcd", dcdin, stride=1, waitfor=-1)
        # number of frames (PDB+DCD)
        nframes = molecule.numframes(molid)
        # nframes = 100
        _ls = []
        [_ls.append(i) for i in range(1, nframes)]
        # pool = Pool(ncpus)
        pool = Pool(16)
        sasa = pool.map(_calc_interface_nopolyg_noil3, _ls)
        pool.close()
        pool.join()

    return sasa


def _calc_interface_noployg(t):
    '''
    if using above seletion, there will be memeroy leaking
    I did some cutoff test. cutoff > 6 show difference in error  <0.1 Angstrom^2
    '''
    cutoff = 8
    sel_AR = atomsel('(protein and (chain A B) and (within ' + str(cutoff) + ' of chain R)) or (protein and (not resid '+str(polyg_vmd)+') and (chain R) and (within ' + str(cutoff) + ' of chain A B))')
    sel_A = atomsel('protein and (chain A B) and (within ' + str(cutoff) + ' of chain R)')
    sel_R = atomsel('protein and (not resid '+str(polyg_vmd)+') and (chain R) and (within ' + str(cutoff) + ' of chain A B)')
    sel_AR.frame = t
    sel_AR.update()
    sel_A.frame = t
    sel_A.update()
    sel_R.frame = t
    sel_R.update()
    sasa_AR = sel_AR.sasa(srad=1.4, points=False)
    sasa_A = sel_A.sasa(srad=1.4, points=False)
    sasa_R = sel_R.sasa(srad=1.4, points=False)
    _interface = round((sasa_A + sasa_R - sasa_AR) / 2, 0)
    return _interface

def get_interface_nopolyg(pdbin, dcdin, Verbose=False):
    molid = molecule.load("pdb", pdbin)
    if (True):
        # multiple core version
        molecule.read(molid, "dcd", dcdin, stride=1, waitfor=-1)
        # number of frames (PDB+DCD)
        nframes = molecule.numframes(molid)
        # nframes = 100
        _ls = []
        [_ls.append(i) for i in range(1, nframes)]
        # pool = Pool(ncpus)
        pool = Pool(16)
        sasa = pool.map(_calc_interface_noployg, _ls)
        pool.close()
        pool.join()

    return sasa



if __name__ == '__main__':

    # argvs
    pdbin = argv[1]
    dcdin = argv[2]
    gpcrin = argv[3]

    # check output
    output_dir = "output"
    check_output_dir(output_dir)

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


    print("interface SASA R/AB (no ploy-G)")
    _interface = get_interface_nopolyg(pdbin, dcdin)
    with open(output_dir + '/gpcr_interface_nopolyg.p', 'wb') as fp:
        pickle.dump(_interface, fp)

    print("interface SASA R/AB (no ploy-G and no il3)")
    _interface = get_interface_nopolyg_noil3(pdbin, dcdin)
    with open(output_dir + '/gpcr_interface_nopolyg_noil3.p', 'wb') as fp:
        pickle.dump(_interface, fp)


    print("interface SASA R/AB (no ploy-G and no il3 and no il2)")
    _interface = get_interface_nopolyg_noil3_noil2(pdbin, dcdin)
    with open(output_dir + '/gpcr_interface_nopolyg_noil3_noil2.p', 'wb') as fp:
        pickle.dump(_interface, fp)

    print("interface SASA R/AB (no ploy-G and no il3 and no il2 and no il1)")
    _interface = get_interface_nopolyg_noil3_noil2(pdbin, dcdin)
    with open(output_dir + '/gpcr_interface_nopolyg_noil3_noil2_noil1.p', 'wb') as fp:
        pickle.dump(_interface, fp)
