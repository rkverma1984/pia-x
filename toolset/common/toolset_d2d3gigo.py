import numpy as np
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline



import seaborn as sns
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
#%matplotlib inline
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Dihedral
## https://www.mdanalysis.org/docs/documentation_pages/analysis/dihedrals.html#MDAnalysis.analysis.dihedrals.Ramachandran.angles

from matplotlib.path import Path
import matplotlib.patches as patches
## https://matplotlib.org/3.3.1/tutorials/advanced/path_tutorial.html#sphx-glr-tutorials-advanced-path-tutorial-py

import warnings
warnings.filterwarnings('ignore')


def get_phi_psi(u, resid):
    
    ## mda
    ref = u
    r = ref.select_atoms("segid PROR and resid "+str(resid))
    R = mda.analysis.dihedrals.Ramachandran(r).run()

    ## assign phi,psi into x,y
    phi,psi=R.angles.T
    return phi,psi


def plot_phi_psi(phi,psi):
    x=phi
    y=psi
    
    ## define plot range
    xedges = np.arange(-180,180,10)
    yedges = xedges
    
    ## 2D-histogram
    ## https://numpy.org/doc/stable/reference/generated/numpy.histogram2d.html
    H, xedges, yedges = np.histogram2d(x.tolist()[0],y.tolist()[0],bins=(xedges, yedges))

    ## to do...
    ## may want to convert it into free energy surface
    #pi=H/len(x.tolist()[0])
    pi=H

    ## prepare for plotting
    x = 0.5 * (xedges[:-1] + xedges[1:])
    y = 0.5 * (yedges[:-1] + yedges[1:])
    pi = pi.T
    vmin = pi.min()
    vmax = pi.max()


    ## draw a path around alpha-helix
    ## http://www.cryst.bbk.ac.uk/PPS2/course/section8/ss-960531_6.html#HEADING5
    ##
    verts = [
       (-71., -34.),  # left, bottom
       (-71., -57.),  # left, top
       (-48., -57.),  # right, top
       (-48., -34.),  # right, bottom
       (-71., -34.),  # ignored
    ]
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    path = Path(verts, codes)




    ## plotting
    nrow=1
    ncol=1
    figheight=4
    figwidth=4
    linewidth=1
    maxylim=8

    fig, (ax1) = plt.subplots(nrow, ncol,figsize=(figwidth,figheight))
    fig.subplots_adjust(hspace=1.5) # make a little extra space between the subplots


    ## counter surface heat map
    mappable = ax1.contourf(x, y, pi, vmin=vmin, vmax = vmax, cmap = plt.get_cmap('hot'), levels=None)
    ax1.set_xlabel('phi')
    ax1.set_ylabel('psi')


    ## 3-10 helix
    # http://www.cryst.bbk.ac.uk/PPS2/course/section8/ss-960531_7.html
    # https://swissmodel.expasy.org/course/text/chapter1.htm
    # https://swissmodel.expasy.org/course/text/chapter1.htm
    ## pure alpha helix
    ## http://www.cryst.bbk.ac.uk/PPS2/course/section8/ss-960531_6.html#HEADING5

    ## 
    ax1.scatter(-74.0, -4.0,marker='^',color='cyan',label='3-10 Helix' )
    ax1.text(-60.0, -4, r'3-10 Helix', {'color': 'cyan', 'fontsize': 16})

    ## 
    ax1.scatter(-57.8, -47.0,marker='s',color='cyan',label='alpha Helix' )
    ax1.text(-40.0, -47, r'pure alpha Helix', {'color': 'cyan', 'fontsize': 16})

    patch = patches.PathPatch(path, facecolor='none', lw=1)
    ax1.add_patch(patch)


    ## 
    axins = inset_axes(ax1,
                       width="5%",  # width = 5% of parent_bbox width
                       height="100%",  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1.05, 0., 1, 1),
                       bbox_transform=ax1.transAxes,
                       borderpad=0)
    cbar = fig.colorbar(mappable, cax=axins)
    cbar.set_label('counts')

    plt.show()
    return()

def plot_phi_psi_v2(fig,axin,phi,psi,titlein,resid,type):
    x=phi
    y=psi

    ## define plot range
    xedges = np.arange(-180,180,10)
    yedges = xedges
    
    ## 2D-histogram
    ## https://numpy.org/doc/stable/reference/generated/numpy.histogram2d.html
    H, xedges, yedges = np.histogram2d(x.tolist()[0],y.tolist()[0],bins=(xedges, yedges))

    ## to do...
    ## may want to convert it into free energy surface
    #pi=H/len(x.tolist()[0])
    pi=H

    ## prepare for plotting
    x = 0.5 * (xedges[:-1] + xedges[1:])
    y = 0.5 * (yedges[:-1] + yedges[1:])
    pi = pi.T
    vmin = pi.min()
    vmax = pi.max()


    ## draw a path around alpha-helix
    ## http://www.cryst.bbk.ac.uk/PPS2/course/section8/ss-960531_6.html#HEADING5
    ##
    verts = [
       (-71., -34.),  # left, bottom
       (-71., -57.),  # left, top
       (-48., -57.),  # right, top
       (-48., -34.),  # right, bottom
       (-71., -34.),  # ignored
    ]
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    path = Path(verts, codes)



    ## counter surface heat map
    mappable = axin.contourf(x, y, pi, vmin=vmin, vmax = vmax, cmap = plt.get_cmap('hot'), levels=None)
    axin.set_xlabel('phi')
    axin.set_ylabel('psi')
    axin.set_title(titlein)


    ## 3-10 helix
    # http://www.cryst.bbk.ac.uk/PPS2/course/section8/ss-960531_7.html
    # https://swissmodel.expasy.org/course/text/chapter1.htm
    # https://swissmodel.expasy.org/course/text/chapter1.htm
    ## pure alpha helix
    ## http://www.cryst.bbk.ac.uk/PPS2/course/section8/ss-960531_6.html#HEADING5

    ## 
    axin.scatter(-74.0, -4.0,marker='^',color='cyan',label='3-10 Helix' )
    axin.text(-60.0, -4, r'3-10 Helix', {'color': 'cyan', 'fontsize': 16})

    ## 
    axin.scatter(-57.8, -47.0,marker='s',color='cyan',label='alpha Helix' )
    axin.text(-40.0, -47, r'pure alpha Helix', {'color': 'cyan', 'fontsize': 16})

    # reference
    if (type=="D3"):
        x,y=get_phi_psi_7cmu(resid)
        axin.scatter(x,y,marker='s',color='green',label='7cmu' )
        axin.text(x,y, r'7cmu', {'color': 'green', 'fontsize': 16})
        x,y=get_phi_psi_7cmv(resid)
        axin.scatter(x,y,marker='s',color='green',label='7cmv' )
        axin.text(x,y, r'7cmv', {'color': 'green', 'fontsize': 16})
    elif (type=="D2"):
        x,y=get_phi_psi_6vms(resid)
        axin.scatter(x,y,marker='s',color='yellow',label='6vms' )
        axin.text(x,y, r'6vms', {'color': 'yellow', 'fontsize': 16})
    # redrence end--

    patch = patches.PathPatch(path, facecolor='none', lw=1)
    axin.add_patch(patch)


    ## 
    axins = inset_axes(axin,
                       width="5%",  # width = 5% of parent_bbox width
                       height="100%",  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1.05, 0., 1, 1),
                       bbox_transform=axin.transAxes,
                       borderpad=0)
    cbar = fig.colorbar(mappable, cax=axins)
    cbar.set_label('counts')

    #plt.show()
    return()

def wrapper_d3gi(resid):
    allphi=[]
    allpsi=[]
    PDBin ="/home/khlee/desmond/output/d3gi/pd/d3gi_pd.f1/mda_protein.pdb"
    DCDin1="/home/khlee/desmond/output/d3gi/pd/d3gi_pd.f1/d3gi_pd.f1_e0-2400.protein.wrapped.dcd"
    DCDin2="/home/khlee/desmond/output/d3gi/pd/d3gi_pd.f2/d3gi_pd.f2_e0-2400.protein.wrapped.dcd"
    DCDin3="/home/khlee/desmond/output/d3gi/pd/d3gi_pd.f3/d3gi_pd.f3_e0-2400.protein.wrapped.dcd"
    DCDin4="/home/khlee/desmond/output/d3gi/pd/d3gi_pd.f4/d3gi_pd.f4_e0-2400.protein.wrapped.dcd"
    DCDin5="/home/khlee/desmond/output/d3gi/pd/d3gi_pd.f5/d3gi_pd.f5_e0-2400.protein.wrapped.dcd"
    DCDin6="/home/khlee/desmond/output/d3gi/pd/d3gi_pd.f6/d3gi_pd.f6_e0-2400.protein.wrapped.dcd"
    u = mda.Universe(PDBin, DCDin1, DCDin2, DCDin3, DCDin4, DCDin5, DCDin6)
    _phi,_psi=get_phi_psi(u, resid)
    allphi=allphi+list(_phi[0])
    allpsi=allpsi+list(_psi[0])

    PDBin ="/home/khlee/desmond/output/d3gi/prm/d3gi_prm.f7/mda_protein.pdb"
    DCDin1="/home/khlee/desmond/output/d3gi/prm/d3gi_prm.f7/d3gi_prm.f7_e0-2400.protein.wrapped.dcd"
    DCDin2="/home/khlee/desmond/output/d3gi/prm/d3gi_prm.f8/d3gi_prm.f8_e0-2400.protein.wrapped.dcd"
    DCDin3="/home/khlee/desmond/output/d3gi/prm/d3gi_prm.f9/d3gi_prm.f9_e0-2400.protein.wrapped.dcd"
    DCDin4="/home/khlee/desmond/output/d3gi/prm/d3gi_prm.f10/d3gi_prm.f10_e0-2400.protein.wrapped.dcd"
    DCDin5="/home/khlee/desmond/output/d3gi/prm/d3gi_prm.f11/d3gi_prm.f11_e0-2400.protein.wrapped.dcd"
    u = mda.Universe(PDBin, DCDin1, DCDin2, DCDin3, DCDin4, DCDin5)
    _phi,_psi=get_phi_psi(u, resid)
    allphi=allphi+list(_phi[0])
    allpsi=allpsi+list(_psi[0])

    allphi=np.array([allphi])
    allpsi=np.array([allpsi])
    
    return allphi,allpsi

def wrapper_d3go(resid):
    allphi=[]
    allpsi=[]
    PDBin ="/home/khlee/desmond/output/d3go/pd/d3go_pd.f1/mda_protein.pdb"
    DCDin1="/home/khlee/desmond/output/d3go/pd/d3go_pd.f1/d3go_pd.f1_e0-2400.protein.wrapped.dcd"
    DCDin2="/home/khlee/desmond/output/d3go/pd/d3go_pd.f2/d3go_pd.f2_e0-2400.protein.wrapped.dcd"
    DCDin3="/home/khlee/desmond/output/d3go/pd/d3go_pd.f3/d3go_pd.f3_e0-2400.protein.wrapped.dcd"
    DCDin4="/home/khlee/desmond/output/d3go/pd/d3go_pd.f4/d3go_pd.f4_e0-2400.protein.wrapped.dcd"
    DCDin5="/home/khlee/desmond/output/d3go/pd/d3go_pd.f5/d3go_pd.f5_e0-2400.protein.wrapped.dcd"
    DCDin6="/home/khlee/desmond/output/d3go/pd/d3go_pd.f6/d3go_pd.f6_e0-2400.protein.wrapped.dcd"
    u = mda.Universe(PDBin, DCDin1, DCDin2, DCDin3, DCDin4, DCDin5, DCDin6)
    _phi,_psi=get_phi_psi(u, resid)
    allphi=allphi+list(_phi[0])
    allpsi=allpsi+list(_psi[0])

    PDBin ="/home/khlee/desmond/output/d3go/prm/d3go_prm.f1/mda_protein.pdb"
    DCDin1="/home/khlee/desmond/output/d3go/prm/d3go_prm.f1/d3go_prm.f1_e0-2400.protein.wrapped.dcd"
    DCDin2="/home/khlee/desmond/output/d3go/prm/d3go_prm.f2/d3go_prm.f2_e0-2400.protein.wrapped.dcd"
    DCDin3="/home/khlee/desmond/output/d3go/prm/d3go_prm.f3/d3go_prm.f3_e0-2400.protein.wrapped.dcd"
    DCDin4="/home/khlee/desmond/output/d3go/prm/d3go_prm.f4/d3go_prm.f4_e0-2400.protein.wrapped.dcd"
    DCDin5="/home/khlee/desmond/output/d3go/prm/d3go_prm.f5/d3go_prm.f5_e0-2400.protein.wrapped.dcd"
    DCDin6="/home/khlee/desmond/output/d3go/prm/d3go_prm.f6/d3go_prm.f6_e0-2400.protein.wrapped.dcd"
    u = mda.Universe(PDBin, DCDin1, DCDin2, DCDin3, DCDin4, DCDin5, DCDin6)
    phi,psi=get_phi_psi(u, resid)
    _phi,_psi=get_phi_psi(u, resid)
    allphi=allphi+list(_phi[0])
    allpsi=allpsi+list(_psi[0])

    allphi=np.array([allphi])
    allpsi=np.array([allpsi])
    
    return allphi,allpsi

def wrapper_d2gi(resid):
    allphi=[]
    allpsi=[]
    PDBin ="/home/khlee/desmond/output/d2gi/bro/d2gi_bro.f12/mda_protein.pdb"
    DCDin1="/home/khlee/desmond/output/d2gi/bro/d2gi_bro.f12/d2gi_bro.f12_e0-2400.protein.wrapped.dcd"
    DCDin2="/home/khlee/desmond/output/d2gi/bro/d2gi_bro.f13/d2gi_bro.f13_e0-2400.protein.wrapped.dcd"
    DCDin3="/home/khlee/desmond/output/d2gi/bro/d2gi_bro.f14/d2gi_bro.f14_e0-2400.protein.wrapped.dcd"
    DCDin4="/home/khlee/desmond/output/d2gi/bro/d2gi_bro.f15/d2gi_bro.f15_e0-2400.protein.wrapped.dcd"
    u = mda.Universe(PDBin, DCDin1, DCDin2, DCDin3, DCDin4)
    _phi,_psi=get_phi_psi(u, resid)
    allphi=allphi+list(_phi[0])
    allpsi=allpsi+list(_psi[0])

    allphi=np.array([allphi])
    allpsi=np.array([allpsi])
    
    return allphi,allpsi

def wrapper_d2go(resid):
    allphi=[]
    allpsi=[]
    PDBin ="/home/khlee/desmond/output/d2go/bro/d2go_bro.f11/mda_protein.pdb"
    DCDin1="/home/khlee/desmond/output/d2go/bro/d2go_bro.f11/d2go_bro.f11_e0-2400.protein.wrapped.dcd"
    DCDin2="/home/khlee/desmond/output/d2go/bro/d2go_bro.f12/d2go_bro.f12_e0-2400.protein.wrapped.dcd"
    DCDin3="/home/khlee/desmond/output/d2go/bro/d2go_bro.f13/d2go_bro.f13_e0-2400.protein.wrapped.dcd"
    DCDin4="/home/khlee/desmond/output/d2go/bro/d2go_bro.f14/d2go_bro.f14_e0-2400.protein.wrapped.dcd"
    DCDin5="/home/khlee/desmond/output/d2go/bro/d2go_bro.f15/d2go_bro.f15_e0-2400.protein.wrapped.dcd"
    u = mda.Universe(PDBin, DCDin1, DCDin2, DCDin3, DCDin4, DCDin5)
    _phi,_psi=get_phi_psi(u, resid)
    allphi=allphi+list(_phi[0])
    allpsi=allpsi+list(_psi[0])
    
    allphi=np.array([allphi])
    allpsi=np.array([allpsi])
    
    return allphi,allpsi





def get_phi_psi_7cmu(resid):
    PDBin ="pdb_ref/7cmu.pdb"
    u = mda.Universe(PDBin)
    ref = u
    r = ref.select_atoms("index 6634:8761 and resid "+str(resid))
    R = mda.analysis.dihedrals.Ramachandran(r).run()
    ## assign phi,psi into x,y
    phi,psi=R.angles.T
    phi=phi.tolist()[0][0]
    psi=psi.tolist()[0][0]
    return phi,psi

def get_phi_psi_7cmv(resid):
    PDBin ="pdb_ref/7cmv.pdb"
    u = mda.Universe(PDBin)
    ref = u
    r = ref.select_atoms("index 6634:8761 and resid "+str(resid))
    R = mda.analysis.dihedrals.Ramachandran(r).run()
    ## assign phi,psi into x,y
    phi,psi=R.angles.T
    phi=phi.tolist()[0][0]
    psi=psi.tolist()[0][0]
    return phi,psi

def get_phi_psi_6vms(resid):
    PDBin ="pdb_ref/6vms.pdb"
    u = mda.Universe(PDBin)
    ref = u
    r = ref.select_atoms("index 6671:8914 and resid "+str(resid))
    R = mda.analysis.dihedrals.Ramachandran(r).run()
    ## assign phi,psi into x,y
    phi,psi=R.angles.T
    phi=phi.tolist()[0][0]
    psi=psi.tolist()[0][0]
    return phi,psi