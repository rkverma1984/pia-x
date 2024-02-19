import sys
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pandas as pd
import numpy as np
from scipy.interpolate import make_interp_spline, BSpline
import pickle
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



def barplot(df, title, xlabel, ylabel):

    # figure parameters
    nrow = 1
    ncol = 1
    figheight =  4
    figwidth = 2*(round(len(df.index)/10)+1)
    linewidth = 2
    width = 0.8
    fig, ax = plt.subplots(nrow, ncol, figsize=(figwidth, figheight))
    plt.xticks(rotation=90)


    labels=list(df.index)
    n = len(df.columns) # number of column
    x  = n*np.arange(len(labels))

    for i in range(n):
        data = df[df.columns[i]].astype(float)
        ax.bar(x+i , data, width, label=df.columns[i])


    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    #ax.grid(True)
    #fig.tight_layout()
    plt.show()
    

def get_smooth(x,y):
    # 300 represents number of points to make between T.min and T.max
    xnew = np.linspace(x.min(), x.max(), 300) 
    # type: BSpline
    spl = make_interp_spline(x, y, k=1)
    power_smooth = spl(xnew)
    return(xnew, power_smooth)


def bining(datain, nbin, binstr, binend):
    """ do bining for datain

    :param datain:
    :param nbin:
    :param binstr:
    :param binend:
    :return: xbin & hist_n
    """
    # do bining for all data
    hist, bins = np.histogram(datain, bins=nbin, range=[binstr, binend],density=False)
    hist_n = hist/sum(hist)
    
    # manipulate bin
    bin_shift=(bins[1]-bins[0])/2
    nbin_center=len(bins)-1
    xbin=bins+bin_shift
    xbin=xbin[0:nbin_center]
    
    return xbin, hist_n


def bining_last_half(datain, nbin, binstr, binend):
    # do bining for last half
    _datain=list(datain)
    n0 = int(len(_datain)/2)
    _datain = _datain[n0:]

    hist, bins = np.histogram(_datain, bins=nbin, range=[binstr, binend], density=False)
    hist_n = hist / sum(hist)

    # manipulate bin
    bin_shift = (bins[1] - bins[0]) / 2
    nbin_center = len(bins) - 1
    xbin = bins + bin_shift
    xbin = xbin[0:nbin_center]

    return xbin, hist_n



def plot_evol_hist(t,datain,nbin=30,binstr=-180,binend=180,
                   title='',title_1='',title_2='',
                   xlabel_1='',ylabel_1='',xlabel_2='',ylabel_2='',data_label=''):
    #t = np.arange(len(datain))

    xbin,hist_n = bining_last_half(datain, nbin, binstr, binend)

    # smooth plot
    x,y=(xbin,hist_n)

    # figure parameters
    nrow=1
    ncol=2
    figheight=2.5
    figwidth=6
    linewidth=1

    fig, (ax1,ax2) = plt.subplots(nrow, ncol,figsize=(figwidth,figheight))
    fig.suptitle(title, fontsize="x-large")
    fig.subplots_adjust(hspace = 0.05, wspace=0.05)

    gs=gridspec.GridSpec(1, 2, width_ratios=[6, 1]) 

    ax1 = plt.subplot(gs[0])
    ax1.plot(t,datain,'bo',markersize=1,lw=linewidth,label=data_label)
    ax1.set_xlabel(xlabel_1)
    ax1.set_ylabel(ylabel_1)
    ax1.set_title(title_1)
    ax1.grid(True, linestyle='dotted', linewidth=0.2)
    #ax1.legend(loc='upper right', frameon=False)
    ax1.set_ylim(binstr,binend)
    ax1.tick_params(top=True, bottom=True, left=True, right=True,
                    labelleft=True, labelbottom=True)

    ax2 = plt.subplot(gs[1])
    ax2.plot(y,x,lw=linewidth,label=data_label)
    ax2.set_xlabel(xlabel_2)
    ax2.set_ylabel(ylabel_2)
    ax2.set_title(title_2)
    ax2.grid(True, linestyle='dotted', linewidth=0.2)
    #ax2.legend(loc='upper right', frameon=False)
    ax2.set_ylim(binstr,binend)
    ax2.tick_params(top=True, bottom=True, left=True, right=True,
                    labelleft=False, labelbottom=True)

    plt.show()
    #plt.savefig(outfilename)


def calc_np_ave(list):
    _ave=np.mean(np.array(list))
    return _ave

def calc_np_std(list):
    _std=np.std(np.array(list))
    return _std


'''Only one input data'''
def plot_evol_hist_stat(t=None,datain=None,
                        nbin=30, binstr=0, binend=8,
                        title='',
                        title_1='', title_2='', xlabel_1='', ylabel_1='',
                        xlabel_2='', ylabel_2='',
                        data_label='',
                        outfilename=None,
                        colorin="red"):

    t=t
    style="p"
    n = len(datain)
    datain=list(datain)
    #colorin="red"

    #if n > 900:
    npoint = 100
    n0 = int(n / 3)
    n1 = int(n / 3 * 2)
    n2 = int(n) - npoint*2
    n3 = int(n) - npoint
    xx = (1, 2, 3, 4)
    yy = [calc_np_ave(datain[n0:n1]), calc_np_ave(datain[n1:n2]),
              calc_np_ave(datain[n2:n3]), calc_np_ave(datain[n3:n])]
    yyerr = [calc_np_std(datain[n0:n1]), calc_np_std(datain[n1:n2]),
              calc_np_std(datain[n2:n3]), calc_np_std(datain[n3:n])]

    #t = np.arange(len(datain))
    xbin, hist_n = bining_last_half(datain, nbin, binstr, binend)

    # smooth plot
    x, y = get_smooth(xbin, hist_n)

    # figure parameters
    nrow = 1
    ncol = 3
    figheight = 2.5
    figwidth = 8
    linewidth = 1

    fig, axs = plt.subplots(nrow, ncol, figsize=(figwidth, figheight))
    fig.suptitle(title, fontsize="x-large")
    fig.subplots_adjust(hspace=0.05, wspace=0.1)

    gs = gridspec.GridSpec(nrow, ncol, width_ratios=[6, 1, 2])

    ax1 = plt.subplot(gs[0])
    if style=='l':
        ax1.plot(t, datain, markersize=1, lw=linewidth, label=data_label,color=colorin)
    elif style=='p':
        ax1.plot(t,datain,'bo',markersize=1,lw=linewidth,label=data_label,color=colorin)
    ax1.set_xlabel(xlabel_1)
    ax1.set_ylabel(ylabel_1)
    ax1.set_title(title_1)
    ax1.grid(True, linestyle='dotted', linewidth=0.2)
    ax1.legend(loc='upper left', frameon=False)
    ax1.set_ylim(binstr, binend)
    ax1.tick_params(top=True, bottom=True, left=True, right=True,
                    labelleft=True, labelbottom=True)


    ax2 = plt.subplot(gs[1])
    ax2.plot(y, x, lw=linewidth, label=data_label,color=colorin)
    ax2.set_xlabel(xlabel_2)
    ax2.set_ylabel(ylabel_2)
    ax2.set_title(title_2)
    ax2.grid(True, linestyle='dotted', linewidth=0.2)
    # ax2.legend(loc='upper right', frameon=False)
    ax2.set_ylim(binstr, binend)
    ax2.tick_params(top=True, bottom=True, left=True, right=True,
                    labelleft=False, labelbottom=True)

    ax3 = plt.subplot(gs[2])
    #ax3.errorbar(xx, yy, 0, yyerr'bo', fmt='o', lw=linewidth, label=data_label)
    ax3.errorbar(xx, yy, yyerr, 0, fmt='o',color=colorin)
    ax3.set_xlabel('')
    ax3.set_ylabel('')
    ax3.set_title('')
    ax3.grid(True, linestyle='dotted', linewidth=0.2)
    # ax3.legend(loc='upper right', frameon=False)
    ax3.set_ylim(binstr, binend)
    ax3.set_xlim(0, 5)
    ax3.tick_params(top=True, bottom=True, left=True, right=True,
                    labelleft=False, labelbottom=False)

    #plt.show()
    plt.savefig(outfilename+".png", dpi=300, bbox_inches='tight')

'''another version with two input data'''


def plot_evol_hist_stat2(t=None, datain=None, datain2=None,
                         nbin=30, binstr=0, binend=8,
                         title='',
                         title_1='', title_2='', xlabel_1='', ylabel_1='',
                         xlabel_2='', ylabel_2='',
                         data_label='', data_label2='',
                         outfilename=None,
                         colorin="red", colorin2="blue"):
    t = t
    style = "p"
    n = len(datain)
    n2 = len(datain2)
    assert n == n2
    datain = list(datain)
    datain2 = list(datain2)
    # colorin="red"

    # if n > 900:
    npoint = 100
    n0 = int(n / 3)
    n1 = int(n / 3 * 2)
    n2 = int(n) - npoint * 2
    n3 = int(n) - npoint
    xx = (1, 2, 3, 4)
    yy = [calc_np_ave(datain[n0:n1]), calc_np_ave(datain[n1:n2]),
          calc_np_ave(datain[n2:n3]), calc_np_ave(datain[n3:n])]
    yyerr = [calc_np_std(datain[n0:n1]), calc_np_std(datain[n1:n2]),
             calc_np_std(datain[n2:n3]), calc_np_std(datain[n3:n])]

    yy2 = [calc_np_ave(datain2[n0:n1]), calc_np_ave(datain2[n1:n2]),
           calc_np_ave(datain2[n2:n3]), calc_np_ave(datain2[n3:n])]
    yyerr2 = [calc_np_std(datain2[n0:n1]), calc_np_std(datain2[n1:n2]),
              calc_np_std(datain2[n2:n3]), calc_np_std(datain2[n3:n])]

    # t = np.arange(len(datain))
    xbin, hist_n = bining_last_half(datain, nbin, binstr, binend)
    xbin2, hist_n2 = bining_last_half(datain2, nbin, binstr, binend)

    # smooth plot
    x, y = get_smooth(xbin, hist_n)
    x2, y2 = get_smooth(xbin2, hist_n2)

    # figure parameters
    nrow = 1
    ncol = 3
    figheight = 2.5
    figwidth = 8
    linewidth = 1

    fig, axs = plt.subplots(nrow, ncol, figsize=(figwidth, figheight))
    fig.suptitle(title, fontsize="x-large")
    fig.subplots_adjust(hspace=0.05, wspace=0.1)

    gs = gridspec.GridSpec(nrow, ncol, width_ratios=[6, 1, 2])

    ax1 = plt.subplot(gs[0])
    if style == 'l':
        ax1.plot(t, datain, markersize=1, lw=linewidth, label=data_label, color=colorin)
        ax1.plot(t, datain2, markersize=1, lw=linewidth, label=data_label, color=colorin2)
    elif style == 'p':
        ax1.plot(t, datain, 'bo', markersize=1, lw=linewidth, label=data_label, color=colorin)
        ax1.plot(t, datain2, 'bo', markersize=1, lw=linewidth, label=data_label2, color=colorin2)
    ax1.set_xlabel(xlabel_1)
    ax1.set_ylabel(ylabel_1)
    ax1.set_title(title_1)
    ax1.grid(True, linestyle='dotted', linewidth=0.2)
    ax1.legend(loc='upper left', frameon=False)
    ax1.set_ylim(binstr, binend)
    ax1.tick_params(top=True, bottom=True, left=True, right=True,
                    labelleft=True, labelbottom=True)

    ax2 = plt.subplot(gs[1])
    ax2.plot(y, x, lw=linewidth, label=data_label, color=colorin)
    ax2.plot(y2, x2, lw=linewidth, label=data_label, color=colorin2)
    ax2.set_xlabel(xlabel_2)
    ax2.set_ylabel(ylabel_2)
    ax2.set_title(title_2)
    ax2.grid(True, linestyle='dotted', linewidth=0.2)
    # ax2.legend(loc='upper right', frameon=False)
    ax2.set_ylim(binstr, binend)
    ax2.tick_params(top=True, bottom=True, left=True, right=True,
                    labelleft=False, labelbottom=True)

    ax3 = plt.subplot(gs[2])
    # ax3.errorbar(xx, yy, 0, yyerr'bo', fmt='o', lw=linewidth, label=data_label)
    ax3.errorbar(xx, yy, yyerr, 0, fmt='o', color=colorin)
    ax3.errorbar(xx, yy2, yyerr2, 0, fmt='o', color=colorin2)
    ax3.set_xlabel('')
    ax3.set_ylabel('')
    ax3.set_title('')
    ax3.grid(True, linestyle='dotted', linewidth=0.2)
    # ax3.legend(loc='upper right', frameon=False)
    ax3.set_ylim(binstr, binend)
    ax3.set_xlim(0, 5)
    ax3.tick_params(top=True, bottom=True, left=True, right=True,
                    labelleft=False, labelbottom=False)

    # plt.show()
    plt.savefig(outfilename + ".png", dpi=300, bbox_inches='tight')


def get_df_plus(df, t, block_size=50):
    for i in range(len(df.columns)):
        _name = df.columns[i] + "_run"
        df[_name] = df[df.columns[i]].rolling(window=block_size).mean()
    df['t'] = t
    return df


def plot2_gpcr_rmsd(t, output_dir):
    df_rmsdbb_alin_all_backbone = pd.read_hdf(output_dir + '/rmsdbb_alin_all_backbone.hdf', 'df')
    df_rmsdbb_alin_prot_A_backbone = pd.read_hdf(output_dir + '/rmsdbb_alin_prot_A_backbone.hdf', 'df')
    df_rmsdbb_alin_prot_R_backbone = pd.read_hdf(output_dir + '/rmsdbb_alin_prot_R_backbone.hdf', 'df')
    df_rmsdbb_alin_OBS_backbone = pd.read_hdf(output_dir + '/rmsdbb_alin_OBS_backbone.hdf', 'df')

    df_rmsdbb_alin_all_backbone = get_df_plus(df_rmsdbb_alin_all_backbone,t )
    df_rmsdbb_alin_prot_A_backbone = get_df_plus(df_rmsdbb_alin_prot_A_backbone,t)
    df_rmsdbb_alin_prot_R_backbone = get_df_plus(df_rmsdbb_alin_prot_R_backbone,t)
    df_rmsdbb_alin_OBS_backbone = get_df_plus(df_rmsdbb_alin_OBS_backbone,t)

    list_df = [df_rmsdbb_alin_prot_R_backbone,
               df_rmsdbb_alin_prot_A_backbone,
               df_rmsdbb_alin_all_backbone,
               df_rmsdbb_alin_OBS_backbone]
    list_names = ["Rec", "Ga", "All", "OBS"]
    list_printligand = [False, False, True, True]

    nrow = 2
    ncol = 2
    figheight = 6
    figwidth = 18
    linewidth = 1
    maxylim = 10

    fig, axs = plt.subplots(nrow, ncol, figsize=(figwidth, figheight))
    fig.subplots_adjust(hspace=0.5, wspace=0.1)
    axs = axs.ravel()

    for i in range(4):
        _df = list_df[i]

        # plot raw data
        if (list_names[i] != list_names[3]):
            axs[i].plot(_df['t'], _df['all'], 'bo', markersize=0.3, lw=linewidth, label='', color='gray')
            axs[i].plot(_df['t'], _df['prot_A'], 'bo', markersize=0.3, lw=linewidth, label='', color='lightblue')
        axs[i].plot(_df['t'], _df['prot_R'], 'bo', markersize=0.3, lw=linewidth, label='', color='pink')
        if list_printligand[i]:
            axs[i].plot(_df['t'], _df['lig'], 'bo', markersize=0.1, lw=linewidth, label='', color='cyan')
            axs[i].plot(_df['t'], _df['lig_last'], 'bo', markersize=0.1, lw=linewidth, label='', color='brown')

        # plot running average
        if (list_names[i] != list_names[3]):
            axs[i].plot(_df['t'], _df['all_run'], lw=linewidth, label='All', color='black')
            axs[i].plot(_df['t'], _df['prot_A_run'], lw=linewidth, label='Ga', color='blue')
        axs[i].plot(_df['t'], _df['prot_R_run'], lw=linewidth, label='Rec', color='red')
        if list_printligand[i]:
            axs[i].plot(_df['t'], _df['lig_run'], lw=linewidth, label='Ligand', color='green')
            axs[i].plot(_df['t'], _df['lig_last_run'], lw=linewidth, label='Ligand(last)', color='orange')

        axs[i].set_xlabel('Time (ns)')
        axs[i].set_ylabel('Backbone RMSD')
        axs[i].grid(True, linestyle='dotted', linewidth=0.2)
        axs[i].legend(loc='upper right')
        axs[i].set_title('Alignment by ' + list_names[i])
        if list_printligand[i]:
            axs[i].set_ylim(0, 6)
        else:
            axs[i].set_ylim(0, maxylim)


def plot2_gpcr_rmsd2(t, output_dir):
    data_label = output_dir.replace("/output", "").split("/")[-1]

    df_rmsdbb_alin_all_backbone = pd.read_hdf(output_dir + '/rmsdbb_alin_all_backbone.hdf', 'df')
    df_rmsdbb_alin_prot_A_backbone = pd.read_hdf(output_dir + '/rmsdbb_alin_prot_A_backbone.hdf', 'df')
    df_rmsdbb_alin_prot_R_backbone = pd.read_hdf(output_dir + '/rmsdbb_alin_prot_R_backbone.hdf', 'df')
    df_rmsdbb_alin_OBS_backbone = pd.read_hdf(output_dir + '/rmsdbb_alin_OBS_backbone.hdf', 'df')

    df_rmsdbb_alin_all_backbone = get_df_plus(df_rmsdbb_alin_all_backbone,t )
    df_rmsdbb_alin_prot_A_backbone = get_df_plus(df_rmsdbb_alin_prot_A_backbone,t)
    df_rmsdbb_alin_prot_R_backbone = get_df_plus(df_rmsdbb_alin_prot_R_backbone,t)
    df_rmsdbb_alin_OBS_backbone = get_df_plus(df_rmsdbb_alin_OBS_backbone,t)

    list_df = [df_rmsdbb_alin_prot_R_backbone,
               df_rmsdbb_alin_prot_A_backbone,
               df_rmsdbb_alin_all_backbone,
               df_rmsdbb_alin_OBS_backbone]
    list_names = ["Rec", "Ga", "All", "OBS"]
    list_printligand = [False, False, True, True]


    for i in range(4):
        _df = list_df[i]
        outfilename = output_dir + "/rmsdbb_alin_"+list_names[i]+"_backbone"
        nrow = 1
        ncol = 1
        figheight = 2.5
        figwidth = 8
        linewidth = 1
        maxylim = 10

        fig, axs = plt.subplots(nrow, ncol, figsize=(figwidth, figheight))

        # plot raw data
        if (list_names[i] != list_names[3]):
            axs.plot(_df['t'], _df['all'], 'bo', markersize=0.3, lw=linewidth, label='', color='gray')
            axs.plot(_df['t'], _df['prot_A'], 'bo', markersize=0.3, lw=linewidth, label='', color='lightblue')
        axs.plot(_df['t'], _df['prot_R'], 'bo', markersize=0.3, lw=linewidth, label='', color='pink')
        if list_printligand[i]:
            axs.plot(_df['t'], _df['lig'], 'bo', markersize=0.1, lw=linewidth, label='', color='cyan')
            axs.plot(_df['t'], _df['lig_last'], 'bo', markersize=0.1, lw=linewidth, label='', color='brown')

        # plot running average
        if (list_names[i] != list_names[3]):
            axs.plot(_df['t'], _df['all_run'], lw=linewidth, label='All', color='black')
            axs.plot(_df['t'], _df['prot_A_run'], lw=linewidth, label='Ga', color='blue')
        axs.plot(_df['t'], _df['prot_R_run'], lw=linewidth, label='Rec', color='red')
        if list_printligand[i]:
            axs.plot(_df['t'], _df['lig_run'], lw=linewidth, label='Ligand', color='green')
            axs.plot(_df['t'], _df['lig_last_run'], lw=linewidth, label='Ligand(last)', color='orange')

        axs.set_xlabel('Time (ns)')
        axs.set_ylabel('Backbone RMSD')
        axs.grid(True, linestyle='dotted', linewidth=0.2)
        axs.legend(loc='upper right')
        axs.set_title('Alignment by ' + list_names[i]+' ('+data_label+')')
        if list_printligand[i]:
            axs.set_ylim(0, 6)
        else:
            axs.set_ylim(0, maxylim)

        plt.savefig(outfilename+".png", dpi=300, bbox_inches='tight')


def plot2_3_7_lock_switch(t, output_dir):
    # pickle.load
    with open (output_dir+'/gpcr_dist_3.32_7.43.p', 'rb') as fp:
        dist = pickle.load(fp)
    #interface
    datain=dist
    nbin = 10
    binstr = 0
    binend = 10
    title='3-7 Lock Switch - Distance D3.32(OD1/OD2) Y7.43(OH) (Å )'
    title_1=''
    title_2=''
    xlabel_1='Time (ns)'
    ylabel_1='Distane'
    xlabel_2='ratio'
    ylabel_2=''
    data_label=output_dir.replace("/output","").split("/")[-1]
    outfilename=output_dir+'/gpcr_dist_3.32_7.43'

    plot_evol_hist_stat(t,datain,nbin,binstr,binend,
                        title,title_1,title_2,xlabel_1,ylabel_1,xlabel_2,ylabel_2,data_label,outfilename,colorin="blue")

def plot2_3_ligandN(t, output_dir):
    # pickle.load
    with open (output_dir+'/gpcr_dist_3.32_ligandN.p', 'rb') as fp:
        dist = pickle.load(fp)
    #interface
    datain=dist
    nbin = 10
    binstr = 0
    binend = 10
    title='3-L Aminergic - Distance D3.32(OD1/OD2) Ligand(N)  (Å )'
    title_1=''
    title_2=''
    xlabel_1='Time (ns)'
    ylabel_1='Distane'
    xlabel_2='ratio'
    ylabel_2=''
    data_label=output_dir.replace("/output","").split("/")[-1]
    outfilename = output_dir + '/gpcr_dist_3.32_ligandN'

    plot_evol_hist_stat(t,datain,nbin,binstr,binend,
                        title,title_1,title_2,xlabel_1,ylabel_1,xlabel_2,ylabel_2,data_label,outfilename,colorin="blue")

def plot2_5_7_water_bridge(t, output_dir):
    # pickle.load
    with open(output_dir + '/gpcr_dist_5.58_7.53.p', 'rb') as fp:
        dist = pickle.load(fp)
    with open(output_dir + '/gpcr_dist_ca_5.58_7.53.p', 'rb') as fp:
        dist_CA = pickle.load(fp)
    # interface
    datain = dist
    datain2 = dist_CA
    nbin = 10
    binstr = 0
    binend = 20
    title = '5-7 Water Bridge - Distance Y5.58(OH) Y7.53(OH) (Å )'
    title_1 = ''
    title_2 = ''
    xlabel_1 = 'Time (ns)'
    ylabel_1 = 'Distane'
    xlabel_2 = 'ratio'
    ylabel_2 = ''
    data_label = output_dir.replace("/output", "").split("/")[-1]
    data_label2 = output_dir.replace("/output", "").split("/")[-1] + "(CA-CA dist)"
    outfilename = output_dir + '/gpcr_dist_5.58_7.53'

    plot_evol_hist_stat2(t, datain, datain2,
                         nbin, binstr, binend,
                         title, title_1, title_2, xlabel_1, ylabel_1, xlabel_2, ylabel_2, data_label, data_label2,
                         outfilename,
                         colorin="blue", colorin2="red")


def plot2_3_6_ionic_lock_switch(t, output_dir):
    # pickle.load
    with open (output_dir+'/gpcr_dist_ca_3.50_6.30.p', 'rb') as fp:
        dist = pickle.load(fp)
    #interface
    datain=dist
    # binstr = int(round(min(dist),0))
    # binend = int(round(max(dist),0)+1)
    binstr = 0
    binend = 20
    nbin = binend-binstr+1
    title='3-6 Ionic Lock Switch - Distance R3.50(CA) E6.30(CA) (Å )'
    title_1=''
    title_2=''
    xlabel_1='Time (ns)'
    ylabel_1='Distane'
    xlabel_2='ratio'
    ylabel_2=''
    data_label=output_dir.replace("/output","").split("/")[-1]
    outfilename = output_dir + '/gpcr_dist_ca_3.50_6.30'

    plot_evol_hist_stat(t,datain,nbin,binstr,binend,
                        title,title_1,title_2,xlabel_1,ylabel_1,xlabel_2,ylabel_2,data_label,outfilename,colorin="blue")

def plot2_rmsdbb_NPxxY(t, output_dir):
    # pickle.load
    with open (output_dir+'/gpcr_rmsdbb_NPxxY.p', 'rb') as fp:
        dist = pickle.load(fp)
    #interface
    datain=dist
    nbin = 10
    binstr = 0
    binend = 2
    title='NPxxY Backbone RMSD (Å )'
    title_1=''
    title_2=''
    xlabel_1='Time (ns)'
    ylabel_1='Distane'
    xlabel_2='ratio'
    ylabel_2=''
    data_label=output_dir.replace("/output","").split("/")[-1]
    outfilename = output_dir + '/gpcr_rmsdbb_NPxxY'

    plot_evol_hist_stat(t,datain,nbin,binstr,binend,
                        title,title_1,title_2,xlabel_1,ylabel_1,xlabel_2,ylabel_2,data_label,outfilename,colorin="brown")


def plot2_toggle_switch(t,output_dir,rec_df):
    list_chi=['208','383','342']
    #list_chi=["5.58", "7.53", "6.48"]
    #list_bwidx=convert_resid_to_bwidx(list_chi,rec_df)
    list_bwidx =["5.58", "7.53", "6.48"]

    for i in range(len(list_chi)):
        #_sel=list_chi[i]
        _sel = list_bwidx[i]

        # pickle.load
        with open (output_dir+'/gpcr_chi1_'+_sel+'.p', 'rb') as fp:
            datain = pickle.load(fp)
        nbin = 30
        binstr = -180
        binend = 180
        title="Chi1 "+list_bwidx[i]+"("+list_chi[i]+")"
        title_1=''
        title_2=''
        xlabel_1='Time (ns)'
        ylabel_1='angle'
        xlabel_2='ratio'
        ylabel_2=''
        #data_label=_sel
        data_label = output_dir.replace("/output", "").split("/")[-1]
        outfilename = output_dir + '/gpcr_chi1_' + list_bwidx[i]

        plot_evol_hist_stat(t,datain,nbin,binstr,binend,
                        title,title_1,title_2,xlabel_1,ylabel_1,xlabel_2,ylabel_2,data_label,outfilename,colorin="red")


def plot2_chi2(t,output_dir,rec_df):
    list_chi=['342']
    #list_chi=["6.48"]
    #list_bwidx=convert_resid_to_bwidx(list_chi,rec_df)
    list_bwidx=["6.48"]

    for i in range(len(list_chi)):
        #_sel=list_chi[i]
        _sel = list_bwidx[i]

        # pickle.load
        with open (output_dir+'/gpcr_chi2_'+_sel+'.p', 'rb') as fp:
            datain = pickle.load(fp)
        nbin = 30
        binstr = -180
        binend = 180
        title="Chi2 "+list_bwidx[i]+"("+list_chi[i]+")"
        title_1=''
        title_2=''
        xlabel_1='Time (ns)'
        ylabel_1='angle'
        xlabel_2='ratio'
        ylabel_2=''
        #data_label=_sel
        data_label = output_dir.replace("/output", "").split("/")[-1]
        outfilename = output_dir + '/gpcr_chi2_' + list_bwidx[i]

        plot_evol_hist_stat(t,datain,nbin,binstr,binend,
                        title,title_1,title_2,xlabel_1,ylabel_1,xlabel_2,ylabel_2,data_label,outfilename,colorin="pink")

def plot2_euler(t,output_dir):
    df_euler = pd.read_hdf(output_dir+'/gpcr_euler.hdf','df') #df_euler
    for i in range(len(df_euler.columns)):
        _sel=df_euler.columns[i]
        datain=df_euler[_sel]
        nbin = 30
        binstr = -180
        binend = 180
        title=_sel
        title_1=''
        title_2=''
        xlabel_1='frame'
        ylabel_1='angle'
        xlabel_2='ratio'
        ylabel_2=''
        data_label=_sel
        plot_evol_hist_stat(t, datain, nbin, binstr, binend, title, title_1, title_2, xlabel_1, ylabel_1, xlabel_2,
                            ylabel_2, data_label, outfilename)


def plot2_interface(t, output_dir):
    # pickle.load
    with open (output_dir+'/gpcr_interface.p', 'rb') as fp:
        interface = pickle.lget_matrix_diffoad(fp)
    #interface
    datain=interface
    binstr = round(min(interface)/100)*100
    binend = round(max(interface)/100+1)*100
    nbin = 10

    title='Interface Area (Å^2)'
    title_1=''
    title_2=''
    xlabel_1='Time (ns)'
    ylabel_1='angle'
    xlabel_2='ratio'
    ylabel_2=''
    data_label="alpha"
    outfilename = output_dir+'/gpcr_interface'

    plot_evol_hist_stat(t,datain,nbin,binstr,binend,title,title_1,title_2,xlabel_1,ylabel_1,xlabel_2,ylabel_2,data_label,outfilename)


def check_convergence(keyobs, arr_output_dir, xlabel, ylabel, nbin):
    # xlabel="distance"
    # ylabel="probability"
    # title="3_7_lock_switch"
    if keyobs == "3_7_lock_switch":
        obsin = "gpcr_dist_3.32_7.43"
    if keyobs == "ionic_lock_switch":
        obsin = "gpcr_dist_ca_3.50_6.30"
    # if keyobs == "rmsdbb_NPxxY":
    #     obsin = "gpcr_rmsdbb_NPxxY"

    title = keyobs + " (nbin=" + str(nbin)+")"

    obs = []
    df = pd.DataFrame()
    min_obs = 9999
    max_obs = 9999
    for i in range(len(arr_output_dir)):
        output_dir_root = arr_output_dir[i]
        output_dir = output_dir_root + "/output"
        # print(output_dir)
        with open(output_dir + '/' + obsin + '.p', 'rb') as fp:
            obs.append(pickle.load(fp))
        min_obs = min(min_obs, min(obs[i]))
        max_obs = min(max_obs, max(obs[i]))
        df[i] = obs[i]

    nrow, ncol = 1, 2
    figwidth, figheight = 4, 2
    linewidth = 2
    fig, ax = plt.subplots(nrow, ncol, figsize=(figwidth, figheight))

    binstr = min_obs
    binend = max_obs
    # nbin=int(binend-binstr+1)
    # nbin=20
    # print(nbin)

    for i in range(len(df.columns)):
        nf = len(df[i])
        nf_str = int(round(nf / 2, 0))
        # print("xx",i,nf_str,nf)
        datain = list(df[i][nf_str:nf])
        xbin, hist_n = bining_last_half(datain, nbin, binstr, binend)
        # smooth plot
        x, y = get_smooth(xbin, hist_n)
        ax.plot(x, y, lw=linewidth, label=arr_output_dir[i].split("/")[-1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(binstr, binend)
    ax.grid(True, linestyle='dotted', linewidth=0.2)
    ax.legend(loc='upper right', frameon=False)


def check_convergence_wrapper(arr_output_dir):
    arr_output_dir = sorted(arr_output_dir)
    xlabel="Distance"
    ylabel="Probability"
    keyobs="3_7_lock_switch"
    check_convergence(keyobs,arr_output_dir,xlabel,ylabel,nbin=20)
    keyobs="ionic_lock_switch"
    check_convergence(keyobs,arr_output_dir,xlabel,ylabel,nbin=20)
    # keyobs="rmsdbb_NPxxY"
    # xlabel="RMSD"
    # check_convergence(keyobs,arr_output_dir,xlabel,ylabel,nbin=20)


def check_convergence2(output_dir, label):
    xli = [
        "rmsdbb_alin_Ga_backbone.png",
        "rmsdbb_alin_Rec_backbone.png",
        "rmsdbb_alin_All_backbone.png",
        "rmsdbb_alin_OBS_backbone.png",
        "gpcr_dist_3.32_ligandN.png",
        "gpcr_dist_3.32_7.43.png",
        "gpcr_dist_5.58_7.53.png",
        "gpcr_dist_ca_3.50_6.30.png",
        "gpcr_chi1_6.48.png",
        "gpcr_chi2_6.48.png",
        "gpcr_chi1_5.58.png",
        "gpcr_chi1_7.53.png"
    ]
    # "gpcr_rmsdbb_NPxxY.png",

    images = [Image.open(output_dir+"/"+x) for x in xli]
    widths, heights = zip(*(i.size for i in images))

    # total_width = sum(widths)
    # max_height = max(heights)
    # new_im = Image.new('RGB', (total_width, max_height))
    #
    # x_offset = 0
    # for im in images:
    #   new_im.paste(im, (x_offset,0))
    #   x_offset += im.size[0]

    maxw = max(widths)
    maxh = max(heights)
    total_width = maxw * 3
    max_height = maxh * 4
    new_im = Image.new('RGB', (total_width, max_height), color=(255, 255, 255, 0))
    new_im.paste(images[0], (0, 0))
    new_im.paste(images[1], (maxw, 0))
    new_im.paste(images[2], (maxw * 2, 0))
    new_im.paste(images[3], (0, maxh))
    new_im.paste(images[4], (maxw, maxh))
    new_im.paste(images[5], (maxw * 2, maxh))
    new_im.paste(images[6], (0, maxh * 2))
    new_im.paste(images[7], (maxw, maxh * 2))
    new_im.paste(images[8], (maxw * 2, maxh * 2))
    new_im.paste(images[9], (0, maxh * 3))
    new_im.paste(images[10], (maxw, maxh * 3))
    new_im.paste(images[11], (maxw * 2, maxh * 3))

    new_im.save('figures/ck_'+label+'.jpg')
    return new_im