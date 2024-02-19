import matplotlib
import matplotlib.pyplot as plt
import numpy as np 
from matplotlib.ticker import FuncFormatter
from matplotlib import gridspec
import seaborn as sns
from scipy.interpolate import make_interp_spline, BSpline

def bining(datain, nbin, binstr, binend):
    # do bining for all data
    hist, bins = np.histogram(datain, bins=nbin, range=[binstr, binend], density=False)
    hist_n = hist / sum(hist)

    # manipulate bin
    bin_shift = (bins[1] - bins[0]) / 2
    nbin_center = len(bins) - 1
    xbin = bins + bin_shift
    xbin = xbin[0:nbin_center]

    return xbin, hist_n

def get_smooth(x,y):
    # 300 represents number of points to make between T.min and T.max
    xnew = np.linspace(x.min(), x.max(), 300)
    # type: BSpline
    spl = make_interp_spline(x, y, k=1)
    power_smooth = spl(xnew)
    return(xnew, power_smooth)



def get_palette(rgb_str, rgb_end, n, show=False):
    ''' given RGB color and get color palette '''
    rgb_str = [(rgb_str[i])/255 for i in range(3)]
    rgb_end = [(rgb_end[i])/255 for i in range(3)]
    palette = []
    rgb_step = [round((rgb_end[x]-rgb_str[x])/(n-1),2)for x in range(3)]
    if show: fig, ax = plt.subplots(figsize=(10, 2))
    for i in range(n):
        rgb = [round(rgb_str[x]+i*rgb_step[x],2) for x in range(3)]
        palette.append(rgb)
        if show: ax.scatter(i,0.5, c=[rgb],s=800, marker='s')
    if show:
        rgb = rgb_str
        ax.scatter(0,0.4, c=[rgb],s=600)
        rgb = rgb_end
        ax.scatter(i,0.6, c=[rgb],s=600)
        ax.tick_params(top=True, bottom=True, left=True, right=True,labelleft=False, labelbottom=False)
    return palette



def distributionplot(indata, inlabels, incolors, in_xlabel, in_ylabel, in_subtitle, \
                     in_xlim, in_ylim, figheight, figwidth, lagend=True, for_ai=False):

    # define fig & axs
    fig, ax = plt.subplots(figsize=(figwidth, figheight))
    fig.subplots_adjust(hspace=0.5, wspace=0.4)
    for i in range(len(indata)):
        sns.distplot(indata[i], hist=False, color = incolors[i], kde_kws={"shade": False}, label = inlabels[i])
    ax.set_xlabel(in_xlabel)
    ax.set_ylabel(in_ylabel)
    ax.set_title(in_subtitle)  # subtitle
    ax.grid(True, linestyle='dotted', linewidth=0.2)
    ax.legend(loc='upper right', frameon=False)
    ax.set_xlim(in_xlim)
    ax.set_ylim(in_ylim)
    ax.tick_params(top=True, bottom=True, left=True, right=True,
                   labelleft=True, labelbottom=True)

    if for_ai:
        # another plot without labels
        inlabels = []
        for i in range(len(indata)):
            inlabels.append("          ")
        in_xlabel = ""
        in_ylabel = ""
        in_subtitle = ""

        fig, ax = plt.subplots(figsize=(figwidth, figheight))
        fig.subplots_adjust(hspace=0.5, wspace=0.4)
        for i in range(len(indata)):
            if lagend:
                sns.distplot(indata[i], hist=False, color = incolors[i], kde_kws={"shade": False}, label = inlabels[i])
            else:
                sns.distplot(indata[i], hist=False, color = incolors[i], kde_kws={"shade": False})

        # for ix in range(len(refdata)):
        #     ax.axvline(x=refdata[ix],linewidth=1, color=refcolors[ix],label = reflabels[ix])

        ax.set_xlabel(in_xlabel)
        ax.set_ylabel(in_ylabel)
        ax.set_title(in_subtitle)  # subtitle
        ax.grid(True, linestyle='dotted', linewidth=0.2)
        ax.legend(loc='upper right', frameon=False)
        ax.set_xlim(in_xlim)
        ax.set_ylim(in_ylim)
        ax.tick_params(top=True, bottom=True, left=True, right=True,
                       labelleft=False, labelbottom=False)


def distributionplot2(indata,inlabels,incolors,refdata,reflabels,refcolors,in_xlabel,in_ylabel,in_subtitle, \
             in_xlim,in_ylim,figheight,figwidth,save2filename,lagend=True):

    # define fig & axs
    fig, ax = plt.subplots(figsize=(figwidth, figheight))
    fig.subplots_adjust(hspace=0.5, wspace=0.4)
    for i in range(len(indata)):
        sns.distplot(indata[i], hist=False, color = incolors[i], kde_kws={"shade": False}, label = inlabels[i])
    for ix in range(len(refdata)):
        ax.axvline(x=refdata[ix],linewidth=1, color=refcolors[ix],label = reflabels[ix]+"("+str(refdata[ix])+")")
    ax.set_xlabel(in_xlabel)
    ax.set_ylabel(in_ylabel)
    ax.set_title(in_subtitle)  # subtitle
    ax.grid(True, linestyle='dotted', linewidth=0.2)
    ax.legend(loc='upper right', frameon=False)
    ax.set_xlim(in_xlim)
    ax.set_ylim(in_ylim)
    ax.tick_params(top=True, bottom=True, left=True, right=True,
                   labelleft=True, labelbottom=True)

    # another plot without labels
    inlabels = ["          ", "   "]
    in_xlabel = ""
    in_ylabel = ""
    in_subtitle = ""

    fig, ax = plt.subplots(figsize=(figwidth, figheight))
    fig.subplots_adjust(hspace=0.5, wspace=0.4)
    for i in range(len(indata)):
        if lagend:
            sns.distplot(indata[i], hist=False, color = incolors[i], kde_kws={"shade": False}, label = inlabels[i])
        else:
            sns.distplot(indata[i], hist=False, color = incolors[i], kde_kws={"shade": False})

    # for ix in range(len(refdata)):
    #     ax.axvline(x=refdata[ix],linewidth=1, color=refcolors[ix],label = reflabels[ix])

    ax.set_xlabel(in_xlabel)
    ax.set_ylabel(in_ylabel)
    ax.set_title(in_subtitle)  # subtitle
    ax.grid(True, linestyle='dotted', linewidth=0.2)
    ax.legend(loc='upper right', frameon=False)
    ax.set_xlim(in_xlim)
    ax.set_ylim(in_ylim)
    ax.tick_params(top=True, bottom=True, left=True, right=True,
                   labelleft=False, labelbottom=False)
    plt.savefig(save2filename, dpi=300,bbox_inches='tight')


def heatmap1D_diff(pair_a_bw,pair_a_contact,pair_a_name,pair_b_bw,pair_b_contact,pair_b_name,
                   tick_max=1, tick_min=-1, tick_gap=0.5):

    ls_diff = []
    ls_diff_labels = []
    for i in range(len(pair_b_bw)):
        if (pair_a_bw[i] == pair_b_bw[i]):
            ls_diff_labels.append(pair_b_bw[i])
            ls_diff.append(round(pair_a_contact[i]-pair_b_contact[i],2))
        else:
            ls_diff_labels.append(pair_a_bw[i]+"-"+pair_b_bw[i])
            ls_diff.append(round(pair_a_contact[i]-pair_b_contact[i],2))

    # datain
    ls_data = ls_diff
    ls_labels = ls_diff_labels
    title_diff = pair_a_name+"-"+pair_b_name

    # tick
    #tick_max = 1
    #tick_min = -1
    #tick_gap = 0.5
    tick_label_size = 12

    # convert datain from list to np array
    np_data =np.array([ls_data])
    x = np.empty([5,np_data.shape[1]])
    x[:,:] = np_data

    plt.figure(figsize=(20,20))
    #plt.imshow(x, interpolation='nearest', cmap=plt.cm.RdBu, aspect=0.25)
    plt.imshow(x, interpolation='nearest', cmap=plt.cm.bwr, aspect=0.25)
    cbar = plt.colorbar(orientation='vertical', shrink=0.15, pad=0.02, aspect=10,
                        ticks = (np.arange(tick_min,tick_max+tick_gap,tick_gap)))
    cbar.ax.tick_params('both', length=10, width=2, which='major', labelsize=tick_label_size)
    # l=max(np.hstack(x).min(), np.hstack(x).max(), key=abs)
    plt.clim(tick_min,tick_max)
    ax=plt.gca()

    ax.tick_params('both', length=2, width=0, which='major')
    ax.set_xticks(np.arange(0,len(ls_labels),1))
    ax.set_xticklabels(ls_labels, fontsize=14, rotation=90)
    ax.get_yaxis().set_visible(False)
    ax.set_ylim(0,4)

    ax.set_xlabel('', fontsize = 10)
    ax.set_title(title_diff, fontsize = 25)

    grid = np.reshape(ls_data, (-1, len(ls_labels)))
    for (j,i),label in np.ndenumerate(grid):
        ax.text(i,j+2.5,label,ha='center',va='top', fontsize=13)
        #ax.text(i,j,label,ha='center',va='center')
    for black_border_line in np.arange(0.5, 0.5+len(ls_labels), 1):
        plt.axvline(x=black_border_line, color='black')





def heatmap1D(labels, rawdata, title='', tick_max=1, tick_min=0, tick_gap=0.5):
    # datain
    ls_data = rawdata
    ls_labels = labels
    title_diff = title

    tick_label_size = 12

    # convert datain from list to np array
    np_data = np.array([ls_data])
    x = np.empty([5, np_data.shape[1]])
    x[:, :] = np_data

    plt.figure(figsize=(20, 20))
    #plt.imshow(x, interpolation='nearest', cmap=plt.cm.YlOrRd, aspect=0.25)
    plt.imshow(x, interpolation='nearest', cmap=plt.cm.Reds, aspect=0.25)
    cbar = plt.colorbar(orientation='vertical', shrink=0.15, pad=0.02, aspect=10,
                        ticks=(np.arange(tick_min, tick_max + tick_gap, tick_gap)))
    cbar.ax.tick_params('both', length=10, width=2, which='major', labelsize=tick_label_size)
    # l=max(np.hstack(x).min(), np.hstack(x).max(), key=abs)
    plt.clim(tick_min, tick_max)
    ax = plt.gca()

    ax.tick_params('both', length=2, width=0, which='major')
    ax.set_xticks(np.arange(0, len(ls_labels), 1))
    # ax.set_xticklabels(ls_labels, fontsize=14, rotation=90)
    ax.set_xticklabels(ls_labels, fontsize=14, rotation=45)
    ax.get_yaxis().set_visible(False)
    ax.set_ylim(0, 4)

    ax.set_xlabel('', fontsize=10)
    ax.set_title(title_diff, fontsize=25)

    grid = np.reshape(ls_data, (-1, len(ls_labels)))
    for (j, i), label in np.ndenumerate(grid):
        ax.text(i, j + 2.5, label, ha='center', va='top', fontsize=13)
        # ax.text(i,j,label,ha='center',va='center')
    for black_border_line in np.arange(0.5, 0.5 + len(ls_labels), 1):
        plt.axvline(x=black_border_line, color='black')
    

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", cbar_min=0, cbar_max=100, **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # set min/max
    for im in plt.gca().get_images():
        im.set_clim(cbar_min, cbar_max)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.1, pad=0.1, aspect=3, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
   
    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-45, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="black", linestyle='-', linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def heatmap_for_manuscript(data, row_labels, col_labels, ax=None,cbar_kw={}, cbarlabel="", **kwargs):
     if not ax:
         ax = plt.gca()


     # Plot the heatmap
     im = ax.imshow(data, **kwargs)

     # Create colorbar
     cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
     cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

     # We want to show all ticks...
     ax.set_xticks(np.arange(data.shape[1]))
     ax.set_yticks(np.arange(data.shape[0]))
     # ... and label them with the respective list entries.
     #ax.set_xticklabels(['' for _ in np.arange(data.shape[1])])
     ax.set_yticklabels(['' for _ in np.arange(data.shape[0])])
     ax.set_xticklabels(col_labels)


     # Let the horizontal axes labeling appear on top.
     ax.tick_params(top=True, bottom=False,
                    labeltop=True, labelbottom=False)

     # Rotate the tick labels and set their alignment.
     plt.setp(ax.get_xticklabels(), rotation=-45, ha="right",
              rotation_mode="anchor")

     # Turn spines off and create white grid.
     for edge, spine in ax.spines.items():
         spine.set_visible(False)

     ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
     ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
     ax.grid(which="minor", color="white", linestyle='-', linewidth=0.8)
     ax.tick_params(which="minor", bottom=False, left=False)

     return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
        print('threshold',threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            #kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def _truncated_colormap(cmap_name, minval=0.0, maxval=1.0, n=-1):
    cmap = plt.get_cmap(cmap_name)
    if n == -1:
        n = cmap.N
    
    _crange = np.linspace(minval, 1, n)
    cmap_range = []
    for c in _crange:
        if c < maxval:
            cmap_range.append(c)
    cmap_range.append(1.0)
    cmap_range = np.array(cmap_range)
        
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
         'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
         cmap(cmap_range))
    return new_cmap





#
from natsort import natsorted
from statistics import mean, stdev
import matplotlib.pyplot as plt
from matplotlib import gridspec

import sys, os
from pathlib import Path

import logging
logging.getLogger().setLevel(logging.CRITICAL)
import warnings
warnings.filterwarnings('ignore')


#
# home = str(Path.home())
# toolset_dir = home + '/repositories/pia-x/toolset'
# toolset_common = toolset_dir + "/common"
# toolset_gpcr = toolset_dir + "/gpcr"
# sys.path.insert(0, toolset_dir)
# sys.path.insert(0, toolset_common)
# sys.path.insert(0, toolset_gpcr)
# from common_gpcr import *
# from parameters_drd3 import *
# from plot_master import *
# from anal_rawdata import *
#
# # def rgb_to_hex(rgb):
# #     return '%02x%02x%02x' % rgb
#
#
# # pre-define color map
# color_dict = {"d2gi_bro":    "cyan",
#               "d2gi_qnp":    "blue",
#               "d2go_bro":    "yellow",
#               "d2go_qnp":    "gold",
#               "d3gi_qnp":    "lime",
#               "d3gi_prm":    "limegreen",
#               "drd3_prm7":   "limegreen",
#               "d3gi_pd":     "mediumspringgreen",
#               "drd3_pd2":    "mediumspringgreen",
#               "d3go_qnp":    "magenta",
#               "d3go_prm":    "purple",
#               "drd3ao_prm3": "purple",
#               "d3go_pd":     "pink",
#               "drd3ao_pd2":  "pink",
#               "d2gi":        "cyan",
#               "d2go":        "yellow",
#               "d3gi":        "lime",
#               "d3go":        "magenta",
#               }
#
#
# def get_name_mean_stdev(indf, obs, by_receptor=False):
#     '''
#     given dataframe and given observable
#     return the name average and color
#     '''
#     x_name = []
#     x_mean = []
#     x_stdev = []
#     ls_color = []
#
#     if by_receptor:
#         for job in natsorted(set(indf["Receptor"]), key=lambda y: y.lower()):
#             ls_color.append(color_dict[job.split(".")[0]])
#             x_name.append(job)
#             indf2 = indf[indf['Receptor'] == job]
#             x_mean.append(statistics.mean(indf2[obs]))
#             x_stdev.append(statistics.stdev(indf2[obs]))
#
#     else:
#         for job in natsorted(set(indf["Jobname"]), key=lambda y: y.lower()):
#             ls_color.append(color_dict[job.split(".")[0]])
#             x_name.append(job)
#             indf2 = indf[indf['Jobname'] == job]
#             x_mean.append(statistics.mean(indf2[obs]))
#             x_stdev.append(statistics.stdev(indf2[obs]))
#
#     return x_name, x_mean, x_stdev, ls_color
#
#
# def barplot_df_obs (nrow, ncol, df_prod, df_bs, obs, in_ylabel="y_label"):
#     # nrow = 3
#     # ncol = 1
#     # figheight = 12
#     # figwidth = 8
#
#     # fig, axs = plt.subplots(nrow, ncol, figsize=(figwidth, figheight))
#     # fig.subplots_adjust(hspace=0.5, wspace=0.1)
#
#     gs = gridspec.GridSpec(nrow, ncol)
#
#     # rawdata dataframe
#     indf = df_prod
#     x_name, x_mean, x_stdev, x_color = get_name_mean_stdev(indf, obs)
#
#     ax = plt.subplot(gs[0])
#     ax.bar(x_name, x_mean, yerr=x_stdev, color=x_color)
#     ax.set_xlabel("")
#     ax.set_ylabel(in_ylabel)
#     ax.set_title(obs)
#     ax.grid(True, linestyle='dotted', linewidth=0.2)
#     ax.legend(loc='upper right', frameon=False)
#     #ax.set_ylim(binstr, binend)
#     ax.tick_params(top=True, bottom=True, left=True, right=True,
#                     labelleft=True, labelbottom=True)
#     ax.tick_params(axis='x', labelrotation=90)
#
#     # bootstrapping dataframe
#     indf = bootstrap_stats(df_bs, obs)
#     x_name, x_mean, x_stdev, x_color = get_name_mean_stdev(indf, obs)
#
#     ax = plt.subplot(gs[1])
#     ax.bar(x_name, x_mean, yerr=x_stdev, color=x_color)
#     ax.set_xlabel("")
#     ax.set_ylabel(in_ylabel)
#     #ax.set_title(obs)
#     ax.grid(True, linestyle='dotted', linewidth=0.2)
#     ax.legend(loc='upper right', frameon=False)
#     #ax.set_ylim(binstr, binend)
#     ax.tick_params(top=True, bottom=True, left=True, right=True,
#                     labelleft=True, labelbottom=True)
#     ax.tick_params(axis='x', labelrotation=90)
#
#
#     # bootstrapping dataframe
#     indf = bootstrap_stats(df_bs, obs, by_receptor=True)
#     x_name, x_mean, x_stdev, x_color = get_name_mean_stdev(df_bs, obs, by_receptor=True)
#
#     ax = plt.subplot(gs[2])
#     ax.bar(x_name, x_mean, yerr=x_stdev, color=x_color)
#     ax.set_xlabel("")
#     ax.set_ylabel(in_ylabel)
#     #ax.set_title(obs)
#     ax.grid(True, linestyle='dotted', linewidth=0.2)
#     ax.legend(loc='upper right', frameon=False)
#     #ax.set_ylim(binstr, binend)
#     ax.tick_params(top=True, bottom=True, left=True, right=True,
#                     labelleft=True, labelbottom=True)
#     ax.tick_params(axis='x', labelrotation=90)
#
#     return x_name, x_mean, x_stdev, x_color
#
#
#
# def plot_SER(df_prod, df_bs, obs='gpcr_rec_chi1_5.42'):
#
#     # figure parameters
#     nrow = 2
#     ncol = 2
#     figheight = 3 * nrow
#     figwidth = 4 * ncol
#
#     # define fig & axs
#     fig, axs = plt.subplots(nrow, ncol, figsize=(figwidth, figheight))
#     fig.subplots_adjust(hspace=0.5, wspace=0.4)
#     gs = gridspec.GridSpec(nrow, ncol)
#
#     indf = df_prod
#     Receptors = sorted(list(set(list(indf['Receptor']))))
#     for i in range(2):
#         datain = list(indf[indf['Receptor'] == Receptors[i]][obs])
#         xbin, hist_n = bining(datain, 36, -180, 180)
#         ax = plt.subplot(gs[i])
#         ax.plot(xbin, hist_n, 'x', markersize=1, lw=2, label="", color="blue")
#         ax.plot(xbin, hist_n, lw=2, label="", color="blue")
#         ax.set_xlabel(obs)
#         ax.set_ylabel("probability")
#         ax.set_title(Receptors[i] + " (prod)")  # subtitle
#         ax.grid(True, linestyle='dotted', linewidth=0.2)
#         ax.legend(loc='upper right', frameon=False)
#         ax.set_xlim(-180, 180)
#         ax.set_ylim(0, 0.6)
#         ax.tick_params(top=True, bottom=True, left=True, right=True,
#                        labelleft=True, labelbottom=True)
#
#     indf = df_bs
#     Receptors = sorted(list(set(list(indf['Receptor']))))
#     for i in range(2):
#         datain = list(indf[indf['Receptor'] == Receptors[i]][obs])
#         xbin, hist_n = bining(datain, 36, -180, 180)
#         ax = plt.subplot(gs[2 + i])
#         ax.plot(xbin, hist_n, 'x', markersize=1, lw=2, label="", color="red")
#         ax.plot(xbin, hist_n, lw=2, label="", color="red")
#         ax.set_xlabel(obs)
#         ax.set_ylabel("probability")
#         ax.set_title(Receptors[i] + " (bs)")  # subtitle
#         # ax.grid(True, linestyle='dotted', linewidth=0.2)
#         # ax.legend(loc='upper right', frameon=False)
#         # ax.set_xlim(-180, 180)
#         # ax.set_ylim(0, 0.6)
#         # ax.tick_params(top=True, bottom=True, left=True, right=True,
#         #                labelleft=True, labelbottom=True)
#         #
#
#
#

