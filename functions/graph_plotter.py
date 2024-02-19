"""
graph plotting function.

Author:
    Ara Abramyan

"""
import glob

import matplotlib as mpl
import numpy as np

mpl.use('Agg')
import matplotlib.pyplot as plt


def plot_function(mymatrix, my_xnames, my_ynames, pltfile, colors, matrix):
    plt.figure(figsize=(15, 15))
    sz = 25
    
    plt.imshow(np.matrix(mymatrix), interpolation='nearest', cmap=colors)
    
    lim = max(np.hstack(mymatrix).min(), np.hstack(mymatrix).max(), key=abs)
    # ticks_position = np.arange(-lim, lim+1, 5.0)
    # print ticks_position
    print(matrix, 'limit =', lim)
    
    cbar = plt.colorbar(orientation='horizontal', aspect=2.0, shrink=0.2, pad=0.1, fraction=.12, ticks=(-int(lim), 0, int(lim)))
    cbar.ax.tick_params(labelsize=sz)
    # cbar.locator()
    
    plt.clim(-lim, lim)
    ax = plt.gca()
    plt.xticks(np.arange(0, len(my_xnames), 1), size=sz, rotation=90)
    ax.set_xticklabels(my_xnames)
    plt.yticks(np.arange(0, len(my_ynames), 1), size=sz)
    ax.set_yticklabels(my_ynames)
    
    ax.xaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
    ax.yaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    plt.tight_layout()
    plt.savefig(pltfile, format='pdf', transparent=True)


def plot_pia_matrix(metrics, parms, colors):
    arr_diff_matrix_file_names = glob.glob(parms.outdir_av_and_diff_mat + '*_vs_*.csv')
    
    for matrix in metrics:
        for fn in arr_diff_matrix_file_names:
            pltfile = fn[:-4].replace(parms.outdir_av_and_diff_mat, parms.plot_dir) + '.pdf'
            if matrix in fn.split('/')[-1]:
                mymatrix = np.loadtxt(fn, delimiter=",", dtype=float)
                if matrix == parms.matrix1 or matrix == parms.matrix2:
                    plot_function(mymatrix, parms.helix_names[parms.protein_name], parms.helix_names[parms.protein_name], pltfile, colors, matrix)
                if matrix == parms.matrix3:
                    plot_function(mymatrix, parms.binding_BW[parms.protein_name], parms.binding_BW[parms.protein_name], pltfile, colors, matrix)
                if matrix == parms.matrix4:
                    plot_function(mymatrix, parms.binding_BW[parms.protein_name], '', pltfile, colors, matrix)
