import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import allel
from cvtk.utils import integerize


def correction_diagnostic_plot(diag, figsize=None, color=True):
    corr_df, models, xpreds, ypreds = diag
    if figsize is not None:
        fig, ax = plt.subplots(ncols=2, nrows=2, figsize=figsize)
    else:
        fig, ax = plt.subplots(ncols=2, nrows=2)
    labelx, labely = 0.05, 0.85
    before = corr_df[corr_df['correction'] == False]
    after = corr_df[corr_df['correction'] == True]
    if color:
        ax[0, 0].scatter(before['depth'], before['diag'], c=integerize(before['seqid']), s=5)
        ax[0, 1].scatter(before['depth'], after['diag'], c=integerize(before['seqid']), s=5)
    else:
        ax[0, 0].scatter(before['depth'], before['diag'], s=5)
        ax[0, 1].scatter(before['depth'], after['diag'], s=5)
    ax[0, 0].plot(xpreds[False], ypreds[False][0], 'r-')
    ax[0, 1].plot(xpreds[True], ypreds[True][0], 'r-')
    ax[0, 0].annotate('before correction', xy=(labelx, labely), xycoords='axes fraction')
    ax[0, 0].set_ylabel('variance')
    ax[0, 0].set_xlabel('average depth per window')
    ax[0, 1].set_xlabel('average depth per window')
    ax[0, 1].annotate('after correction', xy=(labelx, labely), xycoords='axes fraction')

    ax[1, 0].plot(xpreds[False], ypreds[False][1], 'r-')
    ax[1, 0].annotate('before correction', xy=(labelx, labely), xycoords='axes fraction')
#    ax[1, 0].axhline(y=0, color='99', zorder=1, linestyle='--')
    ax[1, 0].set_ylabel('covariance')
    if color:
        ax[1, 0].scatter(before['depth'], before['offdiag'], c=integerize(before['seqid']),
                         zorder=2, s=5)
        ax[1, 1].scatter(before['depth'], after['offdiag'], c=integerize(before['seqid']),
                         zorder=2, s=5)
    else:
        ax[1, 0].scatter(before['depth'], before['offdiag'], zorder=2, s=5)
        ax[1, 1].scatter(before['depth'], after['offdiag'], zorder=2, s=5)
    ax[1, 1].plot(xpreds[True], ypreds[True][1], 'r-')
    ax[1, 1].annotate('after correction', xy=(labelx, labely), xycoords='axes fraction')
 #   ax[1, 1].axhline(y=0, color='99', zorder=1, linestyle='--')
    ax[1, 0].set_xlabel('average depth per window')
    ax[1, 1].set_xlabel('average depth per window')
    plt.tight_layout()
    return fig, ax


def het_plot(tempfreqs, rep, time):
    p = tempfreqs.freqs[rep, time, :]
    het = 2*p*(1-p)
    midpoints = tempfreqs.gintervals.cummulative_midpoint
    seqids = tempfreqs.gintervals.seqid
    plt.scatter(midpoints, het, c=integerize(seqids))

def row_block_mean(mat, width, axis):
    """
    Diagnostic plot -- coverage density by chunk of consecutive SNPs.
    """
    nrow, ncol = mat.shape
    nblocks = np.ceil(nrow / width)
    block_means = []
    for bl in np.arange(nblocks):
        start, end = int(bl*width), int(min((bl+1)*width, nrow))
        block_means.append(mat[start:end, :].mean(axis=axis))
    return np.stack(block_means)


def rep_plot_pca(df, x=1, y=2, s=300, figsize=None, dpi=None, label=True, cmap=None, fontsize=12):
    l1, l2 = f"pc{int(x)}", f"pc{int(y)}"
    pc1, pc2 = df[l1], df[l2]
    gen, rep = df['gen'], df['rep']
    plt.scatter(pc1, pc2, c=integerize(rep), s=s, zorder=2, cmap=cmap)
    # plot lines between consecutive generations
    for rep in rep.unique():
        this_rep = df[df['rep'] == rep]
        plt.plot(this_rep[l1], this_rep[l2], '--', color='0.8', zorder=1)
    # each marker gets a label, of the current replicate
    label_df = df[[l1, l2, 'rep', 'gen']]
    if label:
        for pc1, pc2, rep, gen in label_df.itertuples(index=False):
            plt.text(pc1, pc2, s=f"{rep}", horizontalalignment='center', verticalalignment='center', fontsize=fontsize)
    plt.xlabel(l1.upper())
    plt.ylabel(l2.upper())
    plt.tight_layout()
    return plt


def rep_plot_pca2(df, x=1, y=2, s=300, figsize=None, dpi=None, label=True, cmap=None,
                  fontsize=12):
    l1, l2 = f"pc{int(x)}", f"pc{int(y)}"
    pc1, pc2 = df[l1], df[l2]
    gen, rep = df['gen'], df['rep']
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.scatter(pc1, pc2, c=integerize(rep), s=s, zorder=2, cmap=cmap)
    # plot lines between consecutive generations
    for r in rep.unique():
        this_rep = df[df['rep'] == r]
        ax.plot(this_rep[l1], this_rep[l2], '--', color='0.8', zorder=-1)
    # each marker gets a label, of the current replicate
    label_df = df[[l1, l2, 'rep', 'gen']]
    if label:
        for pc1, pc2, rep, gen in label_df.itertuples(index=False):
            ax.text(pc1, pc2, s=f"{gen}", horizontalalignment='center',
                    verticalalignment='center', fontsize=fontsize,
                    zorder=rep)
    ax.set_xlabel(l1.upper())
    ax.set_ylabel(l2.upper())
    #plt.tight_layout()
    return fig, ax

def rep_plot_pca3(df, x=1, y=2, s=300, figsize=None, dpi=None, label=True, cmap=None,
                  fontsize=12, variance=None):
    l1, l2 = f"pc{int(x)}", f"pc{int(y)}"
    pc1, pc2 = df[l1], df[l2]
    gen, rep = df['gen'], df['rep']
    ids = df['id']
    selection = df['selection']
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    groups = df.groupby("selection")
    for name, group in groups:
        ax.scatter(group['pc1'], group['pc2'], label=name)
    #ax.scatter(pc1, pc2, c=selection)
    # plot lines between consecutive generations
    for r in rep.unique():
        this_rep = df[df['rep'] == r]
        ax.plot(this_rep[l1], this_rep[l2], '--', color='0.8', zorder=-1)
    # each marker gets a label, of the current replicate
    label_df = df[[l1, l2, 'rep', 'gen', 'id', 'selection']]
    if label:
        for i, point in label_df.iterrows():
                 ax.text(point['pc1'], point['pc2'], str(point['id']),
                                 horizontalalignment='center',verticalalignment='center', fontsize=fontsize)
    if variance is not None:
        pc1=str('%.2f' % variance[0])
        pc2=str('%.2f' % variance[1])
        ax.set_xlabel(l1.upper()+": "+ pc1 )
        ax.set_ylabel(l2.upper()+": "+ pc2 )
    else:
        ax.set_xlabel(l1.upper())
        ax.set_ylabel(l2.upper())
    ax.legend()
    #plt.tight_layout()
    return fig, ax
