#!/usr/bin/env python3

import argparse
import copy
import os
import re
import tkinter as tk

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, cityblock
from scipy.cluster.hierarchy import linkage
import seaborn as sns

VCF_COLS = [
    'CHROM',
    'POS',
    'ID',
    'REF',
    'ALT',
    'QUAL',
    'FILTER',
    'INFO',
    'FORMAT',
    'sample',
]
GT_MAP = {'0/0': 0, '0/1': 1, '1/1': 2, './.': np.nan}
DEF_DIST = 'Euclidean'
PID2SID = {
    'PIDasdf': 'SIDjklö'
}  ### These can be adapted manually to have more control on how the sampled should be named in the output
SID2PAPER = {'SIDjklö': 'S1'}

FILE_EXT = '.png'
FONTSIZE = 30
DPI = 300
sns.set_style('whitegrid')  # darkgrid, whitegrid, dark, white, ticks
sns.set_context(
    'paper',
    rc={
        'lines.linewidth': 2,
        'axes.axisbelow': True,
        'font.size': FONTSIZE,
        'axes.labelsize': 'medium',
        'axes.titlesize': 'large',
        'xtick.labelsize': 'medium',
        'ytick.labelsize': 'medium',
        'legend.fontsize': 'medium',
        'legend.title_fontsize': 'large',
        'axes.labelticksize': 50,
    },
)
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 3
plt.rcParams['xtick.bottom'] = True


def dist_nan(u, v):
    valid = ~np.isnan(u) & ~np.isnan(v)
    # return euclidean(u[valid], v[valid]) / np.sqrt(np.sum(2**2 * valid.sum())) # Euclidean
    # return (u[valid].round() != v[valid].round()).mean() # Hamming
    return cityblock(u[valid], v[valid]) / (2 * valid.sum())  # Manhattan


def get_clone_profiles(in_file):
    df = pd.read_csv(in_file)
    df['REF'].fillna('*', inplace=True)
    df['ALT'].fillna('*', inplace=True)
    idx_str = (
        df['CHR'].astype(str)
        + '_'
        + df['POS'].astype(str)
        + '_'
        + df['REF']
        + '_'
        + df['ALT']
    )
    idx_str = idx_str.str.replace('chr', '')
    gt = df.iloc[:, 7:].map(lambda x: int(x.split(':')[-1]))
    gt.replace(3, np.nan, inplace=True)

    df = pd.DataFrame(gt.mean(axis=1).values, index=idx_str, columns=['Cluster'])
    return df


def plot_heatmap(data, out_file):
    cmap = plt.get_cmap('YlOrRd', 100)
    dist_row = pdist(data.values, dist_nan)
    dist_col = pdist(data.values.T, dist_nan)
    Z_row = linkage(np.nan_to_num(dist_row, 10), 'ward')
    Z_col = linkage(np.nan_to_num(dist_col, 10), 'ward')
    try:
        cm = sns.clustermap(
            data,
            row_linkage=Z_row,
            col_linkage=Z_col,
            row_cluster=True,
            col_cluster=True,
            col_colors=None,
            row_colors=None,
            vmin=0,
            vmax=2,
            cmap=cmap,
            dendrogram_ratio=(0.1, 0.1),
            figsize=(25, 15),
            cbar_kws={
                'ticks': [0, 1, 2],
            },
            # cbar_pos=(0, 0, 0.0, 0.0),
            tree_kws=dict(linewidths=2),
        )
    except tk.TclError:
        import matplotlib

        matplotlib.use('Agg')
        cm = sns.clustermap(
            data,
            row_linkage=Z_row,
            col_linkage=Z_col,
            row_cluster=True,
            col_cluster=True,
            col_colors=None,
            row_colors=None,
            vmin=0,
            vmax=2,
            cmap=cmap,
            dendrogram_ratio=(0.1, 0.1),
            figsize=(25, 15),
            cbar_kws={
                'ticks': [0, 1, 2],
            },
            # cbar_pos=(0, 0, 0.0, 0.0),
            tree_kws=dict(linewidths=2),
        )

    cm.ax_heatmap.set_facecolor('#EAEAEA')
    cm.ax_heatmap.set_ylabel('\nProfiles')
    cm.ax_heatmap.set_xlabel('SNPs')
    labels_pretty = [
        'chr{}:{} {}>{}'.format(*i.get_text().split('_'))
        for i in cm.ax_heatmap.get_xticklabels()
    ]
    cm.ax_heatmap.set_xticklabels(labels_pretty, rotation=45, ha='right', va='top')

    cm.ax_cbar.set_title('Genotype')
    cm.ax_cbar.set_yticklabels(['0|0', '0|1', '1|1'])

    cm.fig.tight_layout()
    if out_file:
        cm.fig.savefig(out_file, dpi=DPI)
    else:
        plt.show()


def plot_distance(df, out_file):
    col_no = 1
    row_no = 1
    fig, ax = plt.subplots(
        nrows=row_no, ncols=col_no, figsize=(col_no * 12, row_no * 12)
    )

    sns.set(font_scale=2.5)
    font_size = 30

    df.sort_index(ascending=False, inplace=True)
    df = df[sorted(df.columns)]

    df.columns = [i.split('_')[0] for i in df.columns]
    df.index = [i.split('_')[0] for i in df.index]

    sns.heatmap(
        df,
        annot=True,
        square=True,
        cmap='viridis_r',
        cbar_kws={'shrink': 0.5, 'label': 'Distance'},
        linewidths=0,
        ax=ax,
    )
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=font_size)
    ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=font_size)
    ax.set_xlabel(df.index.name, fontsize=font_size + 10)
    ax.set_ylabel(df.columns.name, fontsize=font_size + 10)

    fig.tight_layout()
    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()


def assign_clusters(df_in, assignment_output):
    df = copy.deepcopy(df_in)
    # import pdb;pdb.set_trace()
    for i in range(df_in.shape[0]):
        cl_idx, other_idx = np.where(df == df.min().min())
        cl = df.index[cl_idx[0]]
        other = df.columns[other_idx[0]]
        print(f'Assigning: {cl} -> {other}')
        with open(assignment_output, 'a') as file:
            file.write(f'{cl}: {other}\n')

        df.drop(cl, axis=0, inplace=True)
        df.drop(other, axis=1, inplace=True)


def main(args):
    if not args.outfile:
        out_base = os.path.splitext(args.input[0])[0]
        out_end = FILE_EXT
    else:
        out_base, out_end = os.path.splitext(args.outfile)

    names = []
    for i, in_file in enumerate(args.input):

        cluster = re.compile(r'([0-9]*)_variants.csv').search(in_file).group(1)
        name = f'{cluster}_DNA'

        names.append(name)
        df_new = get_clone_profiles(in_file)
        df_new.rename({'Cluster': name}, inplace=True, axis=1)

        match_df = df_new.merge(SNP_df, left_index=True, right_index=True, how='inner')

        if match_df.size == 0:
            print('No overlap in RNA and DNA data.')
            exit()
        try:
            df = df.merge(
                match_df[name], left_index=True, right_index=True, how='outer'
            )
        except NameError:
            df = match_df

    df = df[sorted(df.columns)]
    # remove SNPs that are similar in all DNA clusters
    df = df[(df[names] > 1.95).sum(axis=1) != len(names)]
    df = df[(df[names] < 0.05).sum(axis=1) != len(names)]
    df = df[
        ((df[names] > 0.95).sum(axis=1) != len(names))
        | ((df[names] < 1.05).sum(axis=1) != len(names))
    ]
    # Remove SNPS not called in the scDNA-seq
    df = df[(df[names] >= 0).sum(axis=1) > 0]

    dists = []
    for name in names:
        dists.append(
            np.apply_along_axis(
                dist_nan, 1, df[SNP_df.columns].values.T, df[name].values.T
            )
        )

    dist_df = pd.DataFrame(dists, index=names, columns=SNP_df.columns)

    assign_clusters(dist_df, args.assignment_output)
    print(f'\n{dist_df}\n')

    dist_df.index.name = 'scDNA-seq'
    dist_df.columns.name = label

    dist_out = f'{out_base}_distance{out_end}'
    print(f'Generating plot: {dist_out}')
    plot_distance(dist_df, dist_out)

    hm_out = f'{out_base}_heatmap{out_end}'
    print(f'Generating plot: {hm_out}')
    plot_heatmap(df.T, hm_out)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input', type=str, nargs='+', help='mutation_profiles of the demultiplexed samples'
    )
    parser.add_argument(
        '-o', '--outfile', type=str, default='', help='Output file for heatmap.'
    )
    parser.add_argument(
        '-ao',
        '--assignment_output',
        type=str,
        help='Output yaml-file in which the automated sample matching is saved.',
    )

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)