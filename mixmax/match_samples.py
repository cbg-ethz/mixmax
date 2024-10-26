import argparse
import logging

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import yaml
from pathlib import Path



def pool_format2multiplexing_scheme(pool_scheme):
    demultiplexing_scheme = {}
    for pool, samples in pool_scheme.items():
        for sample in samples:
            if "_split1.loom" in sample:
                sample_name = str(Path(sample).stem)
                sample_name = sample_name.replace('_split1', '')
                number_of_iterations =pool.split('.')[1][:-1]
                pool_ID = pool.split('.')[0][1:]
                if sample_name not in demultiplexing_scheme.keys():
                    demultiplexing_scheme[sample_name] = [int(pool_ID), np.inf, int(number_of_iterations)]
                else:
                    demultiplexing_scheme[sample_name][0] = int(pool_ID)
                    demultiplexing_scheme[sample_name][2] = int(number_of_iterations)
            if "_split2.loom" in sample:
                sample_name = str(Path(sample).stem)
                sample_name = sample_name.replace('_split2', '')
                pool_ID = pool.split('.')[0][1:]
                if sample_name not in demultiplexing_scheme.keys():
                    demultiplexing_scheme[sample_name] = [np.inf, int(pool_ID), np.inf]
                else:
                    demultiplexing_scheme[sample_name][1] = int(pool_ID)
            # Swap keys and items in the dictionary
    swapped_demultiplexing_scheme = {tuple(v): k for k, v in demultiplexing_scheme.items()}
    
    return swapped_demultiplexing_scheme




# Create a custom colormap
cmap = sns.color_palette("viridis", as_cmap=True)
cmap.set_bad(color='grey')


def parse_args():
    parser = argparse.ArgumentParser(description="Compare distances between samples")
    parser.add_argument("--tsv_files", nargs='+' , type=str, help="Path to TSV files containing the genotypes of the clusters. Must be in the same order as the pools in the pooling scheme.")
    parser.add_argument("--pool_scheme", type=str, help="Path to the pool scheme file")
    parser.add_argument("--output_plot", type=str, help="Output heatmap")
    parser.add_argument("--output", type=str, help="sample assignment")
    parser.add_argument("--robust", type=bool, required=False, default = False, help="Robust assignment")
    return parser.parse_args()


def custom_distance(u, v):
    mask = ~np.isnan(u) & ~np.isnan(v)
    if np.sum(mask) == 0:
        return np.nan
    return np.sqrt(np.sum((u[mask] - v[mask]) ** 2) / np.sum(mask))


def compute_ratio(matrix):
    flattened = matrix.flatten()
    sorted_values = np.sort(flattened[~np.isnan(flattened)])
    if len(sorted_values) < 2:
        return np.nan
    return sorted_values[1] / sorted_values[0]


def permute_tsv_files(tsv_files, pool_scheme):
    pools = [f"({Path(file).name.split("_")[1]})" for file in tsv_files]
    
    pool_scheme_keys = list(pool_scheme.keys())
    pool_permutation = [pools.index(pool) for pool in pool_scheme_keys]
    
    if not all((pool1 == pool2 for pool1, pool2 in zip([pools[permuted_idx] for permuted_idx in pool_permutation], pool_scheme_keys))):
        raise ValueError("The pools in the pool scheme do not match the pools in the TSV files")

    return pool_permutation

args = parse_args()

def load_data(args):

    with open(args.pool_scheme, 'r') as file:
        pooling_scheme = yaml.safe_load(file)

    permutation_of_pools = permute_tsv_files(args.tsv_files, pooling_scheme)

    if isinstance(args.tsv_files, str):
        tsv_files = [args.tsv_files]
    tsv_files = list(args.tsv_files)
    tsv_files_permuted = [tsv_files[i] for i in permutation_of_pools]
    
    dfs = [pd.read_csv(f, sep='\t', index_col = 0) for f in tsv_files_permuted]

    return dfs, pooling_scheme

dfs, pooling_scheme = load_data(args)

# Get a set of all column names
all_columns = set()
for df in dfs:
    all_columns.update(df.columns)

print(all_columns)


for idx,df in enumerate(dfs):
    for col in all_columns:
        if col not in df.columns:
            dfs[idx][col] = np.nan

for idx,df in enumerate(dfs):
    dfs[idx] = df[sorted(all_columns)]
    
# Concatenate all DataFrames
concatenated_df = pd.concat(dfs, ignore_index=True)[sorted(all_columns)]


print(concatenated_df)
# Remove columns which are above 0.9 in all row entries and below 0.1 in all row entries

cols_to_remove = concatenated_df.columns[
    (concatenated_df > 0.9).all(axis=0) | (concatenated_df < 0.1).all(axis=0)
]
# Remove columns which have entries between 0.45 and 0.55 for all row entries
cols_to_remove_45_55 = concatenated_df.columns[
    ((concatenated_df >= 0.45) & (concatenated_df <= 0.55)).all(axis=0)
]

concatenated_df.drop(columns=cols_to_remove, inplace=True)
concatenated_df.drop(columns=cols_to_remove_45_55, inplace=True)

for idx,df in enumerate(dfs):
    dfs[idx] = df.drop(columns=cols_to_remove, inplace=False)
    dfs[idx] = df.drop(columns=cols_to_remove_45_55, inplace=False)


# Compute the distance matrix
for df in dfs:
    df.drop(columns=cols_to_remove, inplace=True)
    df.drop(columns=cols_to_remove_45_55, inplace=True)

    # Ensure all DataFrames have the same columns, filling missing columns with NaN
    for df in dfs:
        df.drop(columns=cols_to_remove_45_55, inplace=True)
# Compute the distance matrix considering only non-NaN entries and normalizing


    # Compute the distance matrix for all pairs of DataFrames

distance_matrices = np.empty((len(dfs), len(dfs)), dtype=object)
for data1, df1 in enumerate(dfs):
    for data2, df2 in enumerate(dfs):
        distance_matrix = np.zeros((df1.shape[0], df2.shape[0]))
        for i in range(df1.shape[0]):
            for j in range(df2.shape[0]):
                distance_matrix[i, j] = custom_distance(df1.iloc[i].values, df2.iloc[j].values)
        distance_matrices[data1, data2] = distance_matrix


remove = None
""" if args.robust == True:
    for i in range(len(dfs)):
        distance_matrix = distance_matrices[i, i]
        min_off_diagonal = np.min(distance_matrix + np.eye(distance_matrix.shape[0]) * np.inf)
        minimal_value = np.inf
        logging.warning(min_off_diagonal)
        if (min_off_diagonal < 10) and (min_off_diagonal < minimal_value):
            remove = i
            minimal_value = min_off_diagonal

    if remove is not None:
    #    dfs.pop(remove)
    #    distance_matrices = np.delete(distance_matrices, remove, axis=0)
    #    distance_matrices = np.delete(distance_matrices, remove, axis=1)
        args.robust = False
        logging.warning(f"Discarding dataset {remove} and going into non-robust mode") """

# Plot the heatmaps for each pair of DataFrames
fig, axes = plt.subplots(nrows=len(dfs), ncols=len(dfs), figsize=(20, 20))
for i, df1 in enumerate(dfs):
    for j, df2 in enumerate(dfs):
        """ if remove is not None:
            if i == remove or j == remove:
                continue
        """
        sns.heatmap(distance_matrices[i, j], ax=axes[i, j], cmap=cmap, cbar=False)
        axes[i, j].set_title(f'Distance: DF{i+1} vs DF{j+1}')
plt.tight_layout()

heatmap_plot = Path(args.output_plot)
genotype_plot = heatmap_plot.parent / 'genotype_heatmap.png'


plt.savefig(heatmap_plot)



fig, axes = plt.subplots(nrows=1, ncols=len(dfs), figsize=(15, 5))

for i, df in enumerate(dfs):
    sns.heatmap(df.fillna(0), ax=axes[i], cmap=cmap, cbar=False)
    axes[i].set_title(f'Pool {i}')
plt.savefig(genotype_plot)





ratios = []
for i in range(len(dfs)):
    for j in range(len(dfs)):
        ratio = compute_ratio(distance_matrices[i, j])
        ratios.append((i, j, ratio))

sorted_ratios = sorted(ratios, key=lambda x: x[2], reverse=True)


lowest_value_pairs = []


used_samples = [[]*len(dfs) for _ in range(len(dfs))]
""" if remove is not None:
    for j in range(len(dfs)):
        if j != remove:
            used_samples[remove].extend(range(len(dfs[remove])))
            used_samples[j].extend(range(len(dfs[j]))) """

for i, j, _ in sorted_ratios:
    if i <= j:
        continue
    """ if remove is not None:
        if i == remove or j == remove:
            continue """

    distance_matrix = distance_matrices[i, j]

    # Find the pair (i, j) with the lowest value that hasn't been used yet
    min_value = np.inf
    min_pair = None
    for row in range(distance_matrix.shape[0]):
        if row in used_samples[i]:
            continue
        for col in range(distance_matrix.shape[1]):
            if col in used_samples[j]:
                continue
            if distance_matrix[row, col] < min_value:
                min_value = distance_matrix[row, col]
                min_pair = (row, col)

    if min_pair is not None:
        # Save the pair with the lowest value
        lowest_value_pairs.append((i, j, min_pair[0], min_pair[1]))
        # Mark the row and column as used
        used_samples[i].append(min_pair[0])
        used_samples[j].append(min_pair[1])

        # Recompute the ratios for all unprocessed matrices
        ratios = []
        for x in range(len(dfs)):
            for y in range(len(dfs)):
                if x <= y or (x, y) in [(i, j) for i, j, _, _ in lowest_value_pairs]:
                    continue  # Skip diagonal and already processed matrices
                ratio = compute_ratio(distance_matrices[x, y])
                ratios.append((x, y, ratio))

        # Sort the distance matrices by the recomputed ratio, from largest to smallest
        sorted_ratios = sorted(ratios, key=lambda x: x[2], reverse=True)
    else:
        print(f"skipping assignment in matrix {i},{j}")

# Print the pairs with the lowest values
for i, j, row, col in lowest_value_pairs:
    print(f'Lowest value in Distance Matrix DF{i} vs DF{j}: Row = {row}, Column = {col}')



demultiplexing_scheme = pool_format2multiplexing_scheme(pooling_scheme)    

for i, j, row, col in lowest_value_pairs:
    if i > j:
        i, j = j, i
        row, col = col, row
    print(f'Sample {row} from pool {i} and sample {col} from pool {j}: {demultiplexing_scheme[(i,j,0)]}')

if args.robust == False: 
    for i in range(len(used_samples)):
        unused_sample = next(sample for sample in range(len(dfs[i])) if sample not in used_samples[i])
        print(f"Sample {unused_sample} from pool {i}: {demultiplexing_scheme[(i,i,0)]}")



pooling_scheme_keys = list(pooling_scheme.keys())
sample_assignment = {}
for i, j, row, col in lowest_value_pairs:
    if i > j:
        i, j = j, i
        row, col = col, row
    if pooling_scheme_keys[i] not in sample_assignment.keys():
        sample_assignment[pooling_scheme_keys[i]] = {row: demultiplexing_scheme[(i, j, 0)]}
    else:
        sample_assignment[pooling_scheme_keys[i]][row] = demultiplexing_scheme[(i, j, 0)]
    if pooling_scheme_keys[j] not in sample_assignment.keys():
        sample_assignment[pooling_scheme_keys[j]] = {col: demultiplexing_scheme[(i, j, 0)]}
    else:
        sample_assignment[pooling_scheme_keys[j]][col] = demultiplexing_scheme[(i, j, 0)]

if args.robust == False: 
    for i in range(len(used_samples)):
        unused_sample = next(sample for sample in range(len(dfs[i])) if sample not in used_samples[i])
        if pooling_scheme_keys[i] not in sample_assignment.keys():
            sample_assignment[pooling_scheme_keys[i]] = {unused_sample: demultiplexing_scheme[(i, i, 0)]}
        else:
            sample_assignment[pooling_scheme_keys[i]][unused_sample] = demultiplexing_scheme[(i, i, 0)]

with open(args.output, 'w') as yaml_file:
    yaml.dump(sample_assignment, yaml_file)
