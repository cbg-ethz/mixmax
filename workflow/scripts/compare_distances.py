import argparse

import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist, squareform, cityblock
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
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
    parser.add_argument("--tsv_files", nargs='+' , type=str, help="Path to the directory containing the TSV files")
    parser.add_argument("--pool_scheme", type=str, help="Path to the pool scheme file")
    return parser.parse_args()



args = parse_args()

# Read all TSV files into a list of DataFrames
dfs = [pd.read_csv(f, sep='\t', index_col = 0) for f in args.tsv_files]

# Get a set of all column names
all_columns = set()
for df in dfs:
    all_columns.update(df.columns)

# Ensure all DataFrames have the same columns, filling missing columns with -1
for df in dfs:
    for col in all_columns:
        if col not in df.columns:
            df[col] = np.nan


# Concatenate all DataFrames
concatenated_df = pd.concat(dfs[:2], ignore_index=True)


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


# Compute the distance matrix
for df in dfs:
    df.drop(columns=cols_to_remove, inplace=True)
    df.drop(columns=cols_to_remove_45_55, inplace=True)

    # Ensure all DataFrames have the same columns, filling missing columns with NaN
    for df in dfs:
        df.drop(columns=cols_to_remove_45_55, inplace=True)
# Compute the distance matrix considering only non-NaN entries and normalizing
def custom_distance(u, v):
    mask = ~np.isnan(u) & ~np.isnan(v)
    if np.sum(mask) == 0:
        return np.nan
    return np.sqrt(np.sum((u[mask] - v[mask]) ** 2) / np.sum(mask))

    # Compute the distance matrix for all pairs of DataFrames

distance_matrices = np.empty((len(dfs), len(dfs)), dtype=object)
for data1, df1 in enumerate(dfs):
    for data2, df2 in enumerate(dfs):
        distance_matrix = np.zeros((df1.shape[0], df2.shape[0]))
        for i in range(df1.shape[0]):
            for j in range(df2.shape[0]):
                distance_matrix[i, j] = custom_distance(df1.iloc[i].values, df2.iloc[j].values)
        distance_matrices[data1, data2] = distance_matrix

# Plot the heatmaps for each pair of DataFrames
fig, axes = plt.subplots(nrows=len(dfs), ncols=len(dfs), figsize=(20, 20))
for i, df1 in enumerate(dfs):
    for j, df2 in enumerate(dfs):
        sns.heatmap(distance_matrices[i,j], ax=axes[i, j], cmap=cmap, cbar=False)
        axes[i, j].set_title(f'Distance: DF{i+1} vs DF{j+1}')
plt.tight_layout()
plt.savefig("/cluster/work/bewi/members/jgawron/projects/Demultiplexing/pairwise_distances_heatmaps.png")

distance_matrix = pdist(concatenated_df.fillna(0), metric='euclidean')
distance_matrix = pdist(concatenated_df.values, metric=custom_distance)
distance_matrix = squareform(distance_matrix)

# Convert the distance matrix to a DataFrame for easier plotting
distance_df = pd.DataFrame(distance_matrix)

# Plot the heatmap
plt.figure(figsize=(10, 8))



# Mask the values that are -1
masked_distance_df = distance_df.mask(distance_df == -1)

# Plot the heatmap with the custom colormap
sns.heatmap(masked_distance_df, cmap=cmap, cbar_kws={'label': 'Distance'})
plt.title('Distance Matrix Heatmap')
plt.savefig("/cluster/work/bewi/members/jgawron/projects/Demultiplexing/prototype.png")


fig, axes = plt.subplots(nrows=1, ncols=len(dfs), figsize=(15, 5))

for i, df in enumerate(dfs):
    sns.heatmap(df.fillna(0), ax=axes[i], cmap=cmap, cbar=False)
    axes[i].set_title(f'Heatmap {i+1}')
plt.savefig("/cluster/work/bewi/members/jgawron/projects/Demultiplexing/vectors_heatmaps.png")

# Compute the matrix of counts of non-NaN pairs
non_nan_counts = np.zeros((concatenated_df.shape[0], concatenated_df.shape[0]))

for i in range(concatenated_df.shape[0]):
    for j in range(concatenated_df.shape[0]):
        non_nan_counts[i, j] = np.sum(~np.isnan(concatenated_df.iloc[i]) & ~np.isnan(concatenated_df.iloc[j]))

# Convert the counts matrix to a DataFrame for easier plotting
non_nan_counts_df = pd.DataFrame(non_nan_counts)

# Plot the heatmap of non-NaN pairs
plt.figure(figsize=(10, 8))
sns.heatmap(non_nan_counts_df, cmap='viridis', cbar_kws={'label': 'Number of Non-NaN Pairs'})
plt.title('Number of Non-NaN Pairs Heatmap')
plt.savefig("/cluster/work/bewi/members/jgawron/projects/Demultiplexing/non_nan_pairs_heatmap.png")


# Function to compute the ratio of the smallest and the second smallest value in a matrix
def compute_ratio(matrix):
    flattened = matrix.flatten()
    sorted_values = np.sort(flattened[~np.isnan(flattened)])
    if len(sorted_values) < 2:
        return np.nan
    return sorted_values[1] / sorted_values[0]

# Compute the ratios for all distance matrices
ratios = []
for i in range(len(dfs)):
    for j in range(len(dfs)):
        ratio = compute_ratio(distance_matrices[i, j])
        ratios.append((i, j, ratio))

# Sort the distance matrices by the computed ratio, from largest to smallest
sorted_ratios = sorted(ratios, key=lambda x: x[2], reverse=True)


# Initialize a list to store the pairs with the lowest values
lowest_value_pairs = []


used_samples = [[]*len(dfs) for _ in range(len(dfs))]

# Iterate through the sorted ratios, skipping diagonal matrices
for i, j, _ in sorted_ratios:
    if i <= j:
        continue  # Skip diagonal matrices

    # Get the current distance matrix
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
    print(f'Lowest value in Distance Matrix DF{i+1} vs DF{j+1}: Row = {row}, Column = {col}')




with open(args.pool_scheme, 'r') as file:
    pooling_scheme = yaml.safe_load(file)

demultiplexing_scheme = pool_format2multiplexing_scheme(pooling_scheme)    

print(demultiplexing_scheme)
for i, j, row, col in lowest_value_pairs:
    if i > j:
        i, j = j, i
        row, col = col, row
    print(f'Sample {row} from pool {i} and sample {col} from pool {j}: {demultiplexing_scheme[(i,j,0)]}')

for i in range(len(used_samples)):
    unused_sample = next(sample for sample in range(len(dfs[i])) if sample not in used_samples[i])
    print(f"Sample {unused_sample} from pool {i}: {demultiplexing_scheme[(i,i,0)]}")



# Save the concatenated DataFrame to a new TSV file
#output_file = '/cluster/work/bewi/members/jgawron/projects/Demultiplexing/concatenated.tsv'
#concatenated_df.to_csv(output_file, sep='\t', index=False)

#print(f"Concatenated TSV file saved to {output_file}")