import itertools
import csv

# Define the parameter values
cell_counts = [1000]
doublet_rates = [0.05]
pool_sizes = [4]
robust_values = [True]

# Generate all combinations of parameters
combinations = list(itertools.product(robust_values, pool_sizes, doublet_rates, cell_counts))

# Define the output file path
output_file = 'parameter_space.tsv'

# Write the combinations to a TSV file
with open(output_file, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    # Write the header
    writer.writerow(['robust', 'pool_size', 'doublet_rate', 'cell_count'])
    # Write the parameter combinations
    for combination in combinations:
        writer.writerow(combination)

print(f"Parameter space written to {output_file}")