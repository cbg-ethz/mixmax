from pathlib import Path

import loompy
import pandas as pd
import numpy as np
from scipy.stats import poisson

def process_loom_files(folder_path):
    folder = Path(folder_path)
    files = []
    cell_counts = []
    for filename in Path(folder).rglob("*.loom"):
        with loompy.connect(str(filename)) as ds:
            cell_counts.append(ds.shape[1])
            files.append(filename)
    data = pd.DataFrame({'cell_count': cell_counts}, index = files)
    import matplotlib.pyplot as plt

    # Plot the histogram of cell counts
    plt.hist(cell_counts, bins=30, density=True, alpha=0.6, color='g', label='Cell count histogram')

    # Plot the Poisson distribution density

    # Fit a Normal distribution to the cell counts
    mean_count = np.mean(cell_counts)
    std_dev = np.std(cell_counts)
 
    x = np.arange(0, max(cell_counts) + 1)
    from scipy.stats import norm
    plt.plot(x, norm.pdf(x, mean_count, std_dev), 'r-', label='Normal fit')
    
    plt.xlabel('Cell count')
    plt.ylabel('Density')
    plt.title('Histogram of Cell Counts with Poisson Distribution Fit')
    plt.legend()
    plt.savefig(folder / 'cell_count_distribution.png')

folder_path = '/cluster/work/bewi/members/jgawron/projects/Demultiplexing/AML_data/loom_files_complete'
process_loom_files(folder_path)