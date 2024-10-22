from pathlib import Path
import yaml
import logging
import random

import numpy as np
import pandas as pd
from snakemake.utils import Paramspace

paramspace = Paramspace(pd.read_csv(Path(workflow.basedir) / 'sandbox' / 'parameter_space.tsv', sep='\t'))

output_folder = Path(config["output_folder"])
loom_file_path = Path(config["loom_files"])
number_of_samples = config.get("n_samples", 10)
seed=config["seed"]
doublet_rate = config["doublet_rate"]
cell_count = config["downsampling"]

loom_files = []
for p in loom_file_path.iterdir():
    if p.is_file() and p.suffix == ".loom":
        loom_files.append(p)




def get_variant_list(pool_ID):
    file = output_folder / f"pool_{pool_ID}.txt"
    with open(file, 'r') as f:
        return [line.strip() for line in f]


def get_split_files(path):
    if isinstance(path, str):
        path = Path(path)
    split_files = []
    with open(path, 'r') as f:
        for line in f:
            path2 = Path(line.strip())
            split_files.append(path.parent / "loom" / f"{path2.stem}_split1.loom")
            split_files.append(path.parent / "loom" / f"{path2.stem}_split2.loom")

    return split_files



def sample_mixing_ratios(seed, max_pool_size):
    concentrations = [40/max_pool_size] * max_pool_size
    rng = np.random.default_rng(seed=int(seed))
    ratios = rng.dirichlet(concentrations)
    
    return " ".join(str(item) for item in ratios)


def get_demultiplexed_samples(wildcards):
    pools = yaml.safe_load(open(checkpoints.create_demultiplexing_scheme.get(seed=wildcard.seed).output, 'r'))
    pools = list(pools.keys())
    demultiplexed_samples = [f"pool_{pool}_{wildcard.seed}_demultiplexed.assignments.tsv" for pool in pools]

        
def determine_number_of_different_donors(pool_file, pool_ID):
    samples = yaml.safe_load(open(pool_file, 'r'))[pool_ID]
    samples = [Path(sample) for sample in samples]
    unique_patterns = set()
    for sample in samples:
        if sample.suffix in [".loom"]:
            if "_split1" in sample.stem:
                unique_patterns.add(sample.stem.replace("_split1", ""))
            elif "_split2" in sample.stem:
                unique_patterns.add(sample.stem.replace("_split2", ""))
                
    return len(unique_patterns)


def get_all_demultiplexed_assignments(wildcards, wildcard_pattern):
    demultiplexing_scheme = yaml.safe_load(open(checkpoints.create_demultiplexing_scheme.get(seed=wildcards.seed, pool_size = wildcards.pool_size, robust = wildcards.robust).output.pools, 'r'))
    pool_names = list(demultiplexing_scheme.keys())
    pool_names = [pool_name.strip("()") for pool_name in pool_names]
    
    demultiplexed_files = []
    file_path = output_folder / f"seed_{wildcards.seed}" / wildcard_pattern / 'demultiplexed' 
    for pool in pool_names:
        demultiplexed_files.append(file_path / f"pool_{pool}_{wildcards.seed}_demultiplexed.assignments.tsv")

    return demultiplexed_files


def get_all_demultiplexed_profiles(wildcards, wildcard_pattern):
    demultiplexing_scheme = yaml.safe_load(open(checkpoints.create_demultiplexing_scheme.get(seed=wildcards.seed, pool_size = wildcards.pool_size, robust = wildcards.robust).output.pools, 'r'))
    pool_names = list(demultiplexing_scheme.keys())
    pool_names = [pool_name.strip("()") for pool_name in pool_names]
    
    demultiplexed_files = []
    file_path = output_folder / f"seed_{wildcards.seed}" / wildcard_pattern / 'demultiplexed' 
    for pool in pool_names:
        demultiplexed_files.append(file_path / f"pool_{pool}_{wildcards.seed}_demultiplexed.profiles.tsv")

    return demultiplexed_files



output_files = expand(output_folder / f"seed_{seed}" / "{params}" / 'sample_assignment.yaml', params=paramspace.instance_patterns) + expand(output_folder / f"seed_{seed}" / "{params}" / 'sample_identity.yaml', params=paramspace.instance_patterns)
