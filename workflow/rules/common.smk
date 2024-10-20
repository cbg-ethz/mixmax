from pathlib import Path
import yaml
import logging
import random

import numpy as np

#logging.basicConfig(format="{asctime} - {levelname} - {message}", style="{", datefmt="%Y-%m-%d %H:%M",level=logging.INFO)

output_folder = Path(config["output_folder"])
loom_file_path = Path(config["loom_files"])
number_of_samples = config.get("n_samples", 10)
seed=config["seed"]

loom_files = []
for p in loom_file_path.iterdir():
    if p.is_file() and p.suffix == ".loom":
        loom_files.append(p)



def generate_input_samples(loom_files, number_of_samples):
    random.sample(loom_files, number_of_samples)

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


def get_demultiplexed_samples(pool_file):
    pools = yaml.safe_load(open(checkpoints.create_demultiplexing_scheme.get(seed=w.seed).output, 'r'))
    samples = []

        
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
                
    return len(unique_patterns)+1

    
if number_of_samples > len(loom_files):
    logging.error(f"Number of samples ({number_of_samples}) is greater than the number of loom files ({len(loom_files)}). Continuing with maximal number of samples.")
    input_samples = loom_files
    number_of_samples = len(loom_files)
else:
    input_samples = random.sample(loom_files, number_of_samples)

input_samples = generate_input_samples(loom_files, number_of_samples)

#output_files.append(output_folder / f"seed_{seed}" / f"pools_{seed}.txt")
output_files = [output_folder / f"seed_{seed}" / "demultiplexed" / f"pool_{pool_ID}_{seed}_ground_truth_assignment.tsv" for pool_ID in ["0.0", "1.0", "2.0", "3.0"]]
output_files