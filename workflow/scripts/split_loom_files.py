import argparse
import numpy as np
import logging
from pathlib import Path

import loompy

logging.basicConfig(format="{asctime} - {levelname} - {message}", style="{", datefmt="%Y-%m-%d %H:%M",level=logging.DEBUG)

def split_loom_file(seed, alpha, loom_file, output_dir):

    rng = np.random.default_rng(seed=seed)
    splitting_ratio = rng.beta(alpha, alpha, size=1)[0]
    with loompy.connect(loom_file) as ds:
        n_cells = ds.shape[1]
        split_vector = rng.choice([0,1], size=n_cells, p=[1-splitting_ratio, splitting_ratio])
        indices_first_sample = np.where(split_vector == 1)[0]
        indices_second_sample = np.where(split_vector == 0)[0]
        ds1 = ds.view[:, indices_first_sample]
        ds2 = ds.view[:, indices_second_sample]
        loompy.create(f'{str(output_dir)}/temp/{str(loom_file.stem)}_split1.loom', ds1.layers, ds1.ra, ds1.ca)
        loompy.create(f'{str(output_dir)}/temp{str(loom_file.stem)}_split2.loom', ds2.layers, ds2.ra, ds2.ca)
        
        return None


def parse_args():
    parser = argparse.ArgumentParser(description="Simulate multiplexed single-cell (or multi-sample) mutation data")
    parser.add_argument("--input", type=str, required=True, nargs='+', help="Input list of loom files")
    parser.add_argument("--seed", type=int, required=False, help="Set a seed for sampling splits")
    parser.add_argument("--alpha", type=float, required=True, help="Dirichlet distribution parameter")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory")

    logging.info("Parsing arguments")
    args = parser.parse_args()
    return args


args = parse_args()

for idx,loom_file in enumerate(args.input):
    logging.info(f"Splitting file {idx}")
    split_loom_file(seed = args.seed + idx, alpha = args.alpha, loom_file =  Path(loom_file), output_dir = Path(args.output_dir))
    logging.info(f"Done splitting file {idx}")