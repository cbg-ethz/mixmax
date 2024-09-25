import argparse
import numpy as np
import logging
from pathlib import Path
from itertools import batched


import loompy

logging.basicConfig(format="{asctime} - {levelname} - {message}", style="{", datefmt="%Y-%m-%d %H:%M",level=logging.DEBUG)


def sample_split(alpha, seed = None):
    if seed:
        rng = np.random.default_rng(seed=seed)
    else:
        rng = np.random.default_rng()
    splitting_ratio = rng.beta(alpha, alpha, size=1)[0]

    return splitting_ratio

def split_loom_file(splitting_ratio, loom_file, output_dir, seed = None):
    if seed:
        rng = np.random.default_rng(seed=seed)
    else:
        rng = np.random.default_rng()
    with loompy.connect(loom_file) as ds:
        n_cells = ds.shape[1]
        split_vector = rng.choice([0,1], size=n_cells, p=[1-splitting_ratio, splitting_ratio])
        
        indices_first_sample = np.where(split_vector == 1)[0]
        indices_second_sample = np.where(split_vector == 0)[0]
        ds1 = ds.view[:, indices_first_sample]
        ds2 = ds.view[:, indices_second_sample]
        loom_file = Path(loom_file)
        (output_dir / 'temp').mkdir(parents=True, exist_ok=True)
        output_1 = (output_dir / 'temp' / f'{loom_file.stem}_split1.loom').as_posix()
        logging.info(f'Writing first split to {output_1}')
        loompy.create(output_1, ds1.layers, ds1.ra, ds1.ca)
        
        output_2 = (output_dir / 'temp' / f'{loom_file.stem}_split2.loom').as_posix()
        logging.info(f'Writing second split to {output_2}')
        loompy.create(output_2, ds2.layers, ds2.ra, ds2.ca)
        
    return None


def parse_args():
    parser = argparse.ArgumentParser(description="Simulate multiplexed single-cell (or multi-sample) mutation data")
    parser.add_argument("--input", type=str, required=True, help="Input list of loom files")
    parser.add_argument("--alpha", type=float, required=True, help="Dirichlet distribution parameter")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory")

    logging.info("Parsing arguments")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    def main(args):
        logging.info(f"Splitting file {args.input}")
        loom_file = args.input
        splitting_ratio = sample_split(alpha=args.alpha)
        split_loom_file(splitting_ratio = splitting_ratio, loom_file =  loom_file, output_dir = Path(args.output_dir))
        logging.info(f"Done splitting file {args.input}")