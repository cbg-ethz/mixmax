import argparse
import random
from pathlib import Path
import logging

logging.basicConfig(format="{asctime} - {levelname} - {message}", style="{", datefmt="%Y-%m-%d %H:%M",level=logging.DEBUG)

def generate_input_samples(loom_files, number_of_samples):
    return random.sample(loom_files, number_of_samples)




def parse_args():
    parser = argparse.ArgumentParser(description="Get a random list of samples to be pooled")
    parser.add_argument("--loom_file_path", type=str, help="Path to a folder containing loom files")
    parser.add_argument("--output", type=str, default=None, help="path to a directory where the sample list will be stored")
    parser.add_argument("--seed", type=int, required=True, default=None, help="random seed for reproducibility of sampling")
    parser.add_argument("--number_of_samples", type=int, default=1, help="number of samples to be selected")

    return parser.parse_args()


def main(args):
    loom_files = []
    for p in Path(args.loom_file_path).iterdir():
        if p.is_file() and p.suffix == ".loom":
            loom_files.append(p)
            
    if args.number_of_samples > len(loom_files):
        logging.error(f"Number of samples ({number_of_samples}) is greater than the number of loom files ({len(loom_files)}). Continuing with maximal number of samples.")
        input_samples = loom_files
        number_of_samples = len(loom_files)
    else:
        random.seed(args.seed)
        input_samples = random.sample(loom_files, args.number_of_samples)

    
    with open(args.output, "w") as f:
        for sample in input_samples:
            f.write(f"{sample}\n")


if __name__ == "__main__":
    args = parse_args()
    main(args)