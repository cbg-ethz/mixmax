import argparse
import logging
import itertools
import yaml

logging.basicConfig(format="{asctime} - {levelname} - {message}", style="{", datefmt="%Y-%m-%d %H:%M",level=logging.DEBUG)


def define_demultiplexing_scheme_optimal_case(maximal_number_of_samples, maximal_pool_size, n_samples):
    demultiplexing_scheme = {}
    number_of_iterations = int(n_samples / maximal_number_of_samples)
    unordered_pairs_unique_pairs = list(itertools.combinations(range(maximal_pool_size), 2))
    diagonal = list(zip(range(maximal_pool_size), range(maximal_pool_size)))
    unordered_pairs = unordered_pairs_unique_pairs + diagonal                           

    for idx1 in range(number_of_iterations):
        for idx2, pair in enumerate(unordered_pairs):
            demultiplexing_scheme[int(idx1 * maximal_number_of_samples + idx2 + 1)] = pair + (idx1,) # naming of samples starts with 1
    return demultiplexing_scheme


def find_demultiplexing_scheme(args):
    maximal_number_of_samples = (args.maximal_pool_size * (args.maximal_pool_size+1))/2

    if args.n_samples % maximal_number_of_samples != 0:
        logging.error("So far, only defined for certain cohort sizes")
        raise NotImplementedError
    else:
        logging.info("Creating multiplexing scheme")
        demultiplexing_scheme = define_demultiplexing_scheme_optimal_case(maximal_number_of_samples = maximal_number_of_samples, maximal_pool_size = args.maximal_pool_size, n_samples = args.n_samples)
        logging.info(f'Demultiplexing scheme: {demultiplexing_scheme}')
        return demultiplexing_scheme

def multiplexing_scheme_format2pool_format(demultiplexing_scheme):
    pool_scheme = {}
    total_no_pools_per_repetition = np.max([pool[:-1] for pool in list(demultiplexing_scheme.values())])
    total_no_of_repetitions = np.max([pool[-1] for pool in list(demultiplexing_scheme.values())])

    pools = list(itertools.product(range(total_no_pools_per_repetition + 1), range(total_no_of_repetitions + 1)))
    for pool in pools:
        pool_scheme[pool] = []
    for key, value in demultiplexing_scheme.items():
        pool_scheme[value[0],value[-1]].append(int(key))
        pool_scheme[value[1],value[-1]].append(-int(key))
    # The pool scheme has pool identifier as keys and split libraries as values, where the first split is encoded as a positive integer and the second split as a negative integer

    return pool_scheme


def select_samples_for_pooling(pool_scheme, args):
    pools_summary = {}
    for pool, samples in pool_scheme.items():
    
        logging.info(f"Retrieving list of samples to be pooled for pool {pool}")
    
        loom_files = []
        for sample in samples:
            if sample > 0:
                loom_files.append(args.split_loom_files / sample_list[sample].stem + "_split1.loom")
            else:
                loom_files.append(args.split_loom_files / sample_list[-sample].stem + "_split2.loom")
    
        logging.info(f"Done retrieving list of samples to be pooled for pool {pool}")

        pools_summary[pool] = loom_files

    return pools_summary


def write_output(pools_summary, output):
    with open(output, "w") as f:
        yaml.dump(pools_summary, f, default_flow_style=False)


def parse_args():
    parser = argparse.ArgumentParser(description="Output a multiplexing scheme under pool size constraint for a predefined number of samples")
    parser.add_argument("--robust", type=bool, required=False, default = False, help="Input list of loom files")
    parser.add_argument("-k", "--maximal_pool_size", type=int, required=False, help="The maximal amount of samples to be multiplexed in a pool")
    parser.add_argument("--n_samples", type=int, required=True, help="The number of samples to be sequenced in the cohort")
    parser.add_argument("--output", type=str, required=True, help="Where to store the output file")
    parser.add_argument("--split_loom_files", type=str, required=True, help="Directory to the loom files")

    logging.info("Parsing arguments")
    args = parser.parse_args()
    
    if args.n_samples == 0:
        logging.error("Number of samples must be greater than 0")
        raise ValueError
    
    return args


def main():
    #sanity-checks
    demultiplexing_scheme = {1: (0, 1, 0), 2: (0, 2, 0), 3: (1, 2, 0), 4: (0, 0, 0), 5: (1, 1, 0), 6: (2, 2, 0), 7: (0, 1, 1), 8: (0, 2, 1), 9: (1, 2, 1), 10: (0, 0, 1), 11: (1, 1, 1), 12: (2, 2, 1)}
    if not multiplexing_scheme_format2pool_format(demultiplexing_scheme) == {(0, 0): [1,2, 4,-4], (1, 0): [-1,3,5,-5], (2, 0): [-2, -3, 6, -6], (0, 1): [7, 8, 10, -10], (1,1): [-7, 9, 11, -11], (2, 1): [-8, -9, 12, -12]}:
        raise ValueError("Bug in function multiplexing_scheme_format2pool_format")


    args = parse_args()
    demultiplexing_scheme = find_demultiplexing_scheme(args)
    pool_scheme = multiplexing_scheme_format2pool_format(demultiplexing_scheme)
    pools_summary = select_samples_for_pooling(pool_scheme, args)
    write_output(pools_summary, args.output)
