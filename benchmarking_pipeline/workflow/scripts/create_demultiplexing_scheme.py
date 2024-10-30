import argparse
import logging
import itertools
import yaml
from pathlib import Path

import numpy as np

logging.basicConfig(
    format='{asctime} - {levelname} - {message}',
    style='{',
    datefmt='%Y-%m-%d %H:%M',
    level=logging.INFO,
)


def define_demultiplexing_scheme_optimal_case(
    maximal_number_of_samples, maximal_pool_size, n_samples, robust
):
    if n_samples % maximal_number_of_samples != 0:
        raise ValueError(
            'Number of samples must be a multiple of maximal number of samples to run this function!'
        )
    demultiplexing_scheme = {}
    number_of_iterations = int(n_samples / maximal_number_of_samples)
    if not robust:
        unordered_unique_pairs = list(
            itertools.combinations(range(maximal_pool_size), 2)
        )
        diagonal = list(zip(range(maximal_pool_size), range(maximal_pool_size)))
        unordered_pairs = unordered_unique_pairs + diagonal
    else:
        unordered_pairs = list(itertools.combinations(range(maximal_pool_size + 1), 2))

    for idx1 in range(number_of_iterations):
        for idx2, pair in enumerate(unordered_pairs):
            demultiplexing_scheme[int(idx1 * maximal_number_of_samples + idx2 + 1)] = (
                pair + (idx1,)
            )  # naming of samples starts with 1

    return demultiplexing_scheme


def find_demultiplexing_scheme(maximal_pool_size, n_samples, robust):
    maximal_number_of_samples = (maximal_pool_size * (maximal_pool_size + 1)) / 2

    if n_samples % maximal_number_of_samples != 0:
        logging.error('So far, only defined for certain cohort sizes')
        raise NotImplementedError
    else:
        logging.info('Creating multiplexing scheme')
        demultiplexing_scheme = define_demultiplexing_scheme_optimal_case(
            maximal_number_of_samples=maximal_number_of_samples,
            maximal_pool_size=maximal_pool_size,
            n_samples=n_samples,
            robust=robust,
        )

        logging.info(f'Demultiplexing scheme: {demultiplexing_scheme}')

        return demultiplexing_scheme


def multiplexing_scheme_format2pool_format(demultiplexing_scheme):
    pool_scheme = {}
    total_no_pools_per_repetition = np.max(
        [pool[:-1] for pool in list(demultiplexing_scheme.values())]
    )
    total_no_of_repetitions = np.max(
        [pool[-1] for pool in list(demultiplexing_scheme.values())]
    )

    pools = list(
        itertools.product(
            range(total_no_pools_per_repetition + 1), range(total_no_of_repetitions + 1)
        )
    )
    for pool in pools:
        pool_scheme[pool] = []
    for key, value in demultiplexing_scheme.items():
        pool_scheme[value[0], value[-1]].append(int(key))
        pool_scheme[value[1], value[-1]].append(-int(key))
    # The pool scheme has pool identifier as keys and split libraries as values, where the first split is encoded as a positive integer and the second split as a negative integer

    return pool_scheme


def select_samples_for_pooling(pool_scheme, input_dir, sample_list):
    if isinstance(sample_list[0], str):
        for idx, sample in enumerate(sample_list):
            sample_list[idx] = Path(sample)
    if isinstance(input_dir, str):
        input_dir = Path(input_dir)
    logging.info('Selecting samples for pooling')
    pools_summary = {}
    for pool, samples in pool_scheme.items():
        logging.info(f'Retrieving list of samples to be pooled for pool {pool}')
        logging.debug(f'Samples: {samples}')
        logging.debug(f'Sample list: {sample_list}')
        loom_files = []
        for sample in samples:
            if sample > 0:
                loom_files.append(
                    (input_dir / f'{sample_list[sample-1].stem}_split1.loom').as_posix()
                )
            else:
                loom_files.append(
                    (
                        input_dir / f'{sample_list[-sample-1].stem}_split2.loom'
                    ).as_posix()
                )

        logging.info(f'Done retrieving list of samples to be pooled for pool {pool}')

        pools_summary[f'({pool[0]}.{pool[1]})'] = loom_files

    return pools_summary


def write_output(pools_summary, output):
    with open(output, 'w') as f:
        yaml.dump(pools_summary, f, default_flow_style=False)


def load_input_samples(input_sample_file):
    with open(input_sample_file, 'r') as f:
        input_samples = f.readlines()

    return input_samples


def parse_args():
    parser = argparse.ArgumentParser(
        description='Output a multiplexing scheme under pool size constraint for a predefined number of samples'
    )
    parser.add_argument(
        '--robust',
        type=bool,
        required=False,
        default=False,
        help='Input list of loom files',
    )
    parser.add_argument(
        '-k',
        '--maximal_pool_size',
        type=int,
        required=False,
        help='The maximal amount of samples to be multiplexed in a pool',
    )
    parser.add_argument(
        '--n_samples',
        type=int,
        required=True,
        help='The number of samples to be sequenced in the cohort',
    )
    parser.add_argument(
        '--output', type=str, required=True, help='Where to store the output file'
    )
    parser.add_argument(
        '--input_dir', type=str, required=True, help='Directory to the loom files'
    )
    parser.add_argument(
        '--input_sample_file',
        type=str,
        required=True,
        help='Path to the file containing the list of loom files',
    )

    logging.info('Parsing arguments')
    args = parser.parse_args()

    if args.n_samples == 0:
        logging.error('Number of samples must be greater than 0')
        raise ValueError

    return args


def main(args):
    input_samples = load_input_samples(args.input_sample_file)
    n_samples = len(input_samples)
    demultiplexing_scheme = find_demultiplexing_scheme(
        args.maximal_pool_size, n_samples, args.robust
    )
    pool_scheme = multiplexing_scheme_format2pool_format(demultiplexing_scheme)
    pools_summary = select_samples_for_pooling(
        pool_scheme, args.input_dir, input_samples
    )
    logging.debug(f'Output: {pools_summary.keys()}')
    logging.info(f'Writing output to {args.output}')
    write_output(pools_summary, args.output)
    logging.info('Success.')


if __name__ == '__main__':
    args = parse_args()
    main(args)
