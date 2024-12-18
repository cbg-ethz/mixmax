import argparse
import logging
from pathlib import Path

import yaml
import pandas as pd
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)


def parse_args():
    parser = argparse.ArgumentParser(description='Assign subsample simulation')
    parser.add_argument(
        '--demultiplexing_assignment',
        nargs='+',
        type=str,
        help='A tsv file that specifies which cell is assigned to which cluster',
    )
    parser.add_argument(
        '--output', type=str, help='Output file for the sample assignment'
    )

    return parser.parse_args()


def analyse_demultiplexing_assignment(assignment_file, pool_name):
    data = pd.read_csv(assignment_file, sep='\t', index_col=0, header=None)

    data.iloc[0, :] = data.iloc[0, :].map(
        lambda x: '+'.join(
            [part.split('_')[1].split('.')[0].split('-')[1] for part in x.split('+')]
        )
        if '+' in x
        else x.split('_')[1].split('.')[0].split('-')[1]
    )

    unique_values = data.iloc[1, :].unique()
    colors = plt.cm.tab20.colors
    color_map = {
        label: colors[i % len(colors)]
        for i, label in enumerate(data.iloc[0, :].unique())
    }

    fig, axes = plt.subplots(1, len(unique_values), figsize=(15, 5))

    sample_assignment = {}

    for ax, value in zip(axes, unique_values):
        subset = data.loc[:, data.iloc[1, :] == value]
        counts = subset.iloc[0, :].value_counts()
        logging.info(f'Counts for {value} in pool {pool_name}: {counts}')
        major_label = counts.idxmax()
        sample_assignment[value] = major_label
        if not any(counts > 0.75 * counts.sum()):
            logging.warning('Demultiplexing results in contaminated cell clusters.')
        counts.plot.pie(
            autopct='%1.1f%%',
            startangle=90,
            ax=ax,
            colors=[color_map[label] for label in counts.index],
        )
        ax.set_title(f'Abundancies for {value}')
        ax.set_ylabel('')

    plt.tight_layout()
    plt.savefig(
        Path(assignment_file).parent / f'pool_{pool_name}_demultiplexing_assignment.png'
    )
    return sample_assignment


def main(args):
    sample_assignments = {}
    if isinstance(args.demultiplexing_assignment, str):
        assignment_files = [args.demultiplexing_assignment]
    else:
        assignment_files = args.demultiplexing_assignment
    for assignment_file in assignment_files:
        pool_name = str(Path(assignment_file).stem).split('_')[1]
        sample_assignment = analyse_demultiplexing_assignment(
            assignment_file, pool_name
        )
        sample_assignments[f'({pool_name})'] = sample_assignment

    with open(args.output, 'w') as outfile:
        yaml.dump(sample_assignments, outfile)
    logging.info('Success.')


if __name__ == '__main__':
    args = parse_args()
    main(args)
