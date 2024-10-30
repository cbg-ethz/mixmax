import argparse
import numpy as np
import logging
from pathlib import Path
from itertools import islice
import h5py
from typing import Optional, Dict
import os
import shutil
import tempfile


import loompy


logging.basicConfig(
    format='{asctime} - {levelname} - {message}',
    style='{',
    datefmt='%Y-%m-%d %H:%M',
    level=logging.DEBUG,
)


def batched(iterable, chunk_size):
    iterator = iter(iterable)
    while chunk := tuple(islice(iterator, chunk_size)):
        yield chunk


class FancyLoomConnection(loompy.LoomConnection):
    def __init__(self, filename: str, mode: str = 'r+', *, validate: bool = True):
        super().__init__(filename, mode, validate=validate)

    def add_columns_batched(self, loom_view_manager, col_indices, batch_size=30):
        logging.info('subselecting cells in batched mode')
        batches = list(batched(col_indices, batch_size))
        for idx, batch in enumerate(batches):
            logging.info(f'Adding batch {idx+1} of {len(batches)}')
            ds_batch = loom_view_manager[:, np.array(batch)]
            self.add_columns(
                ds_batch.layers, col_attrs=ds_batch.ca, row_attrs=ds_batch.ra
            )


def fancy_loompy_connect(
    filename: str, mode: str = 'r+', *, validate: bool = True
) -> FancyLoomConnection:
    """
    Establish a connection to a .loom file.

    Args:
            filename:		Path to the Loom file to open
            mode:			Read/write mode, 'r+' (read/write) or 'r' (read-only), defaults to 'r+'
            validate:		Validate the file structure against the Loom file format specification
    Returns:
            A LoomConnection instance.

    Remarks:
            This function should typically be used as a context manager (i.e. inside a ``with``-block):

            .. highlight:: python
            .. code-block:: python

                    import loompy
                    with loompy.connect("mydata.loom") as ds:
                            print(ds.ca.keys())

            This ensures that the file will be closed automatically when the context block ends

            Note: if validation is requested, an exception is raised if validation fails.
    """
    return FancyLoomConnection(filename, mode, validate=validate)


def fancy_loompy_new(
    filename: str, *, file_attrs: Optional[Dict[str, str]] = None
) -> FancyLoomConnection:
    """
    Create an empty Loom file, and return it as a context manager.
    """
    if filename.startswith('~/'):
        filename = os.path.expanduser(filename)
    if file_attrs is None:
        file_attrs = {}

    # Create the file (empty).
    # Yes, this might cause an exception, which we prefer to send to the caller
    f = h5py.File(name=filename, mode='w')
    f.create_group('/attrs')  # v3.0.0
    f.create_group('/layers')
    f.create_group('/row_attrs')
    f.create_group('/col_attrs')
    f.create_group('/row_graphs')
    f.create_group('/col_graphs')
    f.flush()
    f.close()

    ds = fancy_loompy_connect(filename, validate=False)
    for vals in file_attrs:
        if file_attrs[vals] is None:
            ds.attrs[vals] = 'None'
        else:
            ds.attrs[vals] = file_attrs[vals]
    # store creation date
    ds.attrs['CreationDate'] = loompy.timestamp()
    ds.attrs['LOOM_SPEC_VERSION'] = loompy.loom_spec_version
    return ds


def sample_split(alpha, seed=None):
    if seed:
        rng = np.random.default_rng(seed=seed)
    else:
        rng = np.random.default_rng()
    splitting_ratio = rng.beta(alpha, alpha, size=1)[0]

    return splitting_ratio


def split_loom_file(splitting_ratio, loom_file, output_dir, seed=None):
    if seed:
        rng = np.random.default_rng(seed=seed)
    else:
        rng = np.random.default_rng()

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_loom_file = Path(temp_dir) / Path(loom_file).name
        shutil.copy2(loom_file, temp_loom_file)

        with loompy.connect(temp_loom_file) as ds:
            n_cells = ds.shape[1]
            split_vector = rng.choice(
                [0, 1], size=n_cells, p=[1 - splitting_ratio, splitting_ratio]
            )

            indices_first_sample = np.where(split_vector == 1)[0]
            indices_second_sample = np.where(split_vector == 0)[0]

            loom_file = Path(loom_file)
            ds_view = loompy.ViewManager(ds)
            output_1 = (output_dir / f'{loom_file.stem}_split1.loom').as_posix()

            with fancy_loompy_new(output_1) as ds1:
                logging.info(f'Connection to {output_1} established')
                ds1.add_columns_batched(ds_view, indices_first_sample, batch_size=30)
                number_of_cells1 = ds1.shape[1]
            logging.info('Writing file 1 was successful!')

            output_2 = (output_dir / f'{loom_file.stem}_split2.loom').as_posix()
            with fancy_loompy_new(output_2) as ds2:
                logging.info(f'Connection to {output_2} established')
                ds2.add_columns_batched(ds_view, indices_second_sample, batch_size=30)
                number_of_cells2 = ds2.shape[1]
            logging.info('Writing file 2 was successful!')
            if number_of_cells1 + number_of_cells2 == n_cells:
                logging.info('Splitting was successful!')
            else:
                logging.error(
                    'The number of cells in the split files does not match the number of cells in the original file'
                )
                raise ValueError(
                    f'Split 1 has {number_of_cells1} cells and split 2 has {number_of_cells2} cells, but the original file has {n_cells} cells'
                )
    return None


def parse_args():
    parser = argparse.ArgumentParser(
        description='Simulate multiplexed single-cell (or multi-sample) mutation data'
    )
    parser.add_argument(
        '--input', type=str, required=True, help='Input list of loom files'
    )
    parser.add_argument(
        '--alpha', type=float, required=True, help='Dirichlet distribution parameter'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='One of the output files. The other output file is saves to the same directory',
    )

    logging.info('Parsing arguments')
    args = parser.parse_args()
    return args


def main(args):
    logging.info(f'Splitting file {args.input}')
    loom_file = args.input
    splitting_ratio = sample_split(alpha=args.alpha)
    output_dir = Path(args.output).parent
    split_loom_file(
        splitting_ratio=splitting_ratio, loom_file=loom_file, output_dir=output_dir
    )
    logging.info(f'Done splitting file {args.input}')


if __name__ == '__main__':
    args = parse_args()
    main(args)
