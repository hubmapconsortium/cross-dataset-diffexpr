#!/usr/bin/env python3
from argparse import ArgumentParser
from contextlib import contextmanager
from pathlib import Path
from cross_dataset_common import get_rows

import anndata
import matplotlib.pyplot as plt
import pandas as pd



@contextmanager
def new_plot():
    """
    When used in a `with` block, clears matplotlib internal state
    after plotting and saving things. Probably not necessary to be this
    thorough in clearing everything, but extra calls to `plt.clf()` and
    `plf.close()` don't *hurt*

    Intended usage:
        ```
        with new_plot():
            do_matplotlib_things()

            plt.savefig(path)
            # or
            fig.savefig(path)
        ```
    """
    plt.clf()
    try:
        yield
    finally:
        plt.clf()
        plt.close()

def main(h5ad_file: Path):
    adata = anndata.read_h5ad(h5ad_file)

    groupings = ['leiden', 'dataset', 'tissue_type']

    group_rows = get_rows(adata, groupings)

    cell_df = adata.obs.copy()

    group_df = pd.DataFrame(group_rows)

    cell_df.to_csv('rna.csv')
    group_df.to_csv('rna_group.csv')


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.bc_h5ad_file)
