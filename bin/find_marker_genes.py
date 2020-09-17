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

def add_quant_columns(df, adata):
    quant_df = adata.to_df()
    for column in quant_df.columns:
        if(type(quant_df[column].to_numpy()[0]) == float):
            df[column] = pd.Series(quant_df[column].to_numpy())
    return df.copy()

def main(h5ad_file: Path):
    adata = anndata.read_h5ad(h5ad_file)

    groupings = ['leiden', 'dataset', 'tissue_type']

    group_rows = get_rows(adata, groupings)

    cell_df = adata.obs.copy()
    cell_df = add_quant_columns(cell_df, adata)

    group_df = pd.DataFrame(group_rows, dtype=object)

    cell_df.to_csv('rna.csv')
    group_df.to_csv('rna_group.csv')


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.bc_h5ad_file)
