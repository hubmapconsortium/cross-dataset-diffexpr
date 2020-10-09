#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_rows, add_quant_columns

import anndata
import pandas as pd

def get_quant_df(adata:anndata.AnnData)->pd.DataFrame:
    return pd.DataFrame(adata.X, columns=adata.var.index, index=adata.obs.index)

def main(h5ad_file: Path):
    adata = anndata.read_h5ad(h5ad_file)

    groupings = ['leiden', 'dataset', 'tissue_type']

    adata.obs['cell_id'] = adata.obs.index

    group_rows = get_rows(adata, groupings)

    cell_df = adata.obs.copy()

    quant_df = get_quant_df(adata)

    group_df = pd.DataFrame(group_rows, dtype=object)

    with pd.HDFStore('rna.hdf5') as store:
        store['cell'] = cell_df
        store['quant'] = quant_df
        store['group'] = group_df


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.bc_h5ad_file)
