#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_rows, add_quant_columns

import anndata
import pandas as pd

def get_quant_df(adata:anndata.AnnData)->pd.DataFrame:
    print(adata.X.shape)
    return pd.DataFrame(adata.X, columns=adata.var.index, index=adata.obs.index)

def main(h5ad_file: Path):
    adata = anndata.read_h5ad(h5ad_file)
    print(adata.X.shape)

    groupings = ['leiden', 'dataset', 'tissue_type']

    adata.obs['cell_id'] = adata.obs.index

    group_rows = get_rows(adata, groupings)

    cell_df = adata.obs.copy()

    quant_df = get_quant_df(adata)

    group_df = pd.DataFrame(group_rows, dtype=object)

    with pd.HDFStore('rna.hdf5') as store:
        store.put('cell', cell_df, format='t')
        store.put('group', group_df)
        store.put('quant', quant_df)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.bc_h5ad_file)
