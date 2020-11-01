#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_dfs, flatten_quant_df

import anndata
import pandas as pd

def get_quant_df(adata: anndata.AnnData) -> pd.DataFrame:
    print(adata.X.shape)
    return pd.DataFrame(adata.X.todense(), columns=adata.var.index, index=adata.obs.index)


def main(h5ad_file: Path):
    adata = anndata.read_h5ad(h5ad_file)
    adata.obs['cell_id'] = adata.obs.index

    pval_df = get_pval_dfs(adata)

    cell_df = adata.obs.copy()

    quant_df = get_quant_df(adata)

    long_df = flatten_quant_df(quant_df)
    long_df.to_csv('rna.csv')

    with pd.HDFStore('rna.hdf5') as store:
        store.put('cell', cell_df, format='t')
        store.put('p_values', pval_df)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.bc_h5ad_file)
