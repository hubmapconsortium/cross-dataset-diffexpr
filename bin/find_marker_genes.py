#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_and_organ_dfs

import anndata
import pandas as pd

def get_quant_df(adata: anndata.AnnData) -> pd.DataFrame:
    print(adata.X.shape)
    return pd.DataFrame(adata.X.todense(), columns=adata.var.index, index=adata.obs.index)

def make_long_df(quant_df: pd.DataFrame):
    dict_list = [{'cell_id': i, 'gene_id': column, 'value':quant_df.at[i, column], 'modality':'rna'} for i in quant_df.index for column in quant_df.columns]
    return pd.DataFrame(dict_list)

def main(h5ad_file: Path):
    adata = anndata.read_h5ad(h5ad_file)
    adata.obs['cell_id'] = adata.obs.index

    pval_df, organ_df = get_pval_and_organ_dfs(adata)

    cell_df = adata.obs.copy()

    quant_df = get_quant_df(adata)

    long_df = pd.DataFrame()
#    long_df = make_long_df(quant_df)
    long_df.to_csv('long_rna_quant.csv')

    with pd.HDFStore('rna.hdf5') as store:
        store.put('cell', cell_df, format='t')
        store.put('organ', organ_df)
        store.put('quant', quant_df)
        store.put('p_values', pval_df)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.bc_h5ad_file)
