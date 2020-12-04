#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_dfs, make_quant_df, create_minimal_dataset
import anndata
import pandas as pd

def percent_zero(df:pd.DataFrame)->float:
    zero_df = df[df['value'] == 0.0]
    return float(len(zero_df.index) / len(df.index))

def main(h5ad_file: Path, old_cluster_file:Path):
    adata = anndata.read_h5ad(h5ad_file)

    cell_df = adata.obs.copy()
    cell_df = cell_df[['cell_id', 'barcode', 'dataset', 'tissue_type', 'modality', 'dataset_leiden', 'leiden']].astype(str)

    quant_df = make_quant_df(adata)
    quant_df.to_csv('rna.csv')

    organ_df, cluster_df = get_pval_dfs(adata)

    with pd.HDFStore(old_cluster_file) as store:
        old_cluster_df = store.get('cluster')

    adj_organ_df, adj_cluster_df = get_pval_dfs(adata, True)

    print('Unadjusted organ df')
    print(percent_zero(organ_df))
    print('Adjusted organ df')
    print(percent_zero(adj_organ_df))
    print('Unadjusted cluster df')
    print(percent_zero(cluster_df))
    print('Adjusted cluster df')
    print(percent_zero(adj_cluster_df))

    cluster_df = pd.concat([old_cluster_df, cluster_df])

    with pd.HDFStore('rna.hdf5') as store:
        store.put('cell', cell_df, format='t')
        store.put('organ', organ_df)
        store.put('cluster', cluster_df)
        store.put('adj_organ', adj_organ_df)
        store.put('adj_cluster', adj_cluster_df)

    create_minimal_dataset(cell_df, quant_df, organ_df, cluster_df, 'rna')

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    p.add_argument('old_cluster_file', type=Path)
    args = p.parse_args()

    main(args.bc_h5ad_file, args.old_cluster_file)
