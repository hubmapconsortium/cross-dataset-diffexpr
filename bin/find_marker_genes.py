#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_dfs, make_quant_csv
import anndata
import pandas as pd

def main(h5ad_file: Path, old_cluster_file:Path):
    adata = anndata.read_h5ad(h5ad_file)

    cell_df = adata.obs.copy()
    cell_df = cell_df[['cell_id', 'dataset', 'tissue_type', 'modality', 'leiden']].astype(str)

    make_quant_csv(adata, 'rna')

    organ_df, cluster_df = get_pval_dfs(adata)

    with pd.HDFStore(old_cluster_file) as store:
        old_cluster_df = store.get('cluster')

    cluster_df = pd.concat([old_cluster_df, cluster_df])

    with pd.HDFStore('rna.hdf5') as store:
        store.put('cell', cell_df, format='t')
        store.put('organ', organ_df)
        store.put('cluster', cluster_df)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    p.add_argument('old_cluster_file', type=Path)
    args = p.parse_args()

    main(args.bc_h5ad_file)
