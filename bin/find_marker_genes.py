#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_dfs, make_quant_df, create_minimal_dataset
import anndata
import pandas as pd
from hubmap_cell_id_gen_py import get_sequencing_cell_id


def main(h5ad_file: Path, old_cluster_file:Path):
    adata = anndata.read_h5ad(h5ad_file)
    cell_id_list = [get_sequencing_cell_id(adata.obs["dataset"][i], adata.obs["barcode"][i]) for i in adata.obs.index]
    adata.obs["cell_id"] = pd.Series(cell_id_list, index=adata.obs.index)

    #This can be threaded in the cdcommon lib
    organ_df, cluster_df = get_pval_dfs(adata)

    with pd.HDFStore(old_cluster_file) as store:
        old_cluster_df = store.get('cluster')

    cluster_df_list = cluster_df.to_dict(orient='records')
    cluster_df_list.extend(old_cluster_df.to_dict(orient='records'))

#    cluster_df = pd.concat([old_cluster_df, cluster_df])
    cluster_df = pd.DataFrame(cluster_df_list)

    cell_df = adata.obs.copy()
    clusters_list = [",".join([cell_df["leiden"][i], cell_df["dataset_leiden"][i]]) for i in cell_df.index]
    cell_df["clusters"] = pd.Series(clusters_list, index=cell_df.index)
    cell_df = cell_df[['cell_id', 'barcode', 'dataset', 'organ', 'modality', 'clusters']]

    adata.X = adata.raw.X

    quant_df = make_quant_df(adata)

    load_data_to_vms('rna', cell_df, quant_df, organ_df, cluster_df)

#    quant_df.to_csv('rna.csv')

#    with pd.HDFStore('rna.hdf5') as store:
#        store.put('cell', cell_df, format='t')
#        store.put('organ', organ_df)
#        store.put('cluster', cluster_df)

    create_minimal_dataset(cell_df, quant_df, organ_df, cluster_df, 'rna')

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    p.add_argument('old_cluster_file', type=Path)
    args = p.parse_args()

    main(args.bc_h5ad_file, args.old_cluster_file)
