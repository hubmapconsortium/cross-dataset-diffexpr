#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_dfs, make_quant_df, create_minimal_dataset
import anndata
import pandas as pd
from hubmap_cell_id_gen_py import get_sequencing_cell_id
import scipy

import json
from typing import Dict
import numpy as np

GENE_MAPPING_DIRECTORIES = [
    Path(__file__).parent.parent / 'data',
    Path('/opt/data'),
]

GENE_LENGTH_PATHS = [Path('/opt/data/gencode-v35-gene-lengths.json'), Path('/opt/data/salmon-index-v1.2-gene-lengths.json')]

def get_inverted_gene_dict():
    inverted_dict = {}
    gene_mapping = read_gene_mapping()
    for key in gene_mapping:
        if gene_mapping[key] not in inverted_dict:
            inverted_dict[gene_mapping[key]] = []
        inverted_dict[gene_mapping[key]].append(key)
    return inverted_dict

def read_gene_mapping() -> Dict[str, str]:
    """
    Try to find the Ensembl to HUGO symbol mapping, with paths suitable
    for running this script inside and outside a Docker container.
    :return:
    """
    for directory in GENE_MAPPING_DIRECTORIES:
        mapping_file = directory / 'ensembl_to_symbol.json'
        if mapping_file.is_file():
            with open(mapping_file) as f:
                return json.load(f)
    message_pieces = ["Couldn't find Ensembl → HUGO mapping file. Tried:"]
    message_pieces.extend(f'\t{path}' for path in GENE_MAPPING_DIRECTORIES)
    raise ValueError('\n'.join(message_pieces))

def counts_to_rpkm(adata):
    cell_totals = adata.X.sum(axis=1)
    cell_total_recip = scipy.sparse.diags(1 / cell_totals)

    with open(GENE_LENGTH_PATHS[0]) as f:
        gene_lengths = json.load(f)

    with open(GENE_LENGTH_PATHS[1]) as f:
        gene_lengths.update(json.load(f))

    inverted_gene_mapping = get_inverted_gene_dict()

    gene_lengths_list = []
    for var in adata.var.index:
        if var in inverted_gene_mapping:
            ensembl_symbols = inverted_gene_mapping[var]
        else:
            var_prefix = var.split('-')[:-1]
            ensembl_symbols = inverted_gene_mapping[''.join(var_prefix)]
        ensembl_gene_lengths = [gene_lengths[symbol] for symbol in ensembl_symbols if symbol in gene_lengths]
        gene_lengths_list.append(sum(ensembl_gene_lengths))

    length_array = np.array(gene_lengths_list)
    length_recip = scipy.sparse.diags(1 / length_array)

    X = cell_total_recip @ adata.X @ length_recip
    return X * 1e9

def main(h5ad_file: Path, old_cluster_file:Path):
    adata = anndata.read_h5ad(h5ad_file)
    cell_id_list = [get_sequencing_cell_id(adata.obs["dataset"][i], adata.obs["barcode"][i]) for i in adata.obs.index]
    adata.obs["cell_id"] = pd.Series(cell_id_list, index=adata.obs.index)

    #This can be threaded in the cdcommon lib
#    organ_df, cluster_df = get_pval_dfs(adata)

#    old_cluster_df = pd.read_hdf(old_cluster_file, 'cluster')

#    cluster_df_list = cluster_df.to_dict(orient='records')
#    cluster_df_list.extend(old_cluster_df.to_dict(orient='records'))

#    cluster_df = pd.concat([old_cluster_df, cluster_df])
#    cluster_df = pd.DataFrame(cluster_df_list)

    cell_df = adata.obs.copy()
#    clusters_list = [",".join([cell_df["leiden"][i], cell_df["dataset_leiden"][i]]) for i in cell_df.index]
#    cell_df["clusters"] = pd.Series(clusters_list, index=cell_df.index)

    umap_one_list = [coord[0] for coord in adata.obsm['X_umap']]
    umap_two_list = [coord[1] for coord in adata.obsm['X_umap']]

#    cell_df['umap_1'] = pd.Series(umap_one_list, index=cell_df.index)
#    cell_df['umap_2'] = pd.Series(umap_two_list, index=cell_df.index)

#    cell_df = cell_df[['cell_id', 'barcode', 'dataset', 'organ', 'modality', 'clusters', 'umap_1', 'umap_2']]

    adata.X = counts_to_rpkm(adata)

    quant_df = make_quant_df(adata)

#    load_data_to_vms('rna', cell_df, quant_df, organ_df, cluster_df)

    quant_df.to_csv('rna.csv')

    with pd.HDFStore('rna.hdf5') as store:
        store.put('cell', cell_df, format='t')
#        store.put('organ', organ_df)
#        store.put('cluster', cluster_df)

#    create_minimal_dataset(cell_df, quant_df, organ_df, cluster_df, 'rna')
    cell_df.to_csv('mini_rna.csv')
    with pd.HDFStore('mini_rna.hdf5') as store:
        store.put('cell', cell_df, format='t')

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    p.add_argument('old_cluster_file', type=Path)
    args = p.parse_args()

    main(args.bc_h5ad_file, args.old_cluster_file)
