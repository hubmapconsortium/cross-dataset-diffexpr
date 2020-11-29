#!/usr/bin/env python3

# Load multiple .h5ad files from different directories, probably the filtered cluster_marker_genes.h5ad files

# Annotate each file/cell/cluster with the data set/tissue of origin -- for labeling the overall data set, you can add a new key to the .uns mapping in the AnnData object, other annotations would probably go in .obs

# Concatenate all of the AnnData objects, ensuring that they have the same columns (genes, stored in AnnData.var) -- might need to expand each AnnData object if loading the filtered versions
import json
from argparse import ArgumentParser
from functools import reduce
from os import fspath
from pathlib import Path
from typing import List
from cross_dataset_common import get_tissue_type, get_gene_dicts, get_cluster_df, hash_cell_id

import anndata
import pandas as pd
import scanpy as sc

pattern = "*out.h5ad"

def get_cluster_adata(h5ad_file):
    dataset = h5ad_file.parent.stem
    cluster_path = h5ad_file.parent / Path('cluster_marker_genes/cluster_marker_genes.h5ad')

    if cluster_path.exists():
        adata = anndata.read_h5ad(cluster_path)

    else:
        cluster_file = [file for file in h5ad_file.parent.iterdir() if file.stem in ['cluster_marker_genes', 'secondary_analysis']][0]
        adata = anndata.read_h5ad(cluster_file)

    adata.obs['dataset'] = dataset

    return adata

def annotate_file(file: Path, token: str) -> anndata.AnnData:
    # Get the directory
    data_set_dir = fspath(file.parent.stem)
    # And the tissue type
    tissue_type = get_tissue_type(data_set_dir, token)

    # Add both to uns
    adata = anndata.read_h5ad(file)
    adata.obs['barcode'] = adata.obs.index
    adata.obs['dataset'] = data_set_dir
    adata.obs['tissue_type'] = tissue_type
    adata.obs['modality'] = 'rna'
    semantic_cell_ids = data_set_dir + adata.obs.index
    adata.obs['cell_id'] = hash_cell_id(semantic_cell_ids)
    cluster_adata = get_cluster_adata(file)
    adata.obs['leiden'] = cluster_adata.obs['leiden']

    #    return adata
    return adata.copy()


def outer_join(adata_1: anndata.AnnData, adata_2: anndata.AnnData) -> anndata.AnnData:
    print(adata_1.X.shape)
    print(adata_2.X.shape)
    new_adata = adata_1.concatenate(adata_2, join='inner', fill_value=0)
    print(new_adata.X.shape)
    return new_adata


def main(token: str, directories: List[Path], ensembl_to_symbol_path=Path('/opt/ensembl_to_symbol.json'),
         symbol_to_ensembl_path=Path('/opt/symbol_to_ensembl.json')):
    # Load files
    h5ad_files = [directory / Path('out.h5ad') for directory in directories]
    annotated_files = [annotate_file(h5ad_file, token) for h5ad_file in h5ad_files]
    cluster_adatas = [get_cluster_adata(h5ad_file) for h5ad_file in h5ad_files]
    for cluster_adata in cluster_adatas:
        sc.tl.rank_genes_groups(cluster_adata, 'leiden', method='t-test', rankby_abs=True)
    cluster_dfs = [get_cluster_df(adata) for adata in cluster_adatas]
    cluster_df = pd.concat(cluster_dfs)


    with pd.HDFStore('cluster.hdf5') as store:
        store.put('cluster', cluster_df)

    concatenated_file = reduce(outer_join, annotated_files)

    symbol_to_ensembl_dict = {}
    ensembl_to_symbol_dict = {}

    var_columns = concatenated_file.var.index

    if ensembl_to_symbol_path.exists():
        with open(ensembl_to_symbol_path, 'r') as json_file:
            ensembl_to_symbol_dict = json.load(json_file)
        with open(symbol_to_ensembl_path, 'r') as json_file:
            symbol_to_ensembl_dict = json.load(json_file)

    else:
        ensembl_to_symbol_dict, symbol_to_ensembl_dict = get_gene_dicts(var_columns)
        with open('ensembl_to_symbol.json', 'w') as json_file:
            json.dump(ensembl_to_symbol_dict, json_file)
        with open('symbol_to_ensembl.json', 'w') as json_file:
            json.dump(symbol_to_ensembl_dict, json_file)

    keep_vars = [key for key in ensembl_to_symbol_dict.keys()]
    concatenated_file = concatenated_file[:, keep_vars]

    concatenated_file.var.index = [ensembl_to_symbol_dict[ensembl_id] for ensembl_id in concatenated_file.var.index]

    concatenated_file.write('concatenated_annotated_data.h5ad')


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('nexus_token', type=str)
    p.add_argument('data_directories', type=Path, nargs='+')
    args = p.parse_args()

    main(args.nexus_token, args.data_directories)
