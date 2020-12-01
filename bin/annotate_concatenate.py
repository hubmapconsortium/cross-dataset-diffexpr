#!/usr/bin/env python3

# Load multiple .h5ad files from different directories, probably the filtered cluster_marker_genes.h5ad files

# Annotate each file/cell/cluster with the data set/tissue of origin -- for labeling the overall data set, you can add a new key to the .uns mapping in the AnnData object, other annotations would probably go in .obs

# Concatenate all of the AnnData objects, ensuring that they have the same columns (genes, stored in AnnData.var) -- might need to expand each AnnData object if loading the filtered versions
import json
from argparse import ArgumentParser
from functools import reduce
from os import fspath, walk
from pathlib import Path
from typing import List
from cross_dataset_common import get_tissue_type, get_gene_dicts, get_cluster_df, hash_cell_id

import anndata
import pandas as pd
import scanpy as sc
import numpy as np

def find_files(directory, patterns):
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            for pattern in patterns:
                if filepath.match(pattern):
                    return filepath

def find_file_pairs(directory):
    filtered_patterns = ['cluster_marker_genes.h5ad', 'secondary_analysis.h5ad']
    unfiltered_patterns = ['out.h5ad']
    filtered_file = find_files(directory, filtered_patterns)
    unfiltered_file = find_files(directory, unfiltered_patterns)
    return filtered_file, unfiltered_file

def annotate_file(filtered_file: Path, unfiltered_file: Path, token: str) -> anndata.AnnData:
    # Get the directory
    data_set_dir = fspath(unfiltered_file.parent.stem)
    # And the tissue type
    tissue_type = get_tissue_type(data_set_dir, token)

    filtered_adata = anndata.read_h5ad(filtered_file)
    unfiltered_adata = anndata.read_h5ad(unfiltered_file)

    cells = list(filtered_adata.obs.index)
    unfiltered_subset = unfiltered_adata[cells,:].copy()
    unfiltered_subset.obs = filtered_adata.obs

    unfiltered_subset.obs['barcode'] = unfiltered_subset.obs.index
    unfiltered_subset.obs['dataset'] = data_set_dir
    unfiltered_subset.obs['tissue_type'] = tissue_type
    unfiltered_subset.obs['modality'] = 'rna'
    semantic_cell_ids = data_set_dir + unfiltered_subset.obs.index
    unfiltered_subset.obs['cell_id'] = hash_cell_id(semantic_cell_ids)

    #    return adata
    return unfiltered_subset.copy()

def inner_join(adata_1: anndata.AnnData, adata_2: anndata.AnnData) -> anndata.AnnData:
    print(adata_1.X.shape)
    print(adata_2.X.shape)
    new_adata = adata_1.concatenate(adata_2, join='inner', fill_value=0)
    print(new_adata.X.shape)
    return new_adata


def main(token: str, directories: List[Path], ensembl_to_symbol_path=Path('/opt/ensembl_to_symbol.json'),
         symbol_to_ensembl_path=Path('/opt/symbol_to_ensembl.json')):
    # Load files
    file_pairs = [find_file_pairs(directory) for directory in directories]
    annotated_files = [annotate_file(file_pair[0],file_pair[1], token) for file_pair in file_pairs]
    for adata in annotated_files:
        sc.tl.rank_genes_groups(adata, 'leiden', method='t-test', rankby_abs=True)
    cluster_dfs = [get_cluster_df(adata) for adata in annotated_files]
    cluster_df = pd.concat(cluster_dfs)


    with pd.HDFStore('cluster.hdf5') as store:
        store.put('cluster', cluster_df)

    concatenated_file = reduce(inner_join, annotated_files)

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