#!/usr/bin/env python3

# Load multiple .h5ad files from different directories, probably the filtered cluster_marker_genes.h5ad files

# Annotate each file/cell/cluster with the data set/tissue of origin -- for labeling the overall data set, you can add a new key to the .uns mapping in the AnnData object, other annotations would probably go in .obs

# Concatenate all of the AnnData objects, ensuring that they have the same columns (genes, stored in AnnData.var) -- might need to expand each AnnData object if loading the filtered versions
import json
from argparse import ArgumentParser
from os import fspath, walk
from pathlib import Path
from typing import List, Dict, Tuple
from cross_dataset_common import get_tissue_type, get_gene_dicts, get_cluster_df, hash_cell_id

import anndata
import pandas as pd
import scanpy as sc

CELL_GY_GENE_FILENAME = 'cell_by_gene.hdf5'
CELL_CLUSTER_FILENAME = 'umap_coords_clusters.csv'

GENE_MAPPING_DIRECTORIES = [
    Path(__file__).parent.parent / 'data',
    Path('/opt/data'),
]

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

def annotate_file(filtered_file: Path, unfiltered_file: Path, token: str) -> Tuple[anndata.AnnData, anndata.AnnData]:
    # Get the directory
    data_set_dir = fspath(unfiltered_file.parent.stem)
    # And the tissue type
    tissue_type = get_tissue_type(data_set_dir, token)

    filtered_adata = anndata.read_h5ad(filtered_file)
    unfiltered_adata = anndata.read_h5ad(unfiltered_file)

    filtered_adata.obs['barcode'] = filtered_adata.obs.index
    filtered_adata.obs['dataset'] = data_set_dir
    filtered_adata.obs['organ'] = tissue_type
    filtered_adata.obs['modality'] = 'rna'

    cells = list(filtered_adata.obs.index)
    unfiltered_subset = unfiltered_adata[cells,:].copy()
    unfiltered_subset.obs = filtered_adata.obs
    unfiltered_subset.obsm = filtered_adata.obsm
    print(unfiltered_subset.obsm.keys())

    return unfiltered_subset.copy(), filtered_adata

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

def map_gene_ids(adata):
    obsm = adata.obsm
    gene_mapping = read_gene_mapping()
    keep_vars = [gene in gene_mapping for gene in adata.var.index]
    adata = adata[:, keep_vars]
    temp_df = pd.DataFrame(adata.X.todense(), index=adata.obs.index, columns=adata.var.index)
    aggregated = temp_df.groupby(level=0, axis=1).sum()
    adata = anndata.AnnData(aggregated, obs=adata.obs)
    adata.var.index = [gene_mapping[var] for var in adata.var.index]
    adata.obsm = obsm
    # This introduces duplicate gene names, use Pandas for aggregation
    # since anndata doesn't have that functionality
    return adata


def get_old_cluster_df(annotated_filtered_files):
    for adata in annotated_filtered_files:
        dataset_leiden_list = [
            f"leiden-UMAP-{adata.obs['dataset'][i]}-{adata.obs['leiden'][i]}" for i in adata.obs.index]
        adata.obs["leiden"] = pd.Series(dataset_leiden_list, index=adata.obs.index)
        sc.tl.rank_genes_groups(adata, 'leiden', method='t-test', rankby_abs=True, n_genes=len(adata.var.index))
    cluster_dfs = [get_cluster_df(adata) for adata in annotated_filtered_files]
    cluster_df = pd.concat(cluster_dfs)

    with pd.HDFStore('cluster.hdf5') as store:
        store.put('cluster', cluster_df)


def main(token: str, directories: List[Path]):

    token = None if token == "None" else token
    # Load files
    file_pairs = [find_file_pairs(directory) for directory in directories]
    #This can be parallelized, though the benefits will likely be minimal
    annotated_files = [annotate_file(file_pair[0],file_pair[1], token) for file_pair in file_pairs]
    annotated_unfiltered_files = [file[0] for file in annotated_files]
    annotated_filtered_files = [file[1] for file in annotated_files]

    #This can be parallelized
    get_old_cluster_df(annotated_filtered_files)

    concatenated_file = annotated_unfiltered_files[0].concatenate(*annotated_unfiltered_files[1:], fill_value=0, index_unique='_')

    concatenated_file = map_gene_ids(concatenated_file)

    dataset_leiden_list = [f"leiden-UMAP-{concatenated_file.obs.at[i, 'dataset']}-{concatenated_file.obs.at[i, 'leiden']}" for i in concatenated_file.obs.index]
    concatenated_file.obs['dataset_leiden'] = pd.Series(dataset_leiden_list, index=concatenated_file.obs.index)

    concatenated_file.write('concatenated_annotated_data.h5ad')


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('nexus_token', type=str)
    p.add_argument('data_directories', type=Path, nargs='+')
    args = p.parse_args()

    main(args.nexus_token, args.data_directories)