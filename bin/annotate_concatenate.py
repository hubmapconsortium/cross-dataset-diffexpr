#!/usr/bin/env python3

# Load multiple .h5ad files from different directories, probably the filtered cluster_marker_genes.h5ad files

# Annotate each file/cell/cluster with the data set/tissue of origin -- for labeling the overall data set, you can add a new key to the .uns mapping in the AnnData object, other annotations would probably go in .obs

# Concatenate all of the AnnData objects, ensuring that they have the same columns (genes, stored in AnnData.var) -- might need to expand each AnnData object if loading the filtered versions
import json
from argparse import ArgumentParser
from collections import defaultdict
from os import fspath, walk
from pathlib import Path
from typing import List, Dict, Sequence, Tuple, Optional
from cross_dataset_common import get_tissue_type, get_gene_dicts, get_cluster_df, hash_cell_id, precompute_dataset_percentages, get_pval_dfs, make_minimal_adata, upload_files_to_s3
from concurrent.futures import ThreadPoolExecutor

import anndata
from anndata import AnnData
import pandas as pd
import scipy.sparse
import scanpy as sc
import numpy as np
import mathplotlib.pyplot as plt

GENE_MAPPING_DIRECTORIES = [
    Path(__file__).parent.parent / 'data',
    Path('/opt/data'),
]

GENE_LENGTH_PATHS = [Path('/opt/data/gencode-v35-gene-lengths.json'), Path('/opt/data/salmon-index-v1.2-gene-lengths.json')]

@contextmanager
def new_plot():
    """
    When used in a `with` block, clears matplotlib internal state
    after plotting and saving things. Probably not necessary to be this
    thorough in clearing everything, but extra calls to `plt.clf()` and
    `plf.close()` don't *hurt*
    Intended usage:
        ```
        with new_plot():
            do_matplotlib_things()
            plt.savefig(path)
            # or
            fig.savefig(path)
        ```
    """
    plt.clf()
    try:
        yield
    finally:
        plt.clf()
        plt.close()

def get_annotation_metadata(filtered_files:List[Path]):
    for filtered_file in filtered_files:
        adata = anndata.read(filtered_file)
        if 'annotation_metadata' in adata.uns.keys() and adata.uns['annotation_metadata']['is_annotated']:
            return {'annotation_metadata':adata.uns['annotation_metadata']}
    return {'annotation_metadata':{'is_annotated':False}}

def get_inverted_gene_dict():
    inverted_dict = defaultdict(list)
    gene_mapping = read_gene_mapping()
    for ensembl, hugo in gene_mapping.items():
        inverted_dict[hugo].append(ensembl)
    return inverted_dict

def counts_to_rpkm(adata):
    cell_totals = adata.X.sum(axis=1).A.flatten()
    cell_total_recip = scipy.sparse.diags(1 / cell_totals)

    with open(GENE_LENGTH_PATHS[0]) as f:
        gene_lengths = json.load(f)

    with open(GENE_LENGTH_PATHS[1]) as f:
        gene_lengths.update(json.load(f))

    inverted_gene_mapping = get_inverted_gene_dict()

    gene_lengths_list = []
    for var in adata.var.index:
        ensembl_symbols = inverted_gene_mapping[var]
        ensembl_gene_lengths = [gene_lengths[symbol] for symbol in ensembl_symbols if symbol in gene_lengths]
        gene_lengths_list.append(sum(ensembl_gene_lengths))
    length_array = np.array(gene_lengths_list)
    length_recip = scipy.sparse.diags(1 / length_array)

    X = cell_total_recip @ adata.X @ length_recip
    return scipy.sparse.csc_matrix(X * 1e9)

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
    unfiltered_patterns = ['out.h5ad', 'expr.h5ad']
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

    if 'predicted.ASCT.celltype' in filtered_adata.obs.columns:
        filtered_adata.obs['cell_type'] = filtered_adata.obs['predicted.ASCT.celltype']
    else:
        filtered_adata.obs['cell_type'] = 'unknown'

    cells = list(filtered_adata.obs.index)
    unfiltered_subset = unfiltered_adata[cells,:].copy()
    unfiltered_subset.obs = filtered_adata.obs
    unfiltered_subset.obsm = filtered_adata.obsm

    cell_ids_list = ["-".join([data_set_dir, barcode]) for barcode in unfiltered_subset.obs['barcode']]
    unfiltered_subset.obs['cell_id'] = pd.Series(cell_ids_list, index=unfiltered_subset.obs.index)
    unfiltered_subset.obs.set_index("cell_id", drop=True, inplace=True)

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

    adata.layers["rpkm"] = counts_to_rpkm(adata)
    # This introduces duplicate gene names, use Pandas for aggregation
    # since anndata doesn't have that functionality
    adata.X = scipy.sparse.csr_matrix(adata.X)
    adata.var_names_make_unique()
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

def batch_correct_umap_cluster(adata):
    print(adata.obsm)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    #    sc.pp.filter_cells(adata, min_genes=200)
    #    sc.pp.filter_genes(adata, min_cells=3)

    adata.raw = adata

    adata.obsm["umap"] = adata.obsm["X_umap"]

    adata.obsm.pop("umap")
    adata.obsm.pop("X_umap")

    adata.obs["n_counts"] = adata.X.sum(axis=1)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.combat(adata, "dataset")
    sc.pp.scale(adata, max_value=10)

    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50)

    sc.tl.umap(adata)

    # leiden clustering
    sc.tl.leiden(adata)

    with new_plot():
        sc.pl.umap(adata, color="leiden", show=False)
        plt.savefig("umap_by_leiden_cluster.pdf", bbox_inches="tight")

    with new_plot():
        sc.pl.umap(adata, color="dataset", show=False)
        plt.savefig("umap_by_dataset.pdf", bbox_inches="tight")

    with new_plot():
        sc.pl.umap(adata, color="organ", show=False)
        plt.savefig("umap_by_organ.pdf", bbox_inches="tight")

    leiden_list = [f"leiden-UMAP-allrna-{adata.obs['leiden'][i]}" for i in adata.obs.index]

    adata.obs['leiden'] = pd.Series(leiden_list, index=adata.obs.index)

    # Write out as h5ad
    output_file = Path('bc_umap_cluster.h5ad')
    print('Saving output to', output_file.absolute())
    print(adata.layers)

#    adata.write_h5ad(output_file)

def annotate_single_gene(param_tuple):
    df = param_tuple[0]
    var_id = param_tuple[1]
    return_dict = {'gene_symbol':var_id}
    gene_df = df[df[var_id] > 0]
    return_dict['num_cells'] = len(gene_df.index)
    return_dict['num_datasets'] = len(gene_df['dataset'].unique())
    return return_dict

def annotate_genes(adata):
    df = adata.to_df()
    df['dataset'] = adata.obs['dataset']
    params_tuples = [(df[['dataset', var_id]], var_id) for var_id in df.columns if var_id != 'dataset']

    with ThreadPoolExecutor(max_workers=20) as e:
        dict_list = e.map(annotate_single_gene, params_tuples)

    return pd.DataFrame(dict_list)

def make_cell_df(adata):

    cell_df = adata.obs.copy()
    clusters_list = [",".join([cell_df["leiden"][i], cell_df["dataset_leiden"][i]]) for i in cell_df.index]
    cell_df["clusters"] = pd.Series(clusters_list, index=cell_df.index)

    umap_one_list = [coord[0] for coord in adata.obsm['X_umap']]
    umap_two_list = [coord[1] for coord in adata.obsm['X_umap']]

    cell_df['umap_1'] = pd.Series(umap_one_list, index=cell_df.index)
    cell_df['umap_2'] = pd.Series(umap_two_list, index=cell_df.index)

    cell_df = cell_df[['cell_id', 'barcode', 'dataset', 'organ', 'modality', 'clusters', 'umap_1', 'umap_2', 'cell_type']]

    return cell_df

def make_grouping_dfs(adata, old_cluster_df):
    #This can be threaded in the cdcommon lib
    #@TODO: Add a cell type df
    organ_df, cluster_df = get_pval_dfs(adata)

    cluster_df_list = cluster_df.to_dict(orient='records')
    cluster_df_list.extend(old_cluster_df.to_dict(orient='records'))

    cluster_df = pd.DataFrame(cluster_df_list)

    return cluster_df, organ_df

def main(token: str, directories: List[Path], access_key_id:str, secret_access_key:str, include_all_outputs: bool):

    token = None if token == "None" else token
    # Load files
    file_pairs = [find_file_pairs(directory) for directory in directories]
    filtered_files = [fp[0] for fp in file_pairs]
    annotation_metadata = get_annotation_metadata(filtered_files)
    annotated_files = [annotate_file(file_pair[0],file_pair[1], token) for file_pair in file_pairs]
    annotated_unfiltered_files = [file[0] for file in annotated_files]
    annotated_filtered_files = [file[1] for file in annotated_files]

    old_cluster_df = get_old_cluster_df(annotated_filtered_files)

    mapped_annotated_unfiltered_files = [map_gene_ids(file) for file in annotated_unfiltered_files]

    with ThreadPoolExecutor(max_workers=len(directories)) as e:
        percentage_dfs = e.map(precompute_dataset_percentages, mapped_annotated_unfiltered_files)

    percentage_df = pd.concat(percentage_dfs)

    concatenated_file = anndata.concat(mapped_annotated_unfiltered_files, fill_value=0, index_unique='_', join='inner')

    dataset_leiden_list = [f"leiden-UMAP-{concatenated_file.obs.at[i, 'dataset']}-{concatenated_file.obs.at[i, 'leiden']}" for i in concatenated_file.obs.index]
    concatenated_file.obs['dataset_leiden'] = pd.Series(dataset_leiden_list, index=concatenated_file.obs.index)
    concatenated_file.uns['annotation_metadata'] = annotation_metadata['annotation_metadata']

    bc_adata = batch_correct_umap_cluster(concatenated_file)

    cluster_df, organ_df = make_grouping_dfs(bc_adata, old_cluster_df)

    gene_df = annotate_genes(bc_adata)

    bc_adata.X = bc_adata.layers["rpkm"]
    cell_df = make_cell_df(bc_adata)
    minimal_adata = make_minimal_adata(bc_adata)
    minimal_adata.write('rna.h5ad')

    with pd.HDFStore('rna.hdf5') as store:
        store.put('cell', cell_df, format='t')
        store.put('organ', organ_df)
        store.put('cluster', cluster_df)
        store.put('gene', gene_df)

    with pd.HDFStore('rna_precompute.hdf5') as store:
        store.put('percentages', percentage_df)

    bc_adata.write_zarr('rna.zarr', chunks=[bc_adata.shape[0], 2])

    files_to_upload = [Path('rna.hdf5'), Path('rna_precompute.hdf5'), Path('rna.h5ad'), Path('rna.zarr')]
    upload_files_to_s3(files_to_upload, access_key_id, secret_access_key)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('nexus_token', type=str)
    p.add_argument('data_directories', type=Path, nargs='+')
    p.add_argument('access_key_id', type=str)
    p.add_argument('secret_access_key', type=str)
    p.add_argument("--enable-manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.nexus_token, args.data_directories, args.access_key_id, args.secret_access_key)
