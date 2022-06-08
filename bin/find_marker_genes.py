#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_dfs, create_minimal_dataset, make_minimal_adata, upload_files_to_s3
import anndata
import pandas as pd
from hubmap_cell_id_gen_py import get_sequencing_cell_id
from concurrent.futures import ThreadPoolExecutor

GENE_MAPPING_DIRECTORIES = [
    Path(__file__).parent.parent / 'data',
    Path('/opt/data'),
]

GENE_LENGTH_PATHS = [Path('/opt/data/gencode-v35-gene-lengths.json'), Path('/opt/data/salmon-index-v1.2-gene-lengths.json')]


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

def make_grouping_dfs(adata, old_cluster_file):
    #This can be threaded in the cdcommon lib
    #@TODO: Add a cell type df
    organ_df, cluster_df = get_pval_dfs(adata)

    old_cluster_df = pd.read_hdf(old_cluster_file, 'cluster')

    cluster_df_list = cluster_df.to_dict(orient='records')
    cluster_df_list.extend(old_cluster_df.to_dict(orient='records'))

    cluster_df = pd.DataFrame(cluster_df_list)

    return cluster_df, organ_df

def main(h5ad_file: Path, old_cluster_file:Path, access_key_id:str, secret_access_key:str):
    adata = anndata.read_h5ad(h5ad_file)
    cell_id_list = [get_sequencing_cell_id(adata.obs["dataset"][i], adata.obs["barcode"][i]) for i in adata.obs.index]
    adata.obs["cell_id"] = pd.Series(cell_id_list, index=adata.obs.index)

    adata.X = adata.layers["rpkm"]
    cell_df = make_cell_df(adata)

    cluster_df, organ_df = make_grouping_dfs(adata, old_cluster_file)

    gene_df = annotate_genes(adata)

    minimal_adata = make_minimal_adata(adata)
    minimal_adata.write('rna.h5ad')

    with pd.HDFStore('rna.hdf5') as store:
        store.put('cell', cell_df, format='t')
        store.put('organ', organ_df)
        store.put('cluster', cluster_df)
        store.put('gene', gene_df)

    files_to_upload = [Path('rna.hdf5'), Path('rna_precompute.hdf5'), Path('rna.h5ad')]
    upload_files_to_s3(files_to_upload, access_key_id, secret_access_key)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    p.add_argument('old_cluster_file', type=Path)
    p.add_argument('access_key_id', type=str)
    p.add_argument('secret_access_key', type=str)
    p.add_argument("--enable-manhole", action="store_true")
    args = p.parse_args()

    if args.enable_manhole:
        import manhole
        manhole.install(activate_on="USR1")

    main(args.concatenated_annotated_file, args.old_cluster_file, args.access_key_id, args.secret_access_key)
