#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from cross_dataset_common import get_pval_dfs, make_quant_df, create_minimal_dataset
import anndata
import pandas as pd
import requests
from hubmap_cell_id_gen_py import get_sequencing_cell_id
from concurrent.futures import ThreadPoolExecutor

SERVERS = ["https://cells.test.hubmapconsortium.org/api/"]
CHUNK_SIZE = 1000

def delete_data_from_servers(modality):
    request_dict = {"modality":modality}
    for server in SERVERS:
        request_url = server + "delete/"
        requests.post(request_url, request_dict)

def make_post_request(params_tuple):
    server = params_tuple[0]
    request_dict = params_tuple[1]
    request_url = server + "insert/"
    requests.post(request_url, request_dict)

def add_data_to_server(kwargs_list, model_name):
    for i in range(len(kwargs_list) // CHUNK_SIZE + 1):
        kwargs_subset = kwargs_list[i * CHUNK_SIZE:(i+1)*CHUNK_SIZE]
        request_dict = {"model_name":model_name, "kwargs_list":kwargs_subset}
        params_tuples = [(server, request_dict) for server in SERVERS]
        with ThreadPoolExecutor(max_workers=len(params_tuples)) as e:
            e.map(make_post_request, params_tuples)

def create_modality(modality):
    add_data_to_server([{"modality_name":modality}], "modality")

def create_datasets(datasets_list, modality):
    kwargs_list = [{"modality":modality, "uuid":dataset} for dataset in datasets_list]
    add_data_to_server(kwargs_list, "dataset")

def create_organs(cell_df):
    organs_list = list(cell_df["organ"].unique())
    kwargs_list = [{"grouping_name":organ} for organ in organs_list]
    add_data_to_server(kwargs_list, "organ")

def create_clusters(cell_df):
    clusters_set = {}
    unique_cluster_lists = [string.split(",") for string in cell_df["clusters"].unique()]
    for cluster_list in unique_cluster_lists:
        for cluster in cluster_list:
            clusters_set.add(cluster)

    cluster_splits = [cluster.split("-") for cluster in clusters_set]
    kwargs_list = [{"cluster_method":cs[0], "cluster_data":cs[1], "dataset":cs[2], "grouping_name":"-".join(cs)}
                   for cs in cluster_splits]
    add_data_to_server(kwargs_list, "cluster")

    return list(clusters_set)

def create_genes(quant_df):
    gene_symbols = list(quant_df["q_var_id"].unique())
    kwargs_list = [{"gene_symbol":gene_symbol} for gene_symbol in gene_symbols]
    add_data_to_server(kwargs_list, "gene")

def create_proteins(quant_df):
    protein_ids = list(quant_df["q_var_id"].unique())
    kwargs_list = [{"protein_id":protein_id} for protein_id in protein_ids]
    add_data_to_server(kwargs_list, "protein")

def create_cells(cell_df):
    cell_df = cell_df[["cell_id", "modality", "dataset", "organ"]]
    kwargs_list = cell_df.to_dict("records")
    add_data_to_server(kwargs_list, "cell")

def create_quants(quant_df, modality):
    model_name = modality + "quant"
    kwargs_list = quant_df.to_dict("records")
    add_data_to_server(kwargs_list, model_name)

def create_pvals(grouping_df):
    kwargs_list = grouping_df.to_dict("records")
    add_data_to_server(kwargs_list, "pvalue")

def set_up_relationships(cell_df, clusters):
    cell_clusters_dict = {}
    for cluster in clusters:
        cell_clusters_dict[cluster] = []
        for i in cell_df.index:
            if cluster in cell_df.at[i, "clusters"].split(","):
                cell_clusters_dict[cluster].append(cell_df.at[i, "cell_id"])

    for server in SERVERS:
        request_url = server + "setuprelationships/"
        requests.post(request_url, cell_clusters_dict)

def load_data_to_vms(modality, cell_df, quant_df, organ_df = None, cluster_df = None):
    delete_data_from_servers(modality)
    datasets_list = list(cell_df["dataset"].unique())
    create_datasets(datasets_list, modality)
    clusters = create_clusters(cell_df)
    create_organs(cell_df)
    if modality in ["rna", "atac"]:
        create_genes(quant_df)
    elif modality in ["codex"]:
        create_proteins(quant_df)

    #Do things up to here first because later things depend on them
    with ThreadPoolExecutor(max_workers=5) as e:
        e.submit(create_cells, cell_df)
        e.submit(create_quants, quant_df, modality)
        if modality in ["rna", "atac"]:
            e.submit(create_pvals, organ_df)
            e.submit(create_pvals,cluster_df)
        e.submit(set_up_relationships, cell_df, clusters)

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
