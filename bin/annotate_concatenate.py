#!/usr/bin/env python3

# Load multiple .h5ad files from different directories, probably the filtered cluster_marker_genes.h5ad files

# Annotate each file/cell/cluster with the data set/tissue of origin -- for labeling the overall data set, you can add a new key to the .uns mapping in the AnnData object, other annotations would probably go in .obs

# Concatenate all of the AnnData objects, ensuring that they have the same columns (genes, stored in AnnData.var) -- might need to expand each AnnData object if loading the filtered versions
import json
from argparse import ArgumentParser
from functools import reduce
from os import fspath, walk
from pathlib import Path
from typing import List, Iterable, Dict

import anndata
import requests
import yaml

pattern = "*out.h5ad"


def get_gene_response(ensembl_ids: List[str]):
    request_url = 'https://mygene.info/v3/gene?fields=symbol'

    chunk_size = 1000
    chunks = (len(ensembl_ids) // chunk_size) + 1

    base_list = []

    for i in range(chunks):
        if i < chunks - 1:
            ensembl_slice = ensembl_ids[i * chunk_size: (i + 1) * chunk_size]
        else:
            ensembl_slice = ensembl_ids[i * chunk_size:]
        request_body = {'ids': ', '.join(ensembl_slice)}
        base_list.extend(requests.post(request_url, request_body).json())

    return base_list


def get_gene_dicts(ensembl_ids: List[str]) -> (Dict, Dict):
    #    temp_forwards_dict = {ensembl_id:ensembl_id.split('.')[0] for ensembl_id in ensembl_ids}
    temp_backwards_dict = {ensembl_id.split('.')[0]: ensembl_id for ensembl_id in ensembl_ids}
    ensembl_ids = [ensembl_id.split('.')[0] for ensembl_id in ensembl_ids]

    json_response = get_gene_response(ensembl_ids)

    forwards_dict = {temp_backwards_dict[item['query']]: item['symbol'] for item in json_response if
                     'symbol' in item.keys()}
    backwards_dict = {item['symbol']: temp_backwards_dict[item['query']] for item in json_response if
                      'symbol' in item.keys()}

    return forwards_dict, backwards_dict


def find_h5ad_files(directory: Path) -> Iterable[Path]:
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            if filepath.match(pattern):
                yield filepath


def get_tissue_type(dataset: str, token: str) -> str:
    organ_dict = yaml.load(open('/opt/organ_types.yaml'), Loader=yaml.BaseLoader)

    dataset_query_dict = {
        "query": {
            "bool": {
                "must": [],
                "filter": [
                    {
                        "match_all": {}
                    },
                    {
                        "exists": {
                            "field": "files.rel_path"
                        }
                    },
                    {
                        "match_phrase": {
                            "uuid": {
                                "query": dataset
                            },
                        }

                    }
                ],
                "should": [],
                "must_not": [
                    {
                        "match_phrase": {
                            "status": {
                                "query": "Error"
                            }
                        }
                    }
                ]
            }
        }
    }

    dataset_response = requests.post(
        'https://search.api.hubmapconsortium.org/search',
        json=dataset_query_dict,
        headers={'Authorization': 'Bearer ' + token})
    hits = dataset_response.json()['hits']['hits']

    for hit in hits:
        for ancestor in hit['_source']['ancestors']:
            if 'organ' in ancestor.keys():
                return organ_dict[ancestor['organ']]['description']


def annotate_file(file: Path, token: str) -> anndata.AnnData:
    # Get the directory
    data_set_dir = fspath(file.parent.stem)
    # And the tissue type
    tissue_type = get_tissue_type(data_set_dir, token)

    # Add both to uns
    adata = anndata.read_h5ad(file)
    adata.obs['dataset'] = data_set_dir
    adata.obs['tissue_type'] = tissue_type
    adata.obs['modality'] = 'rna'

    #    return adata
    return adata.copy()


def outer_join(adata_1: anndata.AnnData, adata_2: anndata.AnnData) -> anndata.AnnData:
    return adata_1.concatenate(adata_2, join='outer', fill_value=0)


def main(token: str, directories: List[Path], ensembl_to_symbol_path=Path('/opt/ensembl_to_symbol.json'),
         symbol_to_ensembl_path=Path('/opt/symbol_to_ensembl.json')):
    # Load files
    h5ad_files = [directory / Path('out.h5ad') for directory in directories]
    annotated_files = [annotate_file(h5ad_file, token) for h5ad_file in h5ad_files]
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
