#!/usr/bin/env python3

# Load multiple .h5ad files from different directories, probably the filtered cluster_marker_genes.h5ad files

# Annotate each file/cell/cluster with the data set/tissue of origin -- for labeling the overall data set, you can add a new key to the .uns mapping in the AnnData object, other annotations would probably go in .obs

# Concatenate all of the AnnData objects, ensuring that they have the same columns (genes, stored in AnnData.var) -- might need to expand each AnnData object if loading the filtered versions

from argparse import ArgumentParser
from functools import reduce
from os import fspath, walk
from pathlib import Path
from typing import List, Iterable

import anndata
import requests
import yaml

pattern = "*out.h5ad"


def ensembl_to_symbol(ensembl_id: str) -> str:
    ensembl_id = ensembl_id.split('.')[0]
    request_url = 'https://mygene.info/v3/gene/' + ensembl_id + '?fields=symbol&dotfield=True'
    r = requests.get(request_url)
    return r.json()['symbol']


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
        'https://search-api.dev.hubmapconsortium.org/search',
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

    new_var_index = [ensembl_to_symbol(gene) for gene in adata.var.index]
    adata.var.index = new_var_index

    #    return adata
    return adata.copy()


def outer_join(adata_1: anndata.AnnData, adata_2: anndata.AnnData) -> anndata.AnnData:
    return adata_1.concatenate(adata_2, join='outer', fill_value=0)


def main(token: str, directories: List[Path]):
    # Load files
    h5ad_files = [directory / Path('out.h5ad') for directory in directories]
    annotated_files = [annotate_file(h5ad_file, token) for h5ad_file in h5ad_files]
    concatenated_file = reduce(outer_join, annotated_files)

    concatenated_file.write('concatenated_annotated_data.h5ad')


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('nexus_token', type=str)
    p.add_argument('data_directories', type=Path, nargs='+')
    args = p.parse_args()

    main(args.nexus_token, args.data_directories)
