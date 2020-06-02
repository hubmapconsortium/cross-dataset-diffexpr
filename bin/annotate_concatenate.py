#Load multiple .h5ad files from different directories, probably the filtered cluster_marker_genes.h5ad files

#Annotate each file/cell/cluster with the data set/tissue of origin -- for labeling the overall data set, you can add a new key to the .uns mapping in the AnnData object, other annotations would probably go in .obs

#Concatenate all of the AnnData objects, ensuring that they have the same columns (genes, stored in AnnData.var) -- might need to expand each AnnData object if loading the filtered versions

from argparse import ArgumentParser
from os import fspath, walk
from pathlib import Path
from subprocess import check_call
from typing import Dict, List, Tuple, Iterable

import anndata
import numpy as np
import pandas as pd

pattern = "*cluster_marker_genes.h5ad"
data_set_spreadsheet = "/home/sean/Documents/code/cross-dataset-diffexpr/spreadsheet.csv"

def find_h5ad_files(directory: Path) -> Iterable[Path]:
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            if filepath.match(pattern):
                yield filepath.relative_to(directory)

def get_tissue_type(file: Path) -> str:
    data_set_dir = fspath(file.parent.stem)

    with pd.read_csv(data_set_spreadsheet) as spreadsheet_df:
    #Open spreadsheet csv and select relevant parts
#        transcriptomics_df = spreadsheet_df[spreadsheet_df['Data Modality'] == 'Transcriptomics'].copy()

#        for i in range(len(transcriptomics_df.index)):
#            if data_set_dir in transcriptomics_df['localPath'][i]:
#                return transcriptomics_df['Organ/Tissue'][i]

        for i in range(len(spreadsheet_df.index)):
            if data_set_dir in spreadsheet_df['localPath'][i]:
                return spreadsheet_df['Organ/Tissue'][i]


def annotate_file(file: Path)->anndata.AnnData:
    #Get the directory
    data_set_dir = fspath(file.parent.stem)
    #And the tissue type
    tissue_type = get_tissue_type(file)

    #Add both to uns
    adata = anndata.read_h5ad(file)
    adata.uns['dataset'] = data_set_dir
    adata.uns['tissue_type'] = tissue_type

    return adata

def main(directory: Path):
    #Load files
    h5ad_files = find_h5ad_files(directory)
    #Annotate
    annotated_h5ad_files = [annotate_file(file) for file in h5ad_files]
    #Concatenate
    annotated_h5ad_files[0].concatenate(annotated_h5ad_files[1:])
    annotated_h5ad_files[0].write('concatenated_annotated_data.h5ad')

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('directory', type=Path)
    args = p.parse_args()

    main(args.directory)
