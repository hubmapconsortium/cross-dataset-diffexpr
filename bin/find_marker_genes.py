#!/usr/bin/env python3
from argparse import ArgumentParser
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, List

import anndata
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import numpy as np

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

def get_gene_rows(gene_groupings:Dict[str, List[str]], marker_gene_groupings:Dict[str, List[str]]):

    gene_rows = []

    for gene in gene_groupings.keys():
        gene_id = gene
        gene_list = gene_groupings[gene]
        if gene in marker_gene_groupings.keys():
            marker_gene_list = marker_gene_groupings[gene]
        else:
            marker_gene_list = []
        gene_rows.append({'gene_id':gene_id, 'groups':gene_list, 'marker_groups':marker_gene_list})

    return gene_rows

def get_rows(adata:anndata.AnnData, groupings:List[str])->List[Dict]:

    cutoff = 0.9
    marker_cutoff = .001

    num_genes = len(adata.var_names)

    cell_df = adata.obs.copy()

    for group_by in groupings:
    #for each thing we want to group by

        sc.tl.rank_genes_groups(adata, group_by, method='t-test', rankby_abs=True, n_genes=num_genes)

        #get the group_ids and then the gene_names and scores for each
        for group_id in cell_df[group_by].unique():

            if type(group_id) == float and np.isnan(group_id):
                continue

            condition = group_by + "==" + str(group_id)

            gene_names = adata.uns['rank_genes_groups']['names'][group_id]
            pvals = adata.uns['rank_genes_groups']['pvals'][group_id]
            names_and_pvals = zip(gene_names, pvals)

            for n_p in names_and_pvals:

                if n_p[1] < cutoff:
                    if n_p[0] not in gene_groupings.keys():
                        gene_groupings[n_p[0]] = []
                    gene_groupings[n_p[0]].append(condition)

                if n_p[1] < marker_cutoff:
                    if n_p[0] not in marker_gene_groupings.keys():
                        marker_gene_groupings[n_p[0]] = []
                    marker_gene_groupings[n_p[0]].append(condition)

            genes = [n_p[0] for n_p in names_and_pvals if np[1] < cutoff]
            marker_genes = [n_p[0] for n_p in names_and_pvals if np[1] < marker_cutoff]

            group_rows.append({'condition':condition, 'genes':genes, 'marker_genes':marker_genes})

    gene_rows = get_gene_rows(gene_groupings, marker_gene_groupings)

    return group_rows, gene_rows

def main(h5ad_file: Path):

    adata = anndata.read_htad(h5ad_file)

    groupings = ['cluster', 'dataset', 'tissue_type']

    group_rows, gene_rows = get_rows(adata, groupings)

    cell_df = adata.obs.copy()
    gene_df = pd.DataFrame(gene_rows)
    group_df = pd.DataFrame(group_rows)

    cell_df.to_csv('rna.csv')
    group_df.to_csv('rna_group.csv')
    gene_df.to_csv('rna_gene.csv')

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.bc_h5ad_file)
