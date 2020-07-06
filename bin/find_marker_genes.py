#!/usr/bin/env python3
from argparse import ArgumentParser
from contextlib import contextmanager
from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import numpy as np
import sqlite3

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

def main(h5ad_file: Path):

    expression_cutoff = .5
    over_expression_cutoff = .001

    groupings = ['leiden', 'dataset', 'tissue_type']

    expressed_genes_dict = {}
    overexpressed_genes_dict = {}

    adata = anndata.read_h5ad(h5ad_file)
    adata.var_names_make_unique()

    #get number of genes
    num_genes = len(adata.var.index)

    df = adata.obs.copy()

    df['expressed_genes'] = pd.Series(dtype=str)
    df['overexpressed_genes'] = pd.Series(dtype=str)

    for group_by in groupings:
    #for each thing we want to group by

        sc.tl.rank_genes_groups(adata, group_by, method='t-test', rankby_abs=True, n_genes=num_genes)

        #get the group_ids and then the gene_names and scores for each
        for group_id in df[group_by].unique():
            df_select = df[df[group_by] == group_id]

            gene_names = adata.uns['rank_genes_groups']['names'][group_id]
            pvals = adata.uns['rank_genes_groups']['pvals'][group_id]
            names_and_pvals = zip(gene_names, pvals)

            expressed_genes = [np[0] for np in names_and_pvals if np[1] < expression_cutoff]
            overexpressed_genes = [np[0] for np in names_and_pvals if np[1] < over_expression_cutoff]

            for index in df_select.index:
                if index not in expressed_genes_dict.keys():
                    expressed_genes_dict[index] = []
                    overexpressed_genes_dict[index] = []
                expressed_genes_dict[index].extend(expressed_genes)
                overexpressed_genes_dict[index].extend(overexpressed_genes)

        #And as pdf
        with new_plot():
            sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
            plt.savefig('marker_genes_by_' + group_by + '_t_test.pdf', bbox_inches='tight')


    for index in expressed_genes_dict.keys():
        expressed_genes_set = set(expressed_genes_dict[index])
        overexpressed_genes_set = set(overexpressed_genes_dict[index])

        expressed_gene_string = ", ".join(expressed_genes_set)
        overexpressed_gene_string = ", ".join(overexpressed_genes_set)

        df.at[index, 'expressed_genes'] = expressed_gene_string
        df.at[index, 'overexpressed_genes'] = overexpressed_gene_string


    print('If you see this first, then conversion to categorical is happening at write-out')
    #Write out as sql db
    output_file = Path('rna.db')
    print('Saving output to', output_file.absolute())
    conn = sqlite3.connect(output_file)
    df.to_sql('rna', conn, if_exists='replace', index=True)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.bc_h5ad_file)
