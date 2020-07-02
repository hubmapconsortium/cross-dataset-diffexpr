#!/usr/bin/env python3
from argparse import ArgumentParser
from contextlib import contextmanager
from pathlib import Path

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

def main(h5ad_file: Path):

    groupings = ['leiden', 'dataset', 'tissue_type']

    adata = anndata.read_h5ad(h5ad_file)
    adata.var_names_make_unique()

    #get number of genes
    num_genes = len(adata.var.index)

    for group_by in groupings:
    #for each thing we want to group by

        sc.tl.rank_genes_groups(adata, group_by, method='t-test', rankby_abs=True, n_genes=num_genes)

        #get the group_ids and then the gene_names and scores for each
        for group_id in adata.obs[group_by].unique():
            df_select = adata.obs[adata.obs[group_by] == group_id]

            gene_names = adata.uns['rank_genes_groups']['names'][group_id]
            scores = adata.uns['rank_genes_groups']['pvals'][group_id]
            names_and_scores = zip(gene_names, scores)

            for ns in names_and_scores:
                column_name = ns[0] + '_' + group_by
                df_select[column_name] = ns[1]


        #Write out as h5ad
        output_file = Path('marker_genes_by_' + group_by + '_t_test.h5ad')
        print('Saving output to', output_file.absolute())
        adatas[group_by].write_h5ad(output_file)

        #And as pdf
        with new_plot():
            sc.pl.rank_genes_groups(adatas[group_by], n_genes=25, sharey=False)
            plt.savefig('marker_genes_by_' + group_by + '_t_test.pdf', bbox_inches='tight')


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bc_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.bc_h5ad_file)
