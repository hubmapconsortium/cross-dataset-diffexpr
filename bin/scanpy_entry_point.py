#!/usr/bin/env python3
from argparse import ArgumentParser
from contextlib import contextmanager
from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

import pickledb as pdb

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
    adatas = {}

    adata = anndata.read_h5ad(h5ad_file)
    adata.var_names_make_unique()

    with new_plot():
        sc.pl.umap(adata, color='batch')
        plt.savefig('umap_by_batch.pdf', bbox_inches='tight')

    #Batch correction with bbknn and dimension reduction
    sc.tl.pca(adata)
    sc.external.pp.bbknn(adata, batch_key='batch')

    with new_plot():
        sc.pl.umap(adata, color='batch')
        plt.savefig('umap_batch_corrected_by_batch.pdf', bbox_inches='tight')

    sc.tl.umap(adata)

    #leiden clustering
    sc.tl.leiden(adata)

    with new_plot():
        sc.pl.umap(adata, color='leiden')
        plt.savefig('umap_by_cluster.pdf', bbox_inches='tight')

    marker_gene_lists = {}

    #get number of genes
    num_genes = len(adata.var.index)

    for group_by in groupings:
    #for each thing we want to group by

        adatas[group_by] = adata.copy()
        sc.tl.rank_genes_groups(adatas[group_by], group_by, method='t-test', rankby_abs=True, n_genes=num_genes)

        #get the group_ids and then the gene_names and scores for each
        for group_id in adatas[group_by].obs[group_by].unique():
            cell_ids = adatas[group_by].obs[adatas[group_by].obs[group_by] == group_id].index
            for cell_id in cell_ids:
            #get cell ids for that group_id/grouping
                gene_names = adatas[group_by].uns['rank_genes_groups']['names'][group_id]
                scores = adatas[group_by].uns['rank_genes_groups']['scores'][group_id]
                names_and_scores = zip(gene_names, scores)
                for ns in names_and_scores:
                    if ns[0] not in marker_gene_lists.keys():
                        marker_gene_lists[ns[0]] = []
                    marker_gene_lists[ns[0]].append((cell_id, (group_by, group_id, str(ns[1]))))

        #Write out as h5ad
        output_file = Path('marker_genes_by_' + group_by + '_t_test.h5ad')
        print('Saving output to', output_file.absolute())
        adatas[group_by].write_h5ad(output_file)

        #And as pdf
        with new_plot():
            sc.pl.rank_genes_groups(adatas[group_by], n_genes=25, sharey=False)
            plt.savefig('marker_genes_by_' + group_by + '_t_test.pdf', bbox_inches='tight')

    db = pdb.load('marker_genes.db', False)

    for gene_name in marker_gene_lists.keys():
        db.set(gene_name, marker_gene_lists[gene_name])

    db.dump()

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('concatenated_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.concatenated_h5ad_file)
