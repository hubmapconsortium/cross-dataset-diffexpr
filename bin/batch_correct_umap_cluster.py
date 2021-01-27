#!/usr/bin/env python3
from argparse import ArgumentParser
from contextlib import contextmanager
from pathlib import Path

import anndata
# import matplotlib.pyplot as plt
import scanpy as sc


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
    adata = anndata.read_h5ad(h5ad_file)
    print(adata.X.shape)
    adata.var_names_make_unique()
#    with new_plot():
#        sc.pl.umap(adata, color='batch')
#        plt.savefig('umap_by_batch.pdf', bbox_inches='tight')

    # Batch correction with bbknn and dimension reduction
    sc.tl.pca(adata)
    print(adata.X.shape)
    sc.external.pp.bbknn(adata, batch_key='dataset')
    print(adata.X.shape)
#    with new_plot():
#        sc.pl.umap(adata, color='batch')
#        plt.savefig('umap_batch_corrected_by_batch.pdf', bbox_inches='tight')


    # add the total counts per cell as observations-annotation to adata
    adata.obs["n_counts"] = adata.X.sum(axis=1)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=5, min_disp=0)

    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50)

    sc.tl.umap(adata)

    adata.obs['dataset_leiden'] = adata.obs['leiden']
    print('leiden' in adata.obs.columns)
    adata.obs = adata.obs.drop('leiden', axis=1, inplace=False)
    print('leiden' in adata.obs.columns)

    print(adata.X.shape)
    # leiden clustering
    sc.tl.leiden(adata)
    print('leiden' in adata.obs.columns)
    print(adata.X.shape)

#    with new_plot():
#        sc.pl.umap(adata, color='leiden')
#        plt.savefig('umap_by_cluster.pdf', bbox_inches='tight')

    # Write out as h5ad
    output_file = Path('bc_umap_cluster.h5ad')
    print('Saving output to', output_file.absolute())
    adata.write_h5ad(output_file)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('concatenated_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.concatenated_h5ad_file)
