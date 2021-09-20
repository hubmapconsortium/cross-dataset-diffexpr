#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import anndata
import scanpy as sc
import pandas as pd
from contextlib import contextmanager

import matplotlib.pyplot as plt

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
    print(adata.obsm)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.raw = adata

    adata.obsm["umap"] = adata.obsm["X_umap"]

    adata.obsm.pop("umap")
    adata.obsm.pop("X_umap")

    adata.obs["n_counts"] = adata.X.sum(axis=1)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.combat(adata, "dataset")
    sc.pp.scale(adata, max_value=10)

    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50)

    sc.tl.umap(adata)

    # leiden clustering
    sc.tl.leiden(adata)

    with new_plot():
        sc.pl.umap(adata, color="leiden", show=False)
        plt.savefig("umap_by_leiden_cluster.pdf", bbox_inches="tight")

    with new_plot():
        sc.pl.umap(adata, color="dataset", show=False)
        plt.savefig("umap_by_dataset.pdf", bbox_inches="tight")

    with new_plot():
        sc.pl.umap(adata, color="organ", show=False)
        plt.savefig("umap_by_organ.pdf", bbox_inches="tight")

    leiden_list = [f"leiden-UMAP-allrna-{adata.obs['leiden'][i]}" for i in adata.obs.index]

    adata.obs['leiden'] = pd.Series(leiden_list, index=adata.obs.index)

    # Write out as h5ad
    output_file = Path('bc_umap_cluster.h5ad')
    print('Saving output to', output_file.absolute())
    print(adata.layers)

    adata.write_h5ad(output_file)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('concatenated_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.concatenated_h5ad_file)
