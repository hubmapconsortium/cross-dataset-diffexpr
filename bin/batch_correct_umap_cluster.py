#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from memory-profiler import profile

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

@profile
def batch_correct_umap_cluster(adata):
    print(adata.obsm)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    #    sc.pp.filter_cells(adata, min_genes=200)
    #    sc.pp.filter_genes(adata, min_cells=3)

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

#    adata.write_h5ad(output_file)

def main(h5ad_file: Path):
    adata = anndata.read_h5ad(h5ad_file)

    dataset_tuples = [(uuid, len(adata.obs[adata.obs['dataset'] == uuid].index)) for uuid in adata.obs["dataset"].unique()]
    dataset_tuples.sort(key=lambda y: y[1])
    current_dataset = adata.obs["dataset"].iloc[0]
    adata.obs["num_in_dataset"] = pd.Series()
    current_count = 0
    for i in adata.obs.index:
        if adata.obs.at[i, "dataset"] != current_dataset:
                current_dataset = adata.obs.at[i, "dataset"]
                current_count = 0
        adata.obs.at[i, "num_in_dataset"] = current_count

    for num_datasets in [2, 4, 8, 16, 32, 64, 128]:
        for cells_per_dataset in [100, 1000, 10000]:
            print(f"{num_datasets} datasets with {cells_per_dataset} cells in each")
            uuids = [dataset_tuple[0] for dataset_tuple in dataset_tuples[0:num_datasets]]
            sub_adata = adata[sub_adata.obs["dataset"].isin(uuids)]
            sub_adata = sub_adata[sub_adata.obs["num_in_dataset"] <= cells_per_dataset]
            print(f"{len(sub_adata.obs.index)} total cells")
            batch_correct_umap_cluster(sub_adata)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('concatenated_h5ad_file', type=Path)
    p.add_argument("--enable-manhole", action="store_true")
    args = p.parse_args()

    if args.enable_manhole:
        import manhole
        manhole.install(activate_on="USR1")

    main(args.concatenated_h5ad_file)
