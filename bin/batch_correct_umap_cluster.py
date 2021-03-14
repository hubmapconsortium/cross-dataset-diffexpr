#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import anndata
import scanpy as sc
import numpy as np
import pandas as pd

def main(h5ad_file: Path):
    adata = anndata.read_h5ad(h5ad_file)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.raw = adata

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

    leiden_list = [f"leiden-UMAP-allrna-{adata.obs['leiden'][i]}" for i in adata.obs.index]

    adata.obs['leiden'] = pd.Series(leiden_list, index=adata.obs.index)

    # Write out as h5ad
    output_file = Path('bc_umap_cluster.h5ad')
    print('Saving output to', output_file.absolute())
    adata.write_h5ad(output_file)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('concatenated_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.concatenated_h5ad_file)
