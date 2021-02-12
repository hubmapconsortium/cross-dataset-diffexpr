#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import anndata
import scanpy as sc

def main(h5ad_file: Path):
    adata = anndata.read_h5ad(h5ad_file)
    adata.var_names_make_unique()

    sc.pp.combat(adata, key='dataset')

    adata.obs["n_counts"] = adata.X.sum(axis=1)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=5, min_disp=0)

    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50)

    sc.tl.umap(adata)

    adata.obs['dataset_leiden'] = adata.obs['leiden']
    adata.obs = adata.obs.drop('leiden', axis=1, inplace=False)

    # leiden clustering
    sc.tl.leiden(adata)

    # Write out as h5ad
    output_file = Path('bc_umap_cluster.h5ad')
    print('Saving output to', output_file.absolute())
    adata.write_h5ad(output_file)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('concatenated_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.concatenated_h5ad_file)
