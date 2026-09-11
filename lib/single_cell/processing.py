import pandas as pd
import numpy as np
import scanpy as sc


def scanpy_workflow(
    adata,
    normalize_kw=None,
    hvg_kw=None,
    filter_hvg=None,
    regress_kw=None,
    scale_kw=None,
    pca_kw=None,
    nn_kw=None,
    umap_kw=None,
    resolutions=None,
    skip_normalize=False
):
    """Filters adata object by cell and gene level filters
    Parameters
    ----------
    adata : anndata, the adata object to filter
    Returns
    -------
    """
    # add defaults
    if normalize_kw is None:
        normalize_kw = {"target_sum": 1e4}

    if hvg_kw is None:
        hvg_kw = {"flavor": "seurat_v3", "layer": "counts", "n_top_genes": 2000}

    # seurat_v3 requires raw counts 
    if hvg_kw.get("flavor") == "seurat_v3" and hvg_kw.get("layer") != "counts":
        raise ValueError("flavor seurat_v3 requires a layer == 'counts'")

    if scale_kw is None:
        scale_kw = {"max_value": 10}

    if pca_kw is None:
        pca_kw = {'use_highly_variable': None}

    if nn_kw is None:
        nn_kw = {"n_neighbors": 20, "n_pcs": 20, "metric": "euclidean"} # this is the Seurat default

    if umap_kw is None:
        umap_kw = {"min_dist": 0.3}

    if resolutions is None:
        resolutions = [0.1, 0.2, 0.5, 0.9]

    if not skip_normalize:
        # begin by storing the raw form of the data
        adata.layers["counts"] = adata.X.copy()
        
        # normalize the data
        print("normalizing...")
        sc.pp.normalize_total(adata, **normalize_kw)
        adata.layers["normalized"] = adata.X.copy()

        sc.pp.log1p(adata)
        adata.layers["log_normalized"] = adata.X.copy()

    # Save normalized data as the raw attribute
    adata.raw = adata

    # identify highly variable genes
    sc.pp.highly_variable_genes(adata, **hvg_kw)

    # scale the data
    print("scaling...")
    sc.pp.scale(adata, **scale_kw)

    # store the scaled layer
    adata.layers["scaled"] = adata.X.copy()

    # run PCA
    print("running PCA...")
    sc.tl.pca(adata, **pca_kw)

    # run find neighbors first, match seurat's defaults
    print("finding neighbors...")
    sc.pp.neighbors(adata, **nn_kw)

    # compute umap embedding
    print("calculating UMAP...")
    sc.tl.umap(adata, **umap_kw)

    # compute the cluster assignments using leiden
    print("calculating clustering...")
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")

    return adata
