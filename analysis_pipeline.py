import scanpy as sc


def run_pipeline(file_path):

    # Load dataset
    adata = sc.read(file_path)

    # Quality control
    sc.pp.filter_cells(adata, min_genes=200)

    # Normalize
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # Highly variable genes
    sc.pp.highly_variable_genes(adata)
    adata = adata[:, adata.var.highly_variable]

    # Scale data
    sc.pp.scale(adata)

    # PCA
    sc.tl.pca(adata)

    # Neighborhood graph
    sc.pp.neighbors(adata)

    # UMAP embedding
    sc.tl.umap(adata)

    # Clustering
    sc.tl.leiden(adata)

    return adata