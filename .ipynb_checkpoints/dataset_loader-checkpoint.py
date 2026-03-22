import scanpy as sc

def load_demo_dataset(name):

    if name == "PBMC 3K":
        return sc.datasets.pbmc3k()

    elif name == "PBMC 68K Reduced":
        return sc.datasets.pbmc68k_reduced()

    elif name == "PBMC 3K Processed":
        adata = sc.datasets.pbmc3k()
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        return adata

    else:
        raise ValueError("Dataset not supported")