import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import streamlit as st
import scanpy as sc
import pandas as pd
import numpy as np
import plotly.express as px
import matplotlib.pyplot as plt

from dataset_loader import load_demo_dataset

st.set_page_config(layout="wide")

# ---------------- TITLE ----------------

st.title("Single Cell RNA-seq Immune Cell Analysis Dashboard")
st.markdown("###  Interactive scRNA-seq Analysis Platform")

# ---------------- SIDEBAR ----------------

st.sidebar.header("Dataset Loader")

dataset_option = st.sidebar.selectbox(
    "Select Dataset",
    [
        "PBMC 3K",
        "PBMC 68K Reduced",
        "PBMC 3K Processed",
        "Upload Your Own Dataset",
        "Load from URL"
    ]
)

# ---------------- DATA LOADING ----------------

if dataset_option == "Upload Your Own Dataset":

    uploaded_file = st.sidebar.file_uploader("Upload .h5ad file", type=["h5ad"])

    if uploaded_file is not None:
        try:
            adata = sc.read(uploaded_file)
            st.sidebar.success("Dataset loaded successfully ✅")
        except:
            st.sidebar.error("Error loading dataset")
            st.stop()
    else:
        st.stop()

elif dataset_option == "Load from URL":

    url = st.sidebar.text_input("Enter .h5ad URL")

    if url:
        try:
            with st.spinner("Downloading dataset..."):
                adata = sc.read(url)
            st.sidebar.success("Dataset loaded from URL ✅")
        except:
            st.sidebar.error("Failed to load dataset")
            st.stop()
    else:
        st.stop()

else:
    adata = load_demo_dataset(dataset_option)

# ---------------- PIPELINE ----------------

@st.cache_resource
def run_pipeline(_adata):

    sc.pp.normalize_total(_adata)
    sc.pp.log1p(_adata)
    sc.pp.highly_variable_genes(_adata)
    sc.pp.pca(_adata)
    sc.pp.neighbors(_adata)
    sc.tl.leiden(_adata)
    sc.tl.umap(_adata)

    return _adata

adata = run_pipeline(adata)

st.write("Dataset Shape:", adata.shape)

# ---------------- TABS ----------------

tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
    "Cell Clusters",
    "UMAP Explorer",
    "Gene Expression",
    "Marker Genes",
    "Cluster Comparison",
    "Cell Annotation",
    "Trajectory"
])

# ---------------- TAB 1 ----------------

with tab1:
    st.subheader("Cell Clusters")
    st.bar_chart(adata.obs["leiden"].value_counts())

# ---------------- TAB 2 ----------------

with tab2:
    st.subheader("UMAP Visualization")

    umap_df = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"])
    umap_df["cluster"] = adata.obs["leiden"].values

    fig = px.scatter(umap_df, x="UMAP1", y="UMAP2", color="cluster")
    st.plotly_chart(fig, use_container_width=True)

# ---------------- TAB 3 ----------------

with tab3:
    st.subheader("Gene Expression Explorer")

    gene = st.selectbox("Select Gene", adata.var_names[:1000])

    expr = adata[:, gene].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    expr = expr.flatten()

    umap_df = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"])
    umap_df["expression"] = expr

    fig = px.scatter(umap_df, x="UMAP1", y="UMAP2", color="expression")
    st.plotly_chart(fig, use_container_width=True)

# ---------------- TAB 4 ----------------

with tab4:
    st.subheader("Marker Gene Discovery")

    top_genes = st.slider("Marker genes (table/heatmap)", 5, 30, 10)
    dotplot_genes = st.slider("Dotplot genes", 3, 10, 5)

    sc.tl.rank_genes_groups(adata, "leiden", method="t-test")

    marker_df = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(top_genes)

    st.dataframe(marker_df)

    # DOWNLOAD MARKER GENES
    st.download_button(
        "Download Marker Genes CSV",
        marker_df.to_csv(index=False),
        "marker_genes.csv"
    )

    st.subheader("Heatmap")
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=top_genes, show=False)
    st.pyplot(plt.gcf())
    plt.clf()

    st.subheader("Dotplot")
    sc.pl.rank_genes_groups_dotplot(
        adata,
        n_genes=dotplot_genes,
        figsize=(12, 6),
        show=False
    )
    st.pyplot(plt.gcf())
    plt.clf()

# ---------------- TAB 5 ----------------

with tab5:

    st.subheader("Volcano Plot")

    clusters = sorted(adata.obs["leiden"].unique())

    cluster_a = st.selectbox("Cluster A", clusters)
    cluster_b = st.selectbox("Cluster B", clusters)

    label_count = st.slider("Label genes", 3, 15, 5)
    pval_threshold = st.slider("p-value cutoff", 0.001, 0.1, 0.05)
    logfc_threshold = st.slider("log2FC cutoff", 0.5, 2.0, 1.0)

    if cluster_a != cluster_b:

        sc.tl.rank_genes_groups(
            adata,
            groupby="leiden",
            groups=[cluster_a],
            reference=cluster_b,
            method="wilcoxon"
        )

        de_df = sc.get.rank_genes_groups_df(adata, group=cluster_a)

        volcano_df = pd.DataFrame({
            "gene": de_df["names"],
            "log2FC": de_df["logfoldchanges"],
            "pval": de_df["pvals"]
        })

        volcano_df["neglog10p"] = -np.log10(volcano_df["pval"])

        volcano_df["significant"] = (
            (volcano_df["pval"] < pval_threshold) &
            (abs(volcano_df["log2FC"]) > logfc_threshold)
        )

        fig = px.scatter(
            volcano_df,
            x="log2FC",
            y="neglog10p",
            color="significant",
            hover_name="gene"
        )

        # Threshold lines
        fig.add_vline(x=logfc_threshold, line_dash="dash")
        fig.add_vline(x=-logfc_threshold, line_dash="dash")
        fig.add_hline(y=-np.log10(pval_threshold), line_dash="dash")

        # Labels
        top_label_genes = volcano_df[volcano_df["significant"]] \
            .sort_values("pval") \
            .head(label_count)

        for _, row in top_label_genes.iterrows():
            fig.add_annotation(
                x=row["log2FC"],
                y=row["neglog10p"],
                text=row["gene"],
                showarrow=True
            )

        st.plotly_chart(fig, use_container_width=True)

        st.subheader("Top Differentially Expressed Genes")

        top_genes_df = volcano_df.sort_values("pval").head(10)

        st.dataframe(top_genes_df[["gene", "log2FC", "pval"]])

        # DOWNLOAD DE RESULTS
        st.download_button(
            "Download Differential Expression CSV",
            volcano_df.to_csv(index=False),
            "differential_expression.csv"
        )

        # DOWNLOAD SIGNIFICANT GENES ONLY
        sig_df = volcano_df[volcano_df["significant"]]

        st.download_button(
            "Download Significant Genes Only",
            sig_df.to_csv(index=False),
            "significant_genes.csv"
        )

# ---------------- TAB 6 ----------------

with tab6:
    st.subheader("Cell Annotation")

    marker_dict = {
        "B cells": ["MS4A1"],
        "T cells": ["CD3D"],
        "NK cells": ["NKG7"],
        "Monocytes": ["LYZ"]
    }

    cluster = st.selectbox("Cluster", sorted(adata.obs["leiden"].unique()))

    cells = adata[adata.obs["leiden"] == cluster]

    for cell_type, genes in marker_dict.items():
        genes = [g for g in genes if g in adata.var_names]
        if genes:
            score = cells[:, genes].X.mean()
            st.write(cell_type, "score:", float(score))

# ---------------- TAB 7 ----------------

with tab7:
    st.subheader("Trajectory")

    sc.tl.diffmap(adata)
    adata.uns["iroot"] = 0
    sc.tl.dpt(adata)

    df = pd.DataFrame({
        "UMAP1": adata.obsm["X_umap"][:, 0],
        "UMAP2": adata.obsm["X_umap"][:, 1],
        "pseudotime": adata.obs["dpt_pseudotime"]
    })

    fig = px.scatter(df, x="UMAP1", y="UMAP2", color="pseudotime")
    st.plotly_chart(fig)