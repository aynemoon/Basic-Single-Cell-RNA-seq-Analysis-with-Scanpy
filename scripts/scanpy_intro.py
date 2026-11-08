# Searchable Python export of the original educational notebook.
# Source commit: 534bff5529fb8bbacd40f550a8db52cfd0e8aa0c
# Attribution: Single-cell Best Practices QC chapter and Scanpy clustering tutorial.
# See README.md and docs/AUDIT.md for sources, assumptions and limitations.
# Code cells are unchanged and retain their original execution order.
# Use the notebook for explanations, inline plots and displayed tables.
# If run as a script, start from the repository root; no results are saved automatically.

# %% Original notebook cell 8 (one-based)
import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation

sc.settings.verbosity = 0
sc.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

# %% Original notebook cell 9 (one-based)
adata = sc.read_10x_h5(
    filename="filtered_feature_bc_matrix.h5",
    backup_url="https://figshare.com/ndownloader/files/39546196",
)
adata

# %% Original notebook cell 11 (one-based)
adata.var_names_make_unique()
adata

# %% Original notebook cell 15 (one-based)
# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

# %% Original notebook cell 17 (one-based)
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
)
adata

# %% Original notebook cell 19 (one-based)
p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
# sc.pl.violin(adata, 'total_counts')
p2 = sc.pl.violin(adata, "pct_counts_mt")
p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# %% Original notebook cell 21 (one-based)
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

# %% Original notebook cell 23 (one-based)
adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)
adata.obs.outlier.value_counts()

# %% Original notebook cell 25 (one-based)
adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)
adata.obs.mt_outlier.value_counts()

# %% Original notebook cell 27 (one-based)
print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

# %% Original notebook cell 28 (one-based)
p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
# sc.pl.violin(adata, 'total_counts')
p2 = sc.pl.violin(adata, "pct_counts_mt")
p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
adata

# %% Original notebook cell 30 (one-based)
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
adata

# %% Original notebook cell 34 (one-based)
sc.pp.scrublet(adata)

# %% Original notebook cell 35 (one-based)
sc.pl.scrublet_score_distribution(adata)

# %% Original notebook cell 37 (one-based)
# Saving count data
adata.layers["counts"] = adata.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

# %% Original notebook cell 39 (one-based)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# %% Original notebook cell 40 (one-based)
sc.pl.highly_variable_genes(adata)

# %% Original notebook cell 42 (one-based)
sc.tl.pca(adata)

# %% Original notebook cell 45 (one-based)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# %% Original notebook cell 47 (one-based)
sc.pl.pca(
    adata,
    color=["pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3)],
    ncols=2,
    size=2,
)

# %% Original notebook cell 49 (one-based)
sc.pp.neighbors(adata)

# %% Original notebook cell 51 (one-based)
sc.tl.umap(adata)

# %% Original notebook cell 52 (one-based)
sc.pl.umap(adata)

# %% Original notebook cell 54 (one-based)
# Using the igraph implementation and a fixed number of iterations can be significantly faster,
# especially for larger datasets
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

# %% Original notebook cell 55 (one-based)
sc.pl.umap(adata, color=["leiden"])

# %% Original notebook cell 57 (one-based)
sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
)

# %% Original notebook cell 58 (one-based)
sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
)

# %% Original notebook cell 59 (one-based)
for res in [0.02, 0.5, 1.0]:
    sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph")

# %% Original notebook cell 60 (one-based)
sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_1.00"],
    legend_loc="on data",
)

# %% Original notebook cell 62 (one-based)
# Obtain cluster-specific differentially expressed genes
sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")

# %% Original notebook cell 64 (one-based)
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=10)

# %% Original notebook cell 65 (one-based)
sc.get.rank_genes_groups_df(adata, group="2").head(10)

# %% Original notebook cell 66 (one-based)
dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group="2").head(10)["names"]
sc.pl.umap(
    adata,
    color=[*dc_cluster_genes, "leiden_res_0.50"],
    legend_loc="on data",
    frameon=False,
    ncols=3,
)

# %% Original notebook cell 68 (one-based)
df = sc.get.rank_genes_groups_df(adata, group="2")

names_str = ", ".join(df["names"].head(10).astype(str).tolist())

print(names_str)


# %% Original notebook cell 71 (one-based)
cluster_to_name = {
    "0": "CD8+ T",
    "1": "CD4+ T",
    "2": "Erythroid",
    "3": "B cells",
    "4": "B cells",
    "5": "Monocytes",
    "6": "B cells",
}

adata.obs["celltype"] = (
    adata.obs["leiden_res_0.50"].map(cluster_to_name).astype("category")
)

# %% Original notebook cell 72 (one-based)
sc.pl.umap(adata, color="celltype", legend_loc="on data")

