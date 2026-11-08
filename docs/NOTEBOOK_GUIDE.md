# Companion guide to the original notebook

[Return to the project overview](../README.md)

The original notebook is retained at its existing root path, without code, output, or markdown changes. This guide adds navigation and scientific clarifications without rewriting the historical teaching material. Cell numbers below are **one-based positions**, including markdown cells, in source commit `534bff5529fb8bbacd40f550a8db52cfd0e8aa0c`.

| Cells | Stage | What to inspect |
| --- | --- | --- |
| 1–5 | AnnData and Scanpy concepts | Cells × genes in `.X`; annotations in `.obs` and `.var`; Scanpy's `pp`, `tl`, and `pl` namespaces. |
| 6–11 | Motivation, imports and loading | The file download, duplicate gene symbols, and initial 16,934 × 36,601 shape. |
| 12–28 | QC metrics and MAD filtering | Compare before/after plots; inspect `outlier` and `mt_outlier`. |
| 29–31 | Minimum detection filters | Retain cells with ≥100 genes and genes detected in ≥3 cells. |
| 32–35 | Doublet scoring | Inspect the score distribution; no doublet-removal statement is implemented. |
| 36–40 | Normalization and variable genes | Save filtered counts, normalize to median library size, apply `log1p`, and mark 2,000 variable genes. |
| 41–47 | PCA | Examine the variance-ratio plot and mitochondrial coloring. |
| 48–55 | Neighbor graph, UMAP and initial Leiden clustering | Clustering uses the graph; UMAP is a visualization of its relationships. |
| 56–60 | QC reassessment and clustering resolution | Compare doublet/QC overlays and resolutions 0.02, 0.5, and 1.0. |
| 61–66 | Marker ranking | Wilcoxon cluster-versus-rest tests, top-10 dotplot, cluster-2 table, and gene-expression UMAPs. |
| 67–72 | Manual interpretation | Marker-list export, external Enrichr screenshots, and a provisional cluster-name dictionary. |
| 73–74 | Statistical caution and trailing empty cell | Read the clarification on replication below. |

## Important clarifications

**Genes, not transcript molecules.** The matrix columns represent gene features. Counts are measured molecules/UMIs summarized by gene, rather than one column per individual transcript. A barcode is not guaranteed to represent exactly one intact cell.

**Zeros have several causes.** A zero can reflect limited sampling or genuinely absent/low expression. It should not automatically be interpreted as a technical dropout requiring imputation.

**QC percentages refer to counts.** `pct_counts_mt` is the percentage of a barcode's counts assigned to mitochondrial genes. The notebook sometimes calls these reads; the implemented metric is a count fraction. The rules here are a teaching example and need biological review for another tissue or protocol.

**MAD logic is two-sided.** `is_outlier()` flags values both below and above the median-based interval. Mitochondrial filtering therefore combines a two-sided 3-MAD rule with an additional >8% upper threshold. Ribosomal and hemoglobin percentages are calculated but are not separate filtering criteria.

**Doublet scoring is not removal.** Scrublet adds predictions and scores. They remain in the data for every later analysis. The introduction discusses removal as a general option; that option is not executed here.

**Variable genes are marked, not physically subset.** The default highly-variable-gene call adds annotations while retaining other genes. In the validation Scanpy version, PCA uses the variable-gene mask by default. Marker ranking can still inspect the full normalized gene matrix. No gene-wise scaling or regression step is implemented.

**PCA explanations are duplicated.** Cells 43 and 44 repeat the same explanation. Cell 47 deliberately repeats `pct_counts_mt` to color two different PC pairs. Neither changes the data. Claims that increasing PC count has no downside should not be generalized: the retained components and neighbor settings should be assessed for each dataset.

**Check cluster identities again after a rerun.** `dc_cluster_genes` in cell 66 is only a variable name inherited from an example. It does not establish a dendritic-cell identity. The code actually fetches cluster 2, which is labeled erythroid in the saved manual dictionary. Numerical cluster IDs and their marker lists can change after changes in dependencies or analysis settings.

**The dotplot colors are gene-wise rescaled.** With `standard_scale="var"`, colors emphasize within-gene differences across clusters; they cannot be used as an absolute comparison of expression between genes. The dot sizes represent fractions of cells with detected expression.

**P-values do not replace biological replication.** The notebook's last paragraph mentions PyDESeq2. Merely switching software does not make cell-level tests replicate-aware. A donor/condition comparison normally needs an appropriate design and independent biological samples, often with count aggregation by sample and cell type. This single-library teaching example does not implement such an analysis.

## Suggested future notebook edits requiring review

Possible small follow-up edits are to add these clarifications next to the relevant cells, remove the repeated PCA paragraph, add a clear attribution block, and rename `dc_cluster_genes` to `cluster_2_genes`. Changes to download paths, explicit random seeds/default parameters, doublet removal, cell-type assignments, or notebook location should be reviewed separately because they can alter behavior or interpretation. No such changes are included in this update.
