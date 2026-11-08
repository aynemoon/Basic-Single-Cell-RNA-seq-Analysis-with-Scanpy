# Read-only audit of the original repository

Audit date: **8 November 2026**. Source revision: [`534bff5529fb8bbacd40f550a8db52cfd0e8aa0c`](https://github.com/aynemoon/Basic-Single-Cell-RNA-seq-Analysis-with-Scanpy/commit/534bff5529fb8bbacd40f550a8db52cfd0e8aa0c).

## Scope and inventory

At that revision the entire tracked repository consisted of one notebook: `Basic Single-Cell RNA-seq Analysis with Scanpy (2).ipynb`, 8,157,270 bytes. No README, dependency specification, datasets, standalone scripts, license, or result files were present.

All **74 cells** were inspected: 37 markdown cells and 37 code cells, including one empty final cell. The 36 nonempty code cells have saved execution counts 1–36 in order. Saved outputs comprise 18 PNG figure displays, eight expression/table displays, and four stream outputs. There are no stored exception outputs. The saved run therefore contains no recorded exception, but that alone does not establish clean-environment reproducibility.

- Original notebook SHA-256: `47a38d879454ab5a21f7e7308119eea153d5d29a30ac9470469b0ea5f70f35b0`.
- Original Git blob: `af8ff5fb83110404ac3a83f910e3b0212068667a`.
- [Cell-by-cell inventory](cell-inventory.json): one-based position, cell type, execution count, source checksum, source preview and saved-output summary for every cell. Embedded base64 content is abbreviated in previews, not removed from the source notebook.
- [Companion guide](NOTEBOOK_GUIDE.md): scientific stages and interpretation notes.

## Confirmed analysis

The implemented pipeline and exact explicit thresholds are listed in the [README workflow table](../README.md#analysis-workflow). They were read from code, not inferred from introductory prose. In particular:

- The loading cell reads a relative filename and provides a Figshare backup download URL. The default RNA-only reader discards ATAC peaks from the multiome file.
- `.X` initially stores sparse counts. The notebook disambiguates gene names, calculates QC covariates, excludes flagged cells and applies minimum detection filters.
- The saved shapes are 16,934 × 36,601 at loading, 14,814 × 36,601 after outlier exclusion, and 14,813 × 28,471 after minimum detection filtering.
- Scrublet scores and predictions are calculated and plotted. There is no code removing predicted doublets.
- The filtered count matrix is copied to `.layers["counts"]`, followed by median-depth total-count normalization and `log1p` transformation.
- Two thousand variable genes are requested; the matrix is not explicitly subset to them. PCA, neighbors and UMAP use unspecified/default settings beyond the explicit plotting parameters.
- Initial Leiden clustering specifies `flavor="igraph", n_iterations=2`. The later resolution loop uses 0.02, 0.5 and 1.0 with `flavor="igraph"` and the default iteration count.
- Marker ranking uses Wilcoxon testing at resolution 0.5. Cluster 2 is inspected and a hard-coded annotation dictionary is applied.
- No scaling, covariate regression, ambient RNA correction, batch integration, trajectory inference, t-SNE, pseudobulk differential expression, automated annotation, or output-data saving is implemented.

## Figures and outputs

The notebook includes count histograms, mitochondrial violin plots and QC scatter plots before/after outlier filtering, a Scrublet score-distribution plot, highly-variable-gene plots, PCA variance and coordinate plots, uncolored and cluster-colored UMAPs, doublet/QC UMAP overlays, a three-resolution UMAP comparison, a marker dotplot, ten cluster-2 gene-expression UMAP panels plus a cluster panel, and a manually labeled cell-type UMAP.

Six markdown cells contain embedded illustrations/screenshots: AnnData structure, AnnData interoperability, Scanpy API, doublet detection, and two Enrichr screenshots. These are distinct from the 18 analysis figure outputs.

The largest stored analysis image is the marker-expression UMAP output in one-based cell 66, approximately 2.56 million base64 characters. The original notebook is approximately 7.78 MiB. The large images serve a teaching purpose; they have not been cleared or recompressed. Selected README plots are exact extractions, with provenance recorded in [figures/README.md](../figures/README.md).

## Environment and portability

| Finding | Evidence and implication |
| --- | --- |
| Historical Python known | Notebook metadata records Python 3.10.19. |
| Historical packages unknown | There is no saved package-version listing or environment file. |
| Local kernel name | Metadata names `xenium`; another user must select an installed kernel. This does not mean Xenium software is needed. |
| Relative input path | The HDF5 file must be in the working directory. There are no absolute filesystem paths in analytical code. |
| Personal paths in saved warnings | Two duplicate-gene warnings and a Numba threading warning contain the original environment's home-directory prefix. They are output text, not required paths. |
| Missing input | The dataset is absent from Git; downloading it is necessary. No other external local input files are referenced. |
| Broken download route in this workspace | Original Figshare route returned 403; the API route for the identical file ID worked. |
| Direct imports | NumPy, Scanpy, Seaborn and SciPy's `median_abs_deviation` are all used. |
| Functional dependencies | AnnData, h5py, Matplotlib, pandas, scikit-learn, UMAP/Numba and igraph support the workflow; scikit-image is needed by automatic Scrublet thresholding in the tested environment. |
| No unnecessary advanced tools | PyTorch is mentioned in explanatory text; PyDESeq2 in the final caution. Neither is imported or needed. |
| Random/default parameters | No global seed or explicit seed arguments are set. Several Scanpy functions have their own seeded defaults; recording versions still matters. |

## Scientific and teaching issues

The companion guide clarifies genes versus transcripts, counts versus reads, two-sided MAD rules, computed versus removed doublets, variable-gene masking, and the meaning of the dotplot scale. It also identifies a duplicated PCA paragraph and the inherited `dc_cluster_genes` variable name. The manual cell-type mapping is tied to historical cluster IDs and has no external validation record. The Enrichr screenshots do not establish a reproducible automated annotation procedure.

The final recommendation to consider PyDESeq2 is incomplete without biological replication and an appropriate design. This one-library exercise supports exploratory cluster-marker interpretation, not replicated condition-level differential-expression claims.

## External tutorial attribution

The opening cell already credits Single-cell Best Practices and identifies the Fall 2025 teaching context. Direct comparison also confirms close correspondence in text/code with:

- [Single-cell Best Practices, Quality Control](https://www.sc-best-practices.org/preprocessing-visualization/quality-control/): source example, QC metrics, the MAD helper and outlier filters.
- [Scanpy, Preprocessing and clustering](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html), including its [1.10-series version](https://scanpy.readthedocs.io/en/1.10.x/tutorials/basics/clustering.html): normalization, feature selection, PCA, graph analysis, Leiden, QC reassessment and marker-ranking examples. Similar variable names and prose establish adaptation; this is not just a generic similarity of methods.

The repository should credit both sources explicitly and describe Aysan's work as an educational adaptation. The new README does so. The current source tutorials have evolved; their newer data-loading and methodological additions are not automatically part of this notebook.

Individual provenance and reuse terms for the embedded diagrams and screenshots are incompletely documented. These images remain untouched in the original notebook and are not newly extracted for the README. No blanket repository license is added because original and third-party rights require separate review.

## Proposed structure and approval boundary

This update adds a root README, requirements and ignore rules, `data/README.md`, `scripts/scanpy_intro.py`, selected original figures, and documentation under `docs/`. The original root notebook is preserved byte-for-byte, and existing Git history is retained.

Moving the notebook to `notebooks/` would require the owner's approval and coordinated path/documentation changes. It offers little benefit for a single-notebook project, so it is deferred. The proposed notebook edits in the reading guide are likewise deferred for review. This update does not replace historical outputs or silently alter analysis behavior.
