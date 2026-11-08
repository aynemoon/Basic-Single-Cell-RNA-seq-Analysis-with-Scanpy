# Reproducibility and validation record

Date: **8 November 2026**. Audited source revision: `534bff5529fb8bbacd40f550a8db52cfd0e8aa0c`.

## What has and has not been reproduced

**The original full-data analysis has not been reproduced end to end.** A clean environment was installed, the exact input was downloaded and checked, and the full-data QC/filtering checkpoints matched the saved notebook. The full-data process was killed during Scrublet. A separate reduced-data compatibility run then executed all 36 nonempty analytical code cells successfully.

| Check | Result | Scope |
| --- | --- | --- |
| Repository and notebook audit | Passed | All 74 cells inspected; no stored exception outputs. |
| Python syntax | Passed | All original nonempty Python code cells parse; exported script has the same Python AST. |
| Clean dependency installation/imports | Passed | Isolated Linux Python 3.12.14 environment; Scanpy 1.11.5 and pinned requirements. |
| Dataset download and metadata | Passed via Figshare API | Same file ID as the original URL; SHA-256, size and library recorded in `data/README.md`. Original ndownloader route returned 403. |
| Native Jupyter execution | Blocked before analysis | Local interface/socket setup returned `Operation not permitted`; the kernel died before replying. No native Jupyter Run All success is claimed. |
| Full-data in-process execution | Partial | Original nonempty cells through one-based cell 30 completed; loading and filtering counts match the historical outputs. Exit 137 occurred during one-based cell 34, `sc.pp.scrublet(adata)`. |
| Reduced real-data execution | Passed | All 36 nonempty code cells ran in order in an in-process IPython shell on a 2,000-barcode input subset. This is compatibility testing, not a full-data reproduction. |
| Figure provenance | Passed | Four README figures decode identically to the original notebook's PNG outputs. |

## Environment

The historical notebook records **Python 3.10.19**, but does not record its package versions. The newly tested environment uses **Python 3.12.14 on Linux**. It is a compatibility baseline and cannot establish exact historical numerical reproducibility.

Direct installation requirements are in [requirements.txt](../requirements.txt), with all 136 installed distribution versions recorded in [environment-freeze.txt](environment-freeze.txt). That file is used as a constraint during installation. The workflow additionally relies on transitive packages such as Numba, UMAP-learn and pandas. No R or GPU libraries were required.

Key tested versions: Scanpy 1.11.5, AnnData 0.12.3, NumPy 2.2.6, SciPy 1.15.3, pandas 2.2.3, Matplotlib 3.10.7, Seaborn 0.13.2, scikit-learn 1.7.2, igraph 0.11.9, scikit-image 0.25.2 and h5py 3.15.1. Additional numerical-library versions are recorded in [validation-results.json](validation-results.json).

The original named kernel `xenium` need not exist on your machine. Install and select the `Scanpy intro` kernel as described in the README. Cross-platform installation and native Windows/macOS execution have not been tested.

## Full-data checkpoints

The exact HDF5 file was placed in the test working directory without changing its bytes or the notebook's loading code. With the new environment:

| Checkpoint | Fresh result | Saved historical result |
| --- | --- | --- |
| Loaded RNA matrix | 16,934 × 36,601 | 16,934 × 36,601 |
| General QC outliers | 869 | 869 |
| Mitochondrial QC outliers | 1,694 | 1,694 |
| After combined QC exclusion | 14,814 cells | 14,814 cells |
| After minimum detection filters | 14,813 × 28,471 | 14,813 × 28,471 |

The two outlier flags overlap, so their counts must not be added to infer total cells removed. Execution counts in the source notebook and one-based cell positions are different numbering systems.

The process was then killed with exit code 137 during Scrublet. No Python traceback identified a methodological error. The kill is consistent with memory exhaustion in a workspace reporting 9.7 GiB RAM and no swap; the precise peak at failure was not captured. Full-data doublet results, embeddings, clusters, markers and labels were therefore **not validated**.

## Reduced real-data compatibility run

The test harness loaded the identical file and selected 2,000 barcodes using:

```python
indices = np.sort(np.random.default_rng(0).choice(16934, 2000, replace=False))
adata = adata[indices, :].copy()
```

Only the test harness's reader return value was reduced. It then executed the original code-cell sources, unchanged and in order, in a fresh in-process IPython shell with inline Matplotlib. This avoided the workspace's Jupyter socket restriction. BLAS, OpenMP and Numba thread counts were each capped at two in the test environment; these were resource controls, not changes to the repository notebook.

- All 36 nonempty code cells completed, including plots, Scrublet scoring, PCA, neighbors, UMAP, Leiden, marker ranking, and the existing manual mapping.
- The filtered subset contained **1,747 cells × 19,965 genes**, with 2,000 variable genes, 50 PCs and two UMAP dimensions.
- Leiden resolution 0.5 yielded **six clusters**, versus seven in the historical full-data figure.
- One predicted doublet remained in the analyzed subset, confirming that no removal is implemented.
- Approximate runtime: **60.2 seconds**, after earlier imports/JIT work in the validation session; maximum process resident memory reported by `resource.getrusage`: **2,131 MiB**, including reading the full source file before subsetting. These are measurements for this smoke test, not estimates of the full-data requirements.

The subset's cluster 2 had a different marker list from historical cluster 2, while the unchanged dictionary still assigned the old label. This directly demonstrates why **successful code execution does not validate the manual annotation**. The subset's embeddings, figures and cell labels are not substituted into the repository's scientific results.

Machine-readable summary and per-cell completion records are available in [validation-results.json](validation-results.json). They are not a new biological results dataset.

## Practical expectations and remaining work

The input download is approximately 163 MiB. A CPU is sufficient, but the full workflow can require substantially more RAM than the reduced run. The full run did not complete in this 9.7 GiB workspace, so no verified full-data memory or runtime requirement can be supplied. Use a machine with additional memory for a full rerun and monitor Scrublet's memory use.

Remaining validation includes a native Jupyter Run All test, full-data completion on a larger-memory machine, comparisons of regenerated clusters and markers against the historical outputs, and biological review of every manual cell-type assignment. Explicit parameter/default pinning and changes to reduce memory use should be reviewed before modifying the original analysis.

The original notebook and its outputs remain unchanged. A new README and reading guide explain limitations rather than silently correcting scientific behavior. In particular, no doublet removal, new marker assignments, or replacement of the original Figshare URL was pushed into notebook code.
