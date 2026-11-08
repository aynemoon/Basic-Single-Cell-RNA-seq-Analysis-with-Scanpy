# Dataset access and provenance

[Project overview](../README.md)

## Biological source

The notebook analyzes the RNA modality of a 10x Multiome library from healthy human bone marrow mononuclear cells, from the NeurIPS 2021 multimodal single-cell integration benchmark. The upstream [Single-cell Best Practices QC chapter](https://www.sc-best-practices.org/preprocessing-visualization/quality-control/) describes this example as sample 4 of donor 8.

The downloaded HDF5 file independently identifies its library as **`site4_donor8`**. Its attributes record `cellranger-arc-2.0.0` and chemistry `ARC-737K-arc-v1`. It contains **190,628 total features × 16,934 barcodes** in 10x storage orientation, with feature types `Gene Expression` and `Peaks`.

The notebook's `sc.read_10x_h5()` call uses the default `gex_only=True`, returning a cell-by-gene AnnData object with **16,934 × 36,601** RNA entries along its two axes. No chromatin-accessibility analysis is performed.

## Retrieve the exact input

The source notebook references Figshare file **39546196**, using:

`https://figshare.com/ndownloader/files/39546196`

On 8 November 2026, that route returned HTTP 403 in the validation workspace. Figshare's public API download route for the **same file ID** succeeded:

`https://api.figshare.com/v2/file/download/39546196`

From the repository root, in the activated Python environment:

```bash
python -c "from urllib.request import urlretrieve; urlretrieve('https://api.figshare.com/v2/file/download/39546196', 'filtered_feature_bc_matrix.h5')"
```

This deliberately saves the file beside the original notebook. The existing loading cell looks in its current working directory; placing it inside this `data/` folder would require a code-path change. Once the local file exists, Scanpy does not need to use the notebook's backup URL.

Verify the downloaded file:

```bash
python -c "from pathlib import Path; import hashlib; p=Path('filtered_feature_bc_matrix.h5'); print(p.stat().st_size); print(hashlib.sha256(p.read_bytes()).hexdigest())"
```

Expected values from the verified download:

| Property | Value |
| --- | --- |
| Size | 170,659,372 bytes, approximately 163 MiB |
| SHA-256 | `aae39275941b2c63adf3a52d8f82127bcdfb67cd6e1abc9f7d2c906a5c1b5e50` |
| Library | `site4_donor8` |

If the checksum or library differs, resolve the discrepancy before running the analysis. Download availability can change. Do not silently substitute PBMC3k, another donor, or the different libraries used by the current Scanpy tutorial.

## Scope and reuse

The parent benchmark includes 12 healthy donors measured at four sites, according to the upstream QC chapter. Only this one library is loaded here, and the notebook does not retain donor/site comparisons as an analysis design. This is not the complete multiome benchmark and does not support conclusions about all donors.

The dataset is downloaded from its provider rather than redistributed in this repository. Follow the original study and hosting provider's applicable terms, and cite Luecken et al. (2021), [A sandbox for prediction and integration of DNA, RNA, and proteins in single cells](https://openreview.net/forum?id=gN35BGa1Rt), together with the teaching sources. A checksum establishes which bytes were tested; it does not establish ownership or a reuse license.
