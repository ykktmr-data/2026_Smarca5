# Spike-in normalized differential analysis of H3K27me3 signal (pyDESeq2)

Differential analysis of H3K27me3 signal between Ctrl and cKO at the union of
MACS2 narrow peaks, normalized with spike-in (dm6) scale factors, with MA plot
and per-category BED output.

## Why spike-in size factors

DESeq2's default median-of-ratios normalization assumes that most regions do
not change. If a mark changes globally between conditions, that assumption
cancels the real change. Here, spike-in-derived scale factors are used
directly as DESeq2 size factors:

```
size_factor = (library size in millions) / scaleFactor
```

This converts the deepTools convention (`bamCoverage --normalizeUsing RPKM
--scaleFactor X`, i.e. signal = RPKM x X) to the DESeq2 convention
(normalized = raw / size_factor).

Scale factors are condition-level (computed on merged replicates), so both
replicates of a condition share the same factor; library size is per
replicate.

## System requirements

- **OS**: developed and run by the authors on macOS Sonoma 14.1. The demo
  below was additionally verified end-to-end on Linux (Ubuntu 24.04,
  aarch64) as part of preparing this repository.
- **Python**: >= 3.9; authors' environment is Python 3.12.3. The demo was
  also verified on Python 3.11.15 with pydeseq2 0.5.4.
- **Python dependencies** (see `requirements.txt`): numpy, pandas,
  matplotlib, pydeseq2.
- **Non-Python dependency**: `bedtools` CLI on PATH (tested with bedtools
  2.31.1).
- **Hardware**: no non-standard hardware required. Runs on a standard
  desktop/laptop CPU.

The script sets size factors through pydeseq2 internals
(`dds.obs["size_factors"]`, `dds.var["_normed_means"]`), which are not public
API. Pin the pydeseq2 version you used in `requirements.txt` and check the
size-factor assertion in the log after upgrading.

## Installation guide

```bash
pip install -r requirements.txt
```

Also install bedtools, e.g. `apt install bedtools` or
`conda install -c bioconda bedtools`.

Typical install time: ~30–60 seconds on a standard broadband connection for
the Python dependencies (measured ~30 s in a clean virtual environment;
pydeseq2 is the largest dependency), plus ~10 s–1 min for bedtools depending
on package manager.

## Demo

```bash
python make_demo_data.py --outdir ./demo_data
python k27_deseq2_dm6_ma.py \
  --peaks-ctrl demo_data/Ctrl_merged_peaks.narrowPeak \
  --peaks-cko  demo_data/cKO_merged_peaks.narrowPeak \
  --beds-ctrl  demo_data/Ctrl_1.rmDup.bed demo_data/Ctrl_2.rmDup.bed \
  --beds-cko   demo_data/cKO_1.rmDup.bed demo_data/cKO_2.rmDup.bed \
  --scale-ctrl 1.0 --scale-cko 1.0 \
  --outdir ./demo_results
```

`make_demo_data.py` generates 60 synthetic peaks on a single fake chromosome,
with 8 peaks given a ~2.5–4x higher count in "cKO" and 8 peaks given a
~0.2–0.4x lower count, so the pipeline has real up/down signal to detect (not
actual H3K27me3 data).

Expected output (in `demo_results/`): `union_peaks.bed`, `raw_counts.tsv`,
`deseq2_results.tsv` (baseMean, log2FoldChange, lfcSE, pvalue, padj per
peak), `up_in_cKO.bed` / `down_in_cKO.bed` / `not_significant.bed`, and
`ma_plot.png` / `.pdf`. On the synthetic data this reproducibly calls close
to 8 peaks up, 8–9 down, and the rest not significant.

Expected run time: ~7 seconds on a standard laptop for 60 peaks / 4
replicates. Runtime scales with the number of union peaks and total fragment
count on real genome-scale data (typically a few minutes).

## Instructions for use — running on your own data

Inputs: narrowPeak files from MACS2 (merged replicates per condition) and
per-replicate fragment BED files (`bedtools bamtobed` on the duplicate-removed
BAMs).

```bash
python3 k27_deseq2_dm6_ma.py \
  --peaks-ctrl Ctrl_merged_peaks.narrowPeak \
  --peaks-cko  cKO_merged_peaks.narrowPeak \
  --beds-ctrl  Ctrl_1.rmDup.bed Ctrl_2.rmDup.bed \
  --beds-cko   cKO_1.rmDup.bed  cKO_2.rmDup.bed \
  --scale-ctrl 0.612 \
  --scale-cko  0.675 \
  --outdir results
```

### Outputs (`--outdir`)

| File | Content |
|---|---|
| `union_peaks.bed` | Union of Ctrl and cKO peaks |
| `raw_counts.tsv` | Raw fragment counts per peak per replicate |
| `deseq2_results.tsv` | baseMean, log2FoldChange, lfcSE, pvalue, padj per peak |
| `up_in_cKO.bed`, `down_in_cKO.bed`, `not_significant.bed` | Peaks by category (padj < threshold and sign of log2FC) |
| `ma_plot.pdf`/`.png` | MA plot (Okabe-Ito colorblind-safe colors) |

Peaks with NaN padj (removed by DESeq2 independent filtering) are excluded
from the three BED files.

## License

Released under the MIT License — see [LICENSE](../LICENSE).
