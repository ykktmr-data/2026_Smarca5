# Spike-in normalized differential analysis of H3K27me3 signal (pyDESeq2)

Differential analysis of H3K27me3 signal between Ctrl and cKO at the union of MACS2 narrow peaks, normalized with spike-in (dm6) scale factors, with MA plot and per-category BED output.

## Why spike-in size factors

DESeq2's default median-of-ratios normalization assumes that most regions do not change. If a mark changes globally between conditions, that assumption cancels the real change. Here, spike-in-derived scale factors are used directly as DESeq2 size factors:

```
size_factor = (library size in millions) / scaleFactor
```

This converts the deepTools convention (`bamCoverage --normalizeUsing RPKM --scaleFactor X`, i.e. signal = RPKM x X) to the DESeq2 convention (normalized = raw / size_factor).

Scale factors are condition-level (computed on merged replicates), so both replicates of a condition share the same factor; library size is per replicate.

## Requirements

- `bedtools` on PATH
- Python >= 3.9: `pip install -r requirements.txt`

The script sets size factors through pydeseq2 internals (`dds.obs["size_factors"]`, `dds.var["_normed_means"]`), which are not public API. Pin the pydeseq2 version you used in `requirements.txt` and check the size-factor assertion in the log after upgrading.

## Usage

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

Inputs: narrowPeak files from MACS2 (merged replicates per condition) and per-replicate fragment BED files (`bedtools bamtobed` on the duplicate-removed BAMs).

## Outputs (`--outdir`)

| File | Content |
|---|---|
| `union_peaks.bed` | Union of Ctrl and cKO peaks |
| `raw_counts.tsv` | Raw fragment counts per peak per replicate |
| `deseq2_results.tsv` | baseMean, log2FoldChange, lfcSE, pvalue, padj per peak |
| `up_in_cKO.bed`, `down_in_cKO.bed`, `not_significant.bed` | Peaks by category (padj < threshold and sign of log2FC) |
| `ma_plot.pdf/.png` | MA plot (Okabe-Ito colorblind-safe colors) |

Peaks with NaN padj (removed by DESeq2 independent filtering) are excluded from the three BED files.
