# signal-overlap-violin-stats

Statistics and figure-generation code for Figure panel C: SMARCA5 / DMRT1 ChIP peak
signal intensity (signalValue), stratified by whether each peak overlaps an
RAR / DMRT1 / SMARCA5 peak.

## What's here

- `scripts/analyze_signal_overlap.py` — loads the per-comparison signalValue data,
  runs a two-sided Mann–Whitney U test for each of the four comparisons, and
  re-draws the violin + box + jittered-scatter figure with the resulting p-values.
- `requirements.txt` — Python dependencies.

## What's *not* here

The upstream step that assigns each peak to an "overlap" / "no overlap" group is
**not** included, because this repo only has access to the already-split
signalValue numbers, not the original peak (BED/narrowPeak) files. That step was:

```
bedtools closest -a <factor>.narrowPeak -b <partner>.narrowPeak -d \
    > <factor>_vs_<partner>.closest.bed
# last column == 0  -> "overlap" group   (peaks share >= 1 bp)
# last column >  0  -> "no overlap" group
```

If you have the original narrowPeak files, add a small script that runs the above
and writes a two-column signalValue table in the same layout as the input workbook
described below, and this analysis script will consume it unchanged.

## Input format

An Excel workbook with one sheet per comparison (`SMARCA5 vs RAR`,
`SMARCA5 vs DMRT1`, `DMRT1 vs RAR`, `DMRT1 vs SMARCA5`), metadata/description text
in rows 1–3, a header row at row 4, and signalValue data from row 5 onward in
columns A (group A) and B (group B). The two columns may have different numbers of
rows.

## Usage

```bash
pip install -r requirements.txt
python scripts/analyze_signal_overlap.py \
    --input violin_plot_underlying_data.xlsx \
    --outdir ./output
```

Produces `output/stats_summary.csv` (U, z, p, and a log10(p)/mantissa/exponent
breakdown for p-values that underflow double precision) and
`output/figure_C.png` / `output/figure_C.pdf`.

## Notes on extremely small p-values

For two of the four comparisons, the Mann–Whitney asymptotic p-value underflows to
exactly `0.0` in double precision (SciPy's normal-approximation p-value computation
loses precision once the z-score gets large enough — here z ≈ -50 and z ≈ -39). This
script recovers a finite estimate by working in log-space with
`scipy.stats.norm.logsf`, then reports the figure annotation as `P < 1e-300` rather
than the literal (and statistically over-precise) mantissa/exponent, since values at
that extreme are more a statement of "the asymptotic approximation says this is
enormously significant" than a trustworthy 3-significant-figure p-value.
