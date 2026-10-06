# signal-overlap-violin-stats

Statistics and figure-generation code for Figure panel C: SMARCA5 / DMRT1 ChIP peak
signal intensity (signalValue), stratified by whether each peak overlaps an
RAR / DMRT1 / SMARCA5 peak.

> This is a standalone, single-panel copy of the same analysis also included
> in the `Fig.S9b` folder of this repository (where it appears alongside the
> Panel B summit-proximity analysis). It's kept here as its own self-contained
> unit so the code for this specific supplementary figure panel can be run
> without the Panel B dependencies (`bedtools`). The two `analyze_signal_overlap.py`
> files are identical.

## System requirements

- **OS**: developed and run by the authors on macOS Sonoma 14.1. The demo
  below was additionally verified end-to-end on Linux (Ubuntu 24.04,
  aarch64) as part of preparing this repository.
- **Python**: 3.12.3 (authors' environment); the demo was also verified on
  Python 3.11.15.
- **Dependencies** (see `requirements.txt`): numpy>=1.24, pandas>=2.0,
  scipy>=1.10, matplotlib>=3.7, openpyxl>=3.1.
- **Hardware**: no non-standard hardware required. Runs on a standard
  desktop/laptop CPU.

## Installation guide

```bash
pip install -r requirements.txt
```

Typical install time: well under 1 minute on a standard broadband connection
(measured ~35 s for these 5 packages in a clean virtual environment).

## Demo

```bash
python examples/make_demo_violin_data.py --outdir ./demo_data
python scripts/analyze_signal_overlap.py \
    --input ./demo_data/violin_plot_underlying_data.xlsx --outdir ./demo_out
```

Expected output: `demo_out/stats_summary.csv` (one row per comparison, with
U, z and p) and `demo_out/figure_C.png` / `.pdf` (4-panel violin/box plot).
The synthetic generator gives the "overlap" group a ~1.5–2x larger effect
size than the "no overlap" group, so most of the four comparisons should come
back significant (P < 0.05) — exact values differ run to run since the data
are randomly generated.

Expected run time: ~5 seconds on a standard laptop.

## What's here

- `scripts/analyze_signal_overlap.py` — loads the per-comparison signalValue data,
  runs a two-sided Mann–Whitney U test for each of the four comparisons, and
  re-draws the violin + box + jittered-scatter figure with the resulting p-values.
- `examples/make_demo_violin_data.py` — generates a small synthetic signalValue
  workbook for the demo above.
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

## Instructions for use — running on your own data

Input format: an Excel workbook with one sheet per comparison (`SMARCA5 vs RAR`,
`SMARCA5 vs DMRT1`, `DMRT1 vs RAR`, `DMRT1 vs SMARCA5`), metadata/description text
in rows 1–3, a header row at row 4, and signalValue data from row 5 onward in
columns A (group A) and B (group B). The two columns may have different numbers of
rows.

```bash
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

## License

Released under the MIT License — see [LICENSE](../LICENSE).
