# smarca5-dmrt1-rar-colocalization

Statistics and figure-generation code for two related panels examining
co-localization of SMARCA5, DMRT1 and RAR binding sites:

- **Figure panel C** — SMARCA5 / DMRT1 ChIP peak signal intensity (signalValue),
  stratified by whether each peak overlaps an RAR / DMRT1 / SMARCA5 peak
  (`scripts/analyze_signal_overlap.py`).
- **Figure panel B** — distance from each peak summit to the nearest summit of
  a partner factor: raw distance histograms (0–5 kb) and fold-enrichment over
  a chromosome-restricted random shuffle (0–20 kb)
  (`scripts/proximity_summit_analysis.py`).

## System requirements

- **OS**: developed and run by the authors on macOS Sonoma 14.1. The demos
  below were additionally verified end-to-end on Linux (Ubuntu 24.04,
  aarch64) as part of preparing this repository. No OS-specific code paths
  are used, so any recent macOS, Linux or WSL install should work.
- **Python**: 3.12.3 (authors' environment). The demos were verified to also
  run correctly on Python 3.11.15.
- **Python dependencies** (see `requirements.txt`): numpy>=1.24, pandas>=2.0,
  scipy>=1.10, matplotlib>=3.7, openpyxl>=3.1.
- **Non-Python dependency**: `proximity_summit_analysis.py` additionally
  requires the `bedtools` CLI on PATH (tested with bedtools 2.31.1;
  https://bedtools.readthedocs.io/).
- **Hardware**: no non-standard hardware required. Runs on a standard
  desktop/laptop CPU; no GPU needed.

## Installation guide

```bash
pip install -r requirements.txt
```

Typical install time: well under 1 minute on a standard broadband connection
(measured ~35 s for these 5 packages in a clean virtual environment). If you
need `proximity_summit_analysis.py`, also install bedtools, e.g.
`apt install bedtools` or `conda install -c bioconda bedtools`
(typically 10 s–1 min depending on package manager).

## Demo

Two independent demos are provided (one per panel), both runnable end-to-end
on synthetic data — no real sequencing files required.

**Panel C — signal intensity vs. overlap**

```bash
python examples/make_demo_violin_data.py --outdir ./demo_data
python scripts/analyze_signal_overlap.py \
    --input ./demo_data/violin_plot_underlying_data.xlsx --outdir ./demo_out_C
```

Expected output: `demo_out_C/stats_summary.csv` (one row per comparison, with
U, z and p) and `demo_out_C/figure_C.png` / `.pdf` (4-panel violin/box plot).
The synthetic generator gives the "overlap" group a ~1.5–2x larger effect
size than the "no overlap" group, so most of the four comparisons should come
back significant (P < 0.05) — exact values differ run to run since the data
are randomly generated.

Expected run time: ~5 seconds on a standard laptop.

**Panel B — summit distance + enrichment**

```bash
python examples/make_synthetic_test_data.py --outdir ./test_data
cd test_data && python ../scripts/proximity_summit_analysis.py
```

Expected output: `test_data/proximity_out/{RAR,SMARCA5,DMRT1}_ref_SUMMIT_50bp_EN.png`
and `.pdf` (one 3-part composite figure per reference factor: two 0–5 kb
summit-distance histograms plus the 0–20 kb enrichment-vs-random curve), plus
per-pair `CLOSEST_*.tsv` distance tables.

Expected run time: ~7 seconds on a standard laptop, for `N_SHUFFLES = 20` on
the synthetic peak sets (~250–1,000 peaks per factor). On real genome-scale
peak files this scales up (one `bedtools shuffle` + `bedtools closest` per
shuffle per factor pair) and can take several minutes.

## Status of each script

| Script | Verified against real data? |
|---|---|
| `analyze_signal_overlap.py` | **Yes.** Run against the actual signalValue workbook; its output (statistics table + figure) was checked and matches the published panel C. |
| `proximity_summit_analysis.py` | **Supplied as the original analysis code**, not reconstructed from the figure. Confirmed here to run end-to-end without errors (see `examples/`) against synthetic stand-in peak files (the real narrowPeak/chrom.sizes files were not available in this environment). Its enrichment-curve *shape* is independently confirmed against `SMARCA5_summit_enrichment_standalone_EN.pdf`, a standalone figure this script actually produced — see **Provenance vs. the published figure** below for what that does and doesn't confirm. |

## proximity_summit_analysis.py — method

1. Extracts each peak's summit (1 bp) from its narrowPeak file using column 10
   (the MACS2 summit offset from `start`).
2. For every peak of factor A, finds the nearest summit of factor B with
   `bedtools closest -d -t first`; peaks with no partner on their chromosome
   (`-d` reports `-1`) are dropped.
3. Builds a null background by shuffling factor B's summits *within their
   chromosome of origin* (`bedtools shuffle -chrom`), `N_SHUFFLES` times, and
   repeating step 2 against each shuffle.
4. For a grid of distance thresholds, reports the cumulative percentage of
   factor-A peaks within that threshold of a factor-B peak, the same
   percentage under the shuffled background, and their ratio
   (`enrichment = observed_pct / random_pct`) — this is a **cumulative**
   (`distance <= threshold`) definition, not a local/windowed density ratio.
5. Renders one summary figure per reference factor (`{REF}_ref_SUMMIT_50bp_EN.png/.pdf`):
   two 50-bp-binned distance histograms (0–5 kb) to its two partner factors,
   plus the enrichment-vs-random curve (0–20 kb, log scale).

Peak-calling parameters encoded in the original file names: RAR and SMARCA5
peaks were called with MACS2 at p < 1e-4, DMRT1 at p < 1e-5; genome build mm10.

This is an **edit-and-run** script, not a CLI: paths and parameters live in
the `CONFIG` block near the top (`OUT_DIR`, `GENOME_SIZES`, `NARROWPEAK`,
`N_SHUFFLES`, `DIST_THRESH_BP`, `BIN_BP`, `COLORS`, etc.) — edit those, then
`python3 proximity_summit_analysis.py`.

### Provenance vs. the published figure (confirmed)

Running this script with the `CONFIG` values as supplied reproduces the
*statistics* behind panel B's enrichment curve (confirmed against
`SMARCA5_summit_enrichment_standalone_EN.pdf`, the standalone enrichment-only
figure this script actually produced: same monotonic decay from ~100-200x
near distance 0 down to ~1.5-2x at 20 kb, DMRT1 enrichment consistently above
RAR across the whole range), but does **not** match the published panel B's
*appearance* out of the box:

- The standalone/original output colors DMRT1 orange and RAR green (as in
  `COLORS` here). The published composite panel recolors them to navy
  (DMRT1) and orange (RAR), matching panel C's palette.
- The standalone output has no dashed reference lines. The published panel
  adds dashed lines at **500 bp and 1 kb**.
- The published panel B is a 3-part composite (schematic + two 0–5 kb summit
  distance histograms + this enrichment curve); this script's direct
  `{REF}_ref_SUMMIT_50bp_EN.png` already *is* a 3-part composite per
  reference factor, but with its own cosmetics (see above) rather than the
  published ones.

In other words: this script's numbers are the real source data for panel B,
but the exact published look (colors, 500 bp/1 kb reference lines, final
layout) was applied in a later manual styling/assembly pass, not by this
script's `CONFIG`.

## What's *not* here

Neither script includes the upstream peak-calling step (MACS2) itself. For
panel C, only the already-split signalValue table was available — not the
original overlap-assignment code — but the rule it used is documented as:

```
bedtools closest -a <factor>.narrowPeak -b <partner>.narrowPeak -d
# last column == 0  -> "overlap" group   (peaks share >= 1 bp)
# last column >  0  -> "no overlap" group
```

If you have the original narrowPeak files, add a small script that runs the
above and writes a two-column signalValue table in the layout
`analyze_signal_overlap.py` expects, and it will consume it unchanged.

## Repo layout

```
scripts/
  analyze_signal_overlap.py        # Panel C: stats + figure (verified against real data)
  proximity_summit_analysis.py     # Panel B: stats + figure (original code; see mismatch note above)
examples/
  make_synthetic_test_data.py      # generates fake narrowPeak files for Panel B's demo
  make_demo_violin_data.py         # generates a fake signalValue workbook for Panel C's demo
requirements.txt
```

## Instructions for use — running on your own data

**Panel C (signal intensity vs. overlap)**

```bash
python scripts/analyze_signal_overlap.py \
    --input violin_plot_underlying_data.xlsx \
    --outdir ./output_C
```

Input: an Excel workbook with one sheet per comparison (`SMARCA5 vs RAR`,
`SMARCA5 vs DMRT1`, `DMRT1 vs RAR`, `DMRT1 vs SMARCA5`), metadata rows 1–3, a
header row at row 4, and signalValue data from row 5 onward in columns A
(group A) and B (group B) — the two columns may have different lengths.

Produces `stats_summary.csv` (U, z, p, and a log10(p)/mantissa/exponent
breakdown for p-values that underflow double precision — see note below) and
`figure_C.png` / `figure_C.pdf`.

**Panel B (summit distance + enrichment)**

```bash
# real run: edit the CONFIG block in proximity_summit_analysis.py to point
# NARROWPEAK / GENOME_SIZES at your real files, then:
python scripts/proximity_summit_analysis.py
```

Inputs: narrowPeak files for SMARCA5, RAR and DMRT1 (column 10 must be a real
summit offset, not -1), and a 2-column `chrom<TAB>size` genome file for
`bedtools shuffle`. The panel B image corresponds to the `SMARCA5_ref_*`
output (SMARCA5 as the reference factor, RAR and DMRT1 as partners); the
script also produces `RAR_ref_*` and `DMRT1_ref_*` as a byproduct of looping
over every reference factor. With `N_SHUFFLES = 20` (default) on real
genome-scale peak sets this can take a while (one `bedtools shuffle` +
`bedtools closest` per shuffle per factor pair).

## Notes on extremely small p-values (Panel C)

For two of the four comparisons, the Mann–Whitney asymptotic p-value
underflows to exactly `0.0` in double precision (SciPy's normal-approximation
p-value computation loses precision once the z-score gets large enough — here
z ≈ -50 and z ≈ -39). `analyze_signal_overlap.py` recovers a finite estimate
by working in log-space with `scipy.stats.norm.logsf`, then reports the
figure annotation as `P < 1e-300` rather than the literal (and statistically
over-precise) mantissa/exponent, since values at that extreme are more a
statement of "the asymptotic approximation says this is enormously
significant" than a trustworthy 3-significant-figure p-value.

## License

Released under the MIT License — see [LICENSE](../LICENSE).
