#!/usr/bin/env python3
"""
analyze_signal_overlap.py

Statistical analysis and figure generation for the "signal intensity vs. binding-site
overlap" violin plots (Figure panel C: SMARCA5/DMRT1 ChIP peak signalValue, grouped by
whether each peak overlaps an RAR / DMRT1 / SMARCA5 peak).

INPUT
-----
An Excel workbook with one sheet per comparison, each sheet holding two columns of
peak signalValue (one column per group), e.g.:

    SMARCA5 vs RAR        : col A = "no RAR overlap"      col B = "overlaps RAR"
    SMARCA5 vs DMRT1       : col A = "no DMRT1 overlap"    col B = "overlaps DMRT1"
    DMRT1 vs RAR          : col A = "no RAR overlap"      col B = "overlaps RAR"
    DMRT1 vs SMARCA5      : col A = "no SMARCA5 overlap"  col B = "overlaps SMARCA5"

Group membership (upstream of this script, NOT reproduced here) was defined as:

    "Overlap" = the two peak regions share >= 1 bp, i.e. `bedtools closest -d`
                reports a distance of 0 between a peak and its nearest partner-factor
                peak. Distance > 0 => "no overlap".

    This script does not call bedtools and does not have access to the original peak
    (BED/narrowPeak) files -- it only consumes the already-split signalValue numbers.
    If you have the original peak files, the overlap step can be reproduced with:

        bedtools closest -a <factor>.narrowPeak -b <partner>.narrowPeak -d \
            > <factor>_vs_<partner>.closest.bed
        # column (last) == 0  -> "overlap" group
        # column (last) >  0  -> "no overlap" group

WHAT THIS SCRIPT DOES
----------------------
1. Loads the signalValue data for each of the 4 comparisons.
2. Runs a two-sided Mann-Whitney U test for each comparison (SciPy, asymptotic/
   normal approximation with tie correction).
3. For p-values that underflow to 0.0 in double precision, recovers a non-zero
   estimate via the log-survival function (mantissa + base-10 exponent), so
   extremely significant results are not silently reported as "n.s."
4. Reproduces the violin + box + jittered-scatter figure (Figure panel C) with the
   corrected p-value / "n.s." annotations.

USAGE
-----
    pip install -r ../requirements.txt
    python analyze_signal_overlap.py --input violin_plot_underlying_data.xlsx \
        --outdir ./output

Outputs (written to --outdir):
    stats_summary.csv          - U, z, p (and log10(p) for underflowed cases) per comparison
    figure_C.png / figure_C.pdf
"""

import argparse
import math
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# Config: one entry per violin-plot panel, in left-to-right display order.
# ----------------------------------------------------------------------------
GREEN, ORANGE, NAVY = "#43974E", "#D95F2E", "#1E2083"

PANELS = [
    dict(
        sheet="SMARCA5 vs RAR",
        ylabel="SMARCA5 signal intensity",
        group_labels=("no RAR overlap", "RAR overlap"),
        colors=(GREEN, ORANGE),
    ),
    dict(
        sheet="SMARCA5 vs DMRT1",
        ylabel="SMARCA5 signal intensity",
        group_labels=("no DMRT1 overlap", "DMRT1 overlap"),
        colors=(GREEN, NAVY),
    ),
    dict(
        sheet="DMRT1 vs RAR",
        ylabel="DMRT1 signal intensity",
        group_labels=("no RAR overlap", "RAR overlap"),
        colors=(NAVY, ORANGE),
    ),
    dict(
        sheet="DMRT1 vs SMARCA5",
        ylabel="DMRT1 signal intensity",
        group_labels=("no SMARCA5 overlap", "SMARCA5 overlap"),
        colors=(NAVY, GREEN),
    ),
]

# p-values more extreme than this are reported as "< EXTREME_P_CUTOFF" rather than
# with their literal (statistically over-precise) mantissa/exponent.
EXTREME_P_CUTOFF = 1e-300


def load_sheet_pair(xlsx_path: str, sheet: str) -> tuple[np.ndarray, np.ndarray]:
    """Load the two signalValue columns of one comparison sheet.

    Skips the metadata/description rows at the top of each sheet (header lives at
    row 4; data starts row 5), and drops any blank cells so the two groups can have
    different lengths.
    """
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb[sheet]
    col_a, col_b = [], []
    for row in ws.iter_rows(min_row=5, values_only=True):
        a, b = row[0], row[1]
        if a is not None:
            col_a.append(a)
        if b is not None:
            col_b.append(b)
    return np.array(col_a, dtype=float), np.array(col_b, dtype=float)


def mannwhitney_with_underflow_guard(a: np.ndarray, b: np.ndarray) -> dict:
    """Two-sided Mann-Whitney U test, with a tie-corrected normal-approximation
    fallback that survives p-values too small for float64 (which SciPy reports as
    an unhelpful literal 0.0).

    Returns a dict with U, the SciPy p-value, and (when the SciPy p-value underflows)
    a (mantissa, exponent) pair such that p ~= mantissa * 10**exponent.
    """
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic")

    n1, n2 = len(a), len(b)
    mu = n1 * n2 / 2
    all_vals = np.concatenate([a, b])
    _, counts = np.unique(all_vals, return_counts=True)
    tie_term = np.sum(counts**3 - counts)
    n_total = n1 + n2
    sigma = math.sqrt(n1 * n2 / 12 * ((n_total + 1) - tie_term / (n_total * (n_total - 1))))
    z = (U - mu) / sigma

    log10_p = (math.log(2) + stats.norm.logsf(abs(z))) / math.log(10)
    mantissa = 10 ** (log10_p - math.floor(log10_p))
    exponent = math.floor(log10_p)

    return dict(U=U, p_scipy=p, z=z, log10_p=log10_p, mantissa=mantissa, exponent=exponent)


def format_p_annotation(result: dict) -> str:
    """Render a p-value for the figure, falling back to '< 10^-300' style notation
    when the exact value underflows double precision."""
    p = result["p_scipy"]
    if p > 0:
        return f"P = {p:.2e}".replace("e-0", "e-").replace("e-", " x 10$^{-") + "}$"
    cutoff_exp = int(math.log10(EXTREME_P_CUTOFF))
    return f"P < 1 x 10$^{{{cutoff_exp}}}$"


def run_stats(xlsx_path: str) -> tuple[dict, pd.DataFrame]:
    data = {}
    rows = []
    for panel in PANELS:
        a, b = load_sheet_pair(xlsx_path, panel["sheet"])
        data[panel["sheet"]] = (a, b)
        res = mannwhitney_with_underflow_guard(a, b)
        rows.append(
            dict(
                comparison=panel["sheet"],
                group_A=panel["group_labels"][0],
                n_A=len(a),
                median_A=np.median(a),
                group_B=panel["group_labels"][1],
                n_B=len(b),
                median_B=np.median(b),
                U=res["U"],
                z=res["z"],
                p_scipy=res["p_scipy"],
                p_mantissa=res["mantissa"],
                p_exponent=res["exponent"],
                p_annotation=format_p_annotation(res),
            )
        )
    return data, pd.DataFrame(rows)


def make_figure(data: dict, summary: pd.DataFrame, out_png: Path, out_pdf: Path, seed: int = 0):
    rng = np.random.default_rng(seed)
    fig, axes = plt.subplots(1, len(PANELS), figsize=(16.5, 5.6))

    for ax, panel, (_, row) in zip(axes, PANELS, summary.iterrows()):
        a, b = data[panel["sheet"]]
        colors = panel["colors"]
        y_max = max(a.max(), b.max())
        y_min = min(a.min(), b.min())

        for pos, group, color in zip((1, 2), (a, b), colors):
            x = rng.normal(loc=pos, scale=0.045, size=len(group))
            ax.scatter(x, group, s=3, color=color, alpha=0.55, linewidths=0, zorder=2)

            vp = ax.violinplot([group], positions=[pos], widths=0.7, showextrema=False)
            for body in vp["bodies"]:
                body.set_facecolor("none")
                body.set_edgecolor(color)
                body.set_linewidth(1.3)
                body.set_zorder(3)

            ax.boxplot(
                [group], positions=[pos], widths=0.18, patch_artist=True, showfliers=False,
                zorder=4,
                medianprops=dict(color=color, linewidth=1.2),
                boxprops=dict(facecolor="white", edgecolor=color, linewidth=1.1),
                whiskerprops=dict(color=color, linewidth=1.1),
                capprops=dict(color=color, linewidth=1.1),
            )
            ax.scatter([pos], [group.mean()], s=18, facecolors="white", edgecolors=color,
                       linewidths=1.2, zorder=5)

        y0, h = y_max * 1.04, y_max * 0.03
        ax.plot([1, 1, 2, 2], [y0, y0 + h, y0 + h, y0], color="black", linewidth=1.0)
        ax.text(1.5, y0 + h * 1.3, row["p_annotation"], ha="center", va="bottom", fontsize=11)

        ax.set_xlim(0.5, 2.5)
        ax.set_ylim(min(0, y_min - 1), y_max * 1.22)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(
            [f"{lab}\n(n={n:,})" for lab, n in zip(panel["group_labels"], (row["n_A"], row["n_B"]))],
            fontsize=10,
        )
        ax.set_ylabel(panel["ylabel"], fontsize=11)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color="#e6e6e6", linewidth=0.6, zorder=0)
        ax.set_axisbelow(True)

    axes[0].text(-0.28, 1.08, "C", transform=axes[0].transAxes, fontsize=20,
                 fontweight="bold", va="top")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="Path to violin_plot_underlying_data.xlsx")
    ap.add_argument("--outdir", default="./output", help="Output directory")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    data, summary = run_stats(args.input)
    summary.to_csv(outdir / "stats_summary.csv", index=False)
    print(summary.to_string(index=False))

    make_figure(data, summary, outdir / "figure_C.png", outdir / "figure_C.pdf")
    print(f"\nWrote {outdir/'stats_summary.csv'}, {outdir/'figure_C.png'}, {outdir/'figure_C.pdf'}")


if __name__ == "__main__":
    main()
