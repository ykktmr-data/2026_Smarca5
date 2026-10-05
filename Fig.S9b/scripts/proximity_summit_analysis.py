#!/usr/bin/env python3
"""
proximity_summit_analysis.py

Genomic proximity analysis between RAR, SMARCA5, and DMRT1 ChIP-seq / CUT&Tag
peak summits.

For each pair of factors (A, B), this script:
  1. Extracts the peak summit (single base pair of maximum signal) from each
     MACS2 narrowPeak file, using the built-in summit-offset column (column 10:
     offset of the summit from the peak start).
  2. For every peak of factor A, finds the nearest peak summit of factor B
     using `bedtools closest` (signed distance is not used; distance is the
     number of bp between the two nearest points, 0 if identical).
  3. Builds a null / background distribution by shuffling factor B's summits
     within their chromosome of origin (`bedtools shuffle -chrom`) N times,
     and repeating the nearest-neighbour search against each shuffled set.
  4. Reports, for a grid of distance thresholds, the fraction of factor-A
     peaks within that threshold of a factor-B peak, the same fraction
     expected under the randomised background, and their ratio
     (observed / random = fold-enrichment).
  5. Renders a summary figure per reference factor: 50-bp-binned distance
     histograms to each partner factor (0-5 kb) plus the enrichment-vs-random
     curve out to 20 kb.

Requirements
------------
  * bedtools >= 2.30 on PATH
  * Python packages: pandas, numpy, matplotlib, openpyxl

Usage
-----
Edit the CONFIG block below to point at your narrowPeak files and a
chromosome-sizes file for your genome build (used by `bedtools shuffle`),
then run:

    python3 proximity_summit_analysis.py

All intermediate BED/TSV files and the final tables/figures are written to
`OUT_DIR`.

Notes on peak-calling thresholds used in this study
----------------------------------------------------
RAR and SMARCA5 peaks were called with MACS2 at p < 1e-4; DMRT1 peaks were
called with MACS2 at p < 1e-5 (as encoded in the source narrowPeak file
names). Genome build: mm10 (confirmed by comparing the maximum observed
peak coordinate per chromosome to the reference mm10 chromosome lengths).
"""

import os
import subprocess
import itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ----------------------------------------------------------------------
# CONFIG - edit these paths for your own data
# ----------------------------------------------------------------------
OUT_DIR = "proximity_out"
GENOME_SIZES = "mm10.chrom.sizes"     # 2-column TSV: chrom<TAB>size
STANDARD_CHROM_RE = r"^chr([0-9]+|X|Y)$"   # excludes chrM and unplaced scaffolds

NARROWPEAK = {
    "RAR":     "trimmed_GS_RAR_RARminus_SRR7500346_SRR7500347_SRR7500348_1e4_peaks.narrowPeak",
    "SMARCA5": "CT_Snf2h_P7_rep_merge_bowtie2_1e4_peaks.narrowPeak",
    "DMRT1":   "SV_CT_P8U_Dmrt1_C_rep_merge_bowtie2.sorted.rmDup.bam_1e5_peaks.narrowPeak",
}

N_SHUFFLES = 20
DIST_THRESH_BP = 5000          # "nearby" cutoff highlighted in figures/tables
ZOOM_MAX_BP = 5000              # histogram x-axis range
BIN_BP = 50                     # histogram bin width
ENRICHMENT_MAX_BP = 20000       # x-axis range of the enrichment curve

COLORS = {"RAR": "#1baf7a", "SMARCA5": "#2a78d6", "DMRT1": "#eb6834"}

os.makedirs(OUT_DIR, exist_ok=True)


# ----------------------------------------------------------------------
# Step 1: summit extraction
# ----------------------------------------------------------------------
def extract_summits(narrowpeak_path, out_bed):
    """Write a 1-bp BED of peak summits (start + column-10 offset)."""
    with open(narrowpeak_path) as fin, open(out_bed, "w") as fout:
        for line in fin:
            f = line.rstrip("\r\n").split("\t")
            chrom, start, name, offset = f[0], int(f[1]), f[3], int(f[9])
            summit = start + offset
            fout.write(f"{chrom}\t{summit}\t{summit + 1}\t{name}\n")
    subprocess.run(f"sort -k1,1 -k2,2n -o {out_bed} {out_bed}", shell=True, check=True)


def filter_standard_chroms(in_bed, out_bed):
    subprocess.run(
        f"awk 'BEGIN{{OFS=\"\\t\"}} $1 ~ /{STANDARD_CHROM_RE}/' {in_bed} "
        f"| sort -k1,1 -k2,2n > {out_bed}",
        shell=True, check=True,
    )


summit_bed = {}
for factor, path in NARROWPEAK.items():
    raw = os.path.join(OUT_DIR, f"{factor}.summit.raw.bed")
    clean = os.path.join(OUT_DIR, f"{factor}.summit.bed")
    extract_summits(path, raw)
    filter_standard_chroms(raw, clean)
    summit_bed[factor] = clean

# ----------------------------------------------------------------------
# Step 2: pairwise nearest-summit distances (bedtools closest)
# ----------------------------------------------------------------------
def bedtools_closest(a_bed, b_bed, out_tsv):
    subprocess.run(
        f"bedtools closest -a {a_bed} -b {b_bed} -d -t first "
        f"| awk -F'\\t' '$NF != -1' > {out_tsv}",   # drop peaks with no partner on their chrom
        shell=True, check=True,
    )


factors = list(NARROWPEAK.keys())
closest_tsv = {}
for a, b in itertools.permutations(factors, 2):
    out = os.path.join(OUT_DIR, f"CLOSEST_{a}_vs_{b}.tsv")
    bedtools_closest(summit_bed[a], summit_bed[b], out)
    closest_tsv[(a, b)] = out


def load_distance(tsv_path):
    return pd.read_csv(tsv_path, sep="\t", header=None).iloc[:, -1].astype(float).values


# ----------------------------------------------------------------------
# Step 3: chromosome-restricted random shuffles + null distributions
# ----------------------------------------------------------------------
def bedtools_shuffle(in_bed, out_bed, seed):
    subprocess.run(
        f"bedtools shuffle -i {in_bed} -g {GENOME_SIZES} -chrom -seed {seed} "
        f"| sort -k1,1 -k2,2n > {out_bed}",
        shell=True, check=True,
    )


shuffle_dir = os.path.join(OUT_DIR, "shuffle")
os.makedirs(shuffle_dir, exist_ok=True)

shuffled_bed = {f: [] for f in factors}
for f in factors:
    for i in range(1, N_SHUFFLES + 1):
        out = os.path.join(shuffle_dir, f"{f}_shuf_{i}.bed")
        bedtools_shuffle(summit_bed[f], out, seed=i)
        shuffled_bed[f].append(out)

# null distance pool: for each (reference, target) pair, concatenate the
# reference-vs-shuffled-target distances across all N shuffles
null_distances = {}
for a, b in itertools.permutations(factors, 2):
    pooled = []
    for shuf in shuffled_bed[b]:
        tmp = os.path.join(shuffle_dir, "tmp_closest.tsv")
        bedtools_closest(summit_bed[a], shuf, tmp)
        pooled.append(load_distance(tmp))
    null_distances[(a, b)] = np.concatenate(pooled)

# ----------------------------------------------------------------------
# Step 4: cumulative observed-vs-random enrichment tables
# ----------------------------------------------------------------------
def enrichment_table(obs, null, edges):
    rows = []
    for t in edges:
        o = 100 * (obs <= t).sum() / len(obs)
        c = 100 * (null <= t).sum() / len(null)
        rows.append({"distance_bp": t, "observed_pct": o, "random_pct": c,
                      "enrichment": (o / c) if c > 0 else np.nan})
    return pd.DataFrame(rows)


edges = list(range(0, DIST_THRESH_BP + BIN_BP, BIN_BP)) + list(
    range(6000, ENRICHMENT_MAX_BP + 1000, 1000)
)

summary_rows = []
enrichment_tables = {}
for a, b in itertools.permutations(factors, 2):
    obs = load_distance(closest_tsv[(a, b)])
    null = null_distances[(a, b)]
    et = enrichment_table(obs, null, edges)
    enrichment_tables[(a, b)] = et
    at_thresh = et.loc[et.distance_bp == DIST_THRESH_BP].iloc[0]
    summary_rows.append({
        "reference": a, "target": b, "n_peaks": len(obs),
        "median_distance_bp": float(np.median(obs)),
        f"pct_within_{DIST_THRESH_BP}bp": round(100 * (obs <= DIST_THRESH_BP).sum() / len(obs), 2),
        f"enrichment_at_{DIST_THRESH_BP}bp": round(at_thresh["enrichment"], 2),
    })

summary = pd.DataFrame(summary_rows)
summary.to_csv(os.path.join(OUT_DIR, "summary_all_pairs.csv"), index=False)
print(summary.to_string(index=False))

with pd.ExcelWriter(os.path.join(OUT_DIR, "proximity_summit_analysis.xlsx")) as xl:
    summary.to_excel(xl, sheet_name="summary", index=False)
    for (a, b), et in enrichment_tables.items():
        et.to_excel(xl, sheet_name=f"{a}_vs_{b}"[:31], index=False)

# ----------------------------------------------------------------------
# Step 5: figures - one per reference factor
# ----------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.edgecolor": "#e3e2dd", "axes.labelcolor": "#0b0b0b",
    "text.color": "#0b0b0b", "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.grid": True, "grid.color": "#e3e2dd", "grid.linewidth": 0.8,
    "figure.facecolor": "#ffffff", "axes.facecolor": "#ffffff",
    "savefig.facecolor": "#ffffff", "pdf.fonttype": 42,
})


def make_figure(ref, targets, out_prefix):
    fig = plt.figure(figsize=(12, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.45, wspace=0.38)
    ax_hists = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    ax_enrich = fig.add_subplot(gs[1, :])

    for target, ax in zip(targets, ax_hists):
        d = load_distance(closest_tsv[(ref, target)])
        dz = d[d <= ZOOM_MAX_BP]
        bins = np.arange(0, ZOOM_MAX_BP + BIN_BP, BIN_BP)
        ax.hist(dz, bins=bins, color=COLORS[target], alpha=0.85,
                edgecolor="white", linewidth=0.15)
        ax.axvline(DIST_THRESH_BP, color="#0b0b0b", linestyle="--", linewidth=1.2)
        pct = 100 * (d <= DIST_THRESH_BP).sum() / len(d)
        ax.set_title(f"{ref} $\\rightarrow$ {target} summit-to-summit distance "
                      f"(0-{ZOOM_MAX_BP//1000} kb, {BIN_BP} bp bins)\n"
                      f"{pct:.1f}% within {DIST_THRESH_BP//1000} kb  (n={len(d):,})",
                      fontsize=10)
        ax.set_xlabel("distance (bp)")
        ax.set_ylabel(f"{ref} peak count")
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x/1000:g}kb" if x > 0 else "0"))
        ax.set_xticks(np.arange(0, ZOOM_MAX_BP + 500, 500))
        ax.spines[["top", "right"]].set_visible(False)

    for target in targets:
        et = enrichment_tables[(ref, target)]
        sub = et[et.distance_bp <= ENRICHMENT_MAX_BP]
        ax_enrich.plot(sub.distance_bp, sub.enrichment, color=COLORS[target],
                        linewidth=2, label=target)
    ax_enrich.axhline(1, color="#52514e", linestyle=":", linewidth=1.2)
    ax_enrich.axvline(DIST_THRESH_BP, color="#0b0b0b", linestyle="--", linewidth=1.2)
    ax_enrich.set_yscale("log")
    ax_enrich.set_xlabel("distance threshold (bp) - cumulative ($\\leq$ x)")
    ax_enrich.set_ylabel("enrichment (observed / random shuffle)")
    ax_enrich.set_title(
        f"Enrichment over chromosome-restricted random shuffle of summit points "
        f"(n={N_SHUFFLES} shuffles, averaged)\n"
        f"dashed = {DIST_THRESH_BP//1000} kb threshold, dotted = no-enrichment line",
        fontsize=11)
    ax_enrich.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{int(x/1000)}kb" if x > 0 else "0"))
    ax_enrich.legend(frameon=False, loc="upper right")
    ax_enrich.spines[["top", "right"]].set_visible(False)

    fig.suptitle(f"Proximity of {ref} peak summits to nearest "
                  f"{' / '.join(targets)} peak summits", fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(os.path.join(OUT_DIR, f"{out_prefix}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT_DIR, f"{out_prefix}.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)


factor_pairs = {
    "RAR": ["SMARCA5", "DMRT1"],
    "SMARCA5": ["RAR", "DMRT1"],
    "DMRT1": ["SMARCA5", "RAR"],
}
for ref, targets in factor_pairs.items():
    make_figure(ref, targets, f"{ref}_ref_SUMMIT_50bp_EN")

print("Done. See", OUT_DIR)
