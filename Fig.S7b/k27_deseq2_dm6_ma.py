#!/usr/bin/env python3
"""
Spike-in (dm6) normalized differential analysis of H3K27me3 CUT&Tag/ChIP signal
at the union of MACS2 narrow peaks, using pyDESeq2, with MA plot and BED export.

Pipeline
--------
1. Union of narrow peaks from the two conditions (bedtools merge).
2. Raw fragment counts per peak per replicate (bedtools intersect -c).
3. DESeq2 size factors derived from spike-in scale factors:
       size_factor = (library size in millions) / scaleFactor
   (bamCoverage convention: signal = RPKM * scaleFactor; DESeq2 convention:
   normalized = raw / size_factor). These are injected directly, bypassing
   DESeq2's median-of-ratios estimation, which assumes most regions are
   unchanged and would cancel a genuine global change in signal.
4. Negative-binomial Wald test (cKO vs Ctrl).
5. MA plot, full results table, and BED files for each category.

Requirements: bedtools on PATH; python packages in requirements.txt.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito colorblind-safe palette
COL_NS, COL_UP, COL_DOWN = "#999999", "#D55E00", "#0072B2"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--peaks-ctrl", required=True, help="MACS2 narrowPeak, Ctrl (merged replicates)")
    p.add_argument("--peaks-cko", required=True, help="MACS2 narrowPeak, cKO (merged replicates)")
    p.add_argument("--beds-ctrl", nargs=2, required=True, metavar=("REP1", "REP2"),
                   help="Per-replicate fragment BED files, Ctrl (bedtools bamtobed output)")
    p.add_argument("--beds-cko", nargs=2, required=True, metavar=("REP1", "REP2"),
                   help="Per-replicate fragment BED files, cKO")
    p.add_argument("--scale-ctrl", type=float, required=True, help="Spike-in scale factor for Ctrl")
    p.add_argument("--scale-cko", type=float, required=True, help="Spike-in scale factor for cKO")
    p.add_argument("--padj", type=float, default=0.05, help="padj threshold (default: 0.05)")
    p.add_argument("--outdir", default="results", help="Output directory (default: results)")
    return p.parse_args()


def run(cmd, stdout=None):
    subprocess.run(cmd, check=True, stdout=stdout)


def count_lines(path):
    with open(path, "rb") as f:
        return sum(1 for _ in f)


def build_union_peaks(peaks_ctrl, peaks_cko, out_bed):
    """Union of both peak sets: cut -f1-3 | sort | bedtools merge."""
    rows = []
    for f in (peaks_ctrl, peaks_cko):
        df = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
        rows.append(df)
    union = pd.concat(rows).sort_values(["chr", "start", "end"])
    tmp = out_bed.with_suffix(".unmerged.bed")
    union.to_csv(tmp, sep="\t", header=False, index=False)
    with open(out_bed, "w") as fh:
        run(["bedtools", "merge", "-i", str(tmp)], stdout=fh)
    tmp.unlink()


def count_fragments(union_bed, frag_bed):
    out = subprocess.run(["bedtools", "intersect", "-c", "-a", str(union_bed), "-b", str(frag_bed)],
                         check=True, capture_output=True, text=True).stdout
    from io import StringIO
    return pd.read_csv(StringIO(out), sep="\t", header=None, names=["chr", "start", "end", "count"])


def main():
    a = parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    samples = ["Ctrl_1", "Ctrl_2", "cKO_1", "cKO_2"]
    beds = dict(zip(samples, a.beds_ctrl + a.beds_cko))
    scale = {"Ctrl_1": a.scale_ctrl, "Ctrl_2": a.scale_ctrl, "cKO_1": a.scale_cko, "cKO_2": a.scale_cko}

    # 1-2. union peaks and raw counts
    union_bed = out / "union_peaks.bed"
    build_union_peaks(a.peaks_ctrl, a.peaks_cko, union_bed)
    cnt = {s: count_fragments(union_bed, beds[s]) for s in samples}
    peaks = cnt[samples[0]][["chr", "start", "end"]]
    for s in samples[1:]:
        assert (cnt[s][["chr", "start", "end"]].values == peaks.values).all(), "peak order mismatch"
    peak_ids = peaks["chr"] + ":" + peaks["start"].astype(str) + "-" + peaks["end"].astype(str)
    counts = pd.DataFrame({s: cnt[s]["count"].values for s in samples}, index=peak_ids)
    counts.to_csv(out / "raw_counts.tsv", sep="\t")
    keep = counts.sum(axis=1) > 0
    counts, peaks, peak_ids = counts[keep], peaks[keep.values], peak_ids[keep.values]
    print(f"Union peaks: {len(keep)}; with >=1 read: {len(counts)}")

    # 3. spike-in derived size factors
    lib_M = {s: count_lines(beds[s]) / 1e6 for s in samples}
    size_factors = np.array([lib_M[s] / scale[s] for s in samples])
    print("Size factors:", dict(zip(samples, size_factors.round(4))))

    # 4. pyDESeq2 with injected size factors
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts_t = counts[samples].T
    meta = pd.DataFrame({"condition": ["Ctrl", "Ctrl", "cKO", "cKO"]}, index=counts_t.index)
    dds = DeseqDataSet(counts=counts_t, metadata=meta, design="~condition")

    # NOTE: written against pydeseq2's AnnData-based internals (size factors in
    # dds.obs, per-peak mean in dds.var["_normed_means"]). These are not public
    # API and may change between versions; see README for the tested version.
    dds.obs["size_factors"] = size_factors
    dds.layers["normed_counts"] = dds.X / dds.obs["size_factors"].values[:, None]
    dds.var["_normed_means"] = dds.layers["normed_counts"].mean(axis=0)
    dds.fit_size_factors = lambda *args, **kwargs: None  # prevent re-estimation

    dds.fit_genewise_dispersions()
    dds.fit_dispersion_trend()
    dds.fit_dispersion_prior()
    dds.fit_MAP_dispersions()
    dds.fit_LFC()
    dds.calculate_cooks()
    dds.refit()
    assert np.allclose(dds.obs["size_factors"].values, size_factors), "size factors were overwritten"

    stat = DeseqStats(dds, contrast=["condition", "cKO", "Ctrl"])
    stat.summary()
    res = stat.results_df.join(peaks.set_axis(peak_ids))
    res = res[["chr", "start", "end", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]].sort_values("padj")
    res.to_csv(out / "deseq2_results.tsv", sep="\t", index=False)

    # 5. categories, BED export, MA plot
    v = res.dropna(subset=["padj", "log2FoldChange", "baseMean"])
    up = (v["padj"] < a.padj) & (v["log2FoldChange"] > 0)
    down = (v["padj"] < a.padj) & (v["log2FoldChange"] < 0)
    ns = ~up & ~down
    for name, mask in [("up_in_cKO", up), ("down_in_cKO", down), ("not_significant", ns)]:
        v[mask].sort_values(["chr", "start"])[["chr", "start", "end"]].to_csv(
            out / f"{name}.bed", sep="\t", header=False, index=False)
        print(f"{name}: {int(mask.sum())}")

    A, M = np.log2(v["baseMean"] + 1), v["log2FoldChange"]
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(A[ns], M[ns], s=6, color=COL_NS, alpha=0.35, label=f"padj>={a.padj} (n={int(ns.sum())})")
    ax.scatter(A[up], M[up], s=6, color=COL_UP, alpha=0.7, label=f"up in cKO (n={int(up.sum())})")
    ax.scatter(A[down], M[down], s=6, color=COL_DOWN, alpha=0.7, label=f"down in cKO (n={int(down.sum())})")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("log2(baseMean + 1), spike-in normalized")
    ax.set_ylabel("log2FoldChange (cKO vs Ctrl)")
    ax.legend(loc="upper right", fontsize=8, markerscale=2)
    plt.tight_layout()
    plt.savefig(out / "ma_plot.pdf")
    plt.savefig(out / "ma_plot.png", dpi=150)
    print(f"Done. Outputs in {out}/")


if __name__ == "__main__":
    sys.exit(main())
