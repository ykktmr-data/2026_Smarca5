#!/usr/bin/env python3
"""
make_demo_data.py

Generates a small SYNTHETIC demo dataset (narrowPeak + fragment BED files) with
the same layout k27_deseq2_dm6_ma.py expects, so the pipeline (union peaks ->
raw counts -> spike-in-normalized DESeq2 -> MA plot/BED export) can be run
end-to-end without any real CUT&Tag data.

This is for demonstrating that the pipeline runs correctly, NOT to reproduce
the numbers in the published figure -- the peaks, counts and scale factors
here are synthetic.

Usage:
    python make_demo_data.py --outdir ./demo_data
"""
import argparse
from pathlib import Path

import numpy as np

N_PEAKS = 60
PEAK_WIDTH = 400
PEAK_SPACING = 2000
CHROM = "chr1"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="./demo_data")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    starts = np.arange(N_PEAKS) * PEAK_SPACING + 10_000
    ends = starts + PEAK_WIDTH

    # Baseline per-peak mean fragment count (same biology in both conditions,
    # except for a handful of up/down peaks injected below).
    base_mean = rng.gamma(shape=3.0, scale=15.0, size=N_PEAKS)
    fc = np.ones(N_PEAKS)
    fc[:8] = rng.uniform(2.5, 4.0, size=8)     # up in cKO
    fc[8:16] = rng.uniform(0.2, 0.4, size=8)   # down in cKO
    rng.shuffle(fc)

    def write_narrowpeak(path, signal):
        with open(path, "w") as fh:
            for i in range(N_PEAKS):
                name = f"peak_{i}"
                summit = PEAK_WIDTH // 2
                fh.write(
                    f"{CHROM}\t{starts[i]}\t{ends[i]}\t{name}\t0\t.\t"
                    f"{signal[i]:.3f}\t-1\t-1\t{summit}\n"
                )

    def write_fragment_bed(path, counts):
        rows = []
        for i in range(N_PEAKS):
            n = max(0, int(counts[i]))
            for _ in range(n):
                frag_start = starts[i] + rng.integers(0, max(1, PEAK_WIDTH - 150))
                frag_end = frag_start + 150
                rows.append((CHROM, frag_start, frag_end))
        rows.sort()
        with open(path, "w") as fh:
            for chrom, s, e in rows:
                fh.write(f"{chrom}\t{s}\t{e}\n")
        return len(rows)

    # Ctrl replicates: baseline counts with Poisson/replicate noise.
    ctrl1 = rng.poisson(base_mean * 1.0)
    ctrl2 = rng.poisson(base_mean * 1.05)
    # cKO replicates: baseline * fold-change, with its own noise.
    cko1 = rng.poisson(base_mean * fc * 0.95)
    cko2 = rng.poisson(base_mean * fc * 1.0)

    write_narrowpeak(outdir / "Ctrl_merged_peaks.narrowPeak", (ctrl1 + ctrl2) / 2)
    write_narrowpeak(outdir / "cKO_merged_peaks.narrowPeak", (cko1 + cko2) / 2)

    n_ctrl1 = write_fragment_bed(outdir / "Ctrl_1.rmDup.bed", ctrl1)
    n_ctrl2 = write_fragment_bed(outdir / "Ctrl_2.rmDup.bed", ctrl2)
    n_cko1 = write_fragment_bed(outdir / "cKO_1.rmDup.bed", cko1)
    n_cko2 = write_fragment_bed(outdir / "cKO_2.rmDup.bed", cko2)

    print(f"Wrote {N_PEAKS} synthetic peaks and fragment BEDs to {outdir}/")
    print(f"  Ctrl_1: {n_ctrl1} fragments, Ctrl_2: {n_ctrl2} fragments")
    print(f"  cKO_1:  {n_cko1} fragments, cKO_2:  {n_cko2} fragments")
    print("\nDemo scale factors (arbitrary, for pipeline smoke-testing only):")
    print("  --scale-ctrl 1.0 --scale-cko 1.0")


if __name__ == "__main__":
    main()
