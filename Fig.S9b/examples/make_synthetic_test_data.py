#!/usr/bin/env python3
"""
make_synthetic_test_data.py

Generates small, fake narrowPeak files + a chrom.sizes file so that
proximity_summit_analysis.py can be sanity-checked end-to-end before pointing
it at real SMARCA5 / RAR / DMRT1 peak calls.

Filenames match the CONFIG block in proximity_summit_analysis.py exactly, so
you can run this generator, copy its output next to the script (or point
GENOME_SIZES / NARROWPEAK at this directory), and do a dry run with no edits.

A fraction of the RAR and DMRT1 peaks are deliberately placed close to a
random SMARCA5 summit (controlled by --cluster-frac-*) so the resulting
histograms/enrichment curve show the expected "excess near distance 0" shape
instead of pure noise. This is test fixture data only -- it proves the code
runs and produces sane shapes, not that it matches the real biology.

Usage:
    python make_synthetic_test_data.py --outdir ./test_data
    cd ./test_data && python ../../scripts/proximity_summit_analysis.py
"""
import argparse
from pathlib import Path

import numpy as np

# Must match the NARROWPEAK dict / GENOME_SIZES in proximity_summit_analysis.py
SMARCA5_FILENAME = "CT_Snf2h_P7_rep_merge_bowtie2_1e4_peaks.narrowPeak"
RAR_FILENAME = "trimmed_GS_RAR_RARminus_SRR7500346_SRR7500347_SRR7500348_1e4_peaks.narrowPeak"
DMRT1_FILENAME = "SV_CT_P8U_Dmrt1_C_rep_merge_bowtie2.sorted.rmDup.bam_1e5_peaks.narrowPeak"
GENOME_SIZES_FILENAME = "mm10.chrom.sizes"


def make_narrowpeak(path, n, chrom_sizes, rng, cluster_near=None, cluster_frac=0.0, cluster_sd=300):
    rows = []
    chroms = list(chrom_sizes.keys())
    for i in range(n):
        chrom = rng.choice(chroms)
        size = chrom_sizes[chrom]
        if cluster_near is not None and rng.random() < cluster_frac:
            c2, pos2 = cluster_near[rng.integers(len(cluster_near))]
            chrom = c2
            start = int(np.clip(pos2 + rng.normal(0, cluster_sd), 0, chrom_sizes[chrom] - 200))
        else:
            start = rng.integers(0, size - 200)
        end = start + 150
        summit_offset = 75
        rows.append(
            f"{chrom}\t{start}\t{end}\tpeak{i}\t{rng.integers(100, 1000)}\t.\t"
            f"{rng.uniform(2, 10):.4f}\t0\t0\t{summit_offset}"
        )
    with open(path, "w") as f:
        f.write("\n".join(rows) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="./test_data")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-smarca5", type=int, default=2000)
    ap.add_argument("--n-rar", type=int, default=500)
    ap.add_argument("--n-dmrt1", type=int, default=900)
    ap.add_argument("--cluster-frac-rar", type=float, default=0.3)
    ap.add_argument("--cluster-frac-dmrt1", type=float, default=0.5)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # a couple of "non-standard" chroms included on purpose, to exercise
    # proximity_summit_analysis.py's STANDARD_CHROM_RE filtering step
    chrom_sizes = {"chr1": 2_000_000, "chr2": 1_500_000, "chrM": 16_000, "chr1_random": 50_000}
    with open(outdir / GENOME_SIZES_FILENAME, "w") as f:
        for c, s in chrom_sizes.items():
            f.write(f"{c}\t{s}\n")

    make_narrowpeak(outdir / SMARCA5_FILENAME, args.n_smarca5, chrom_sizes, rng)

    smarca5_summits = []
    with open(outdir / SMARCA5_FILENAME) as f:
        for line in f:
            p = line.strip().split("\t")
            smarca5_summits.append((p[0], int(p[1]) + int(p[9])))

    make_narrowpeak(outdir / RAR_FILENAME, args.n_rar, chrom_sizes, rng,
                     cluster_near=smarca5_summits, cluster_frac=args.cluster_frac_rar, cluster_sd=400)
    make_narrowpeak(outdir / DMRT1_FILENAME, args.n_dmrt1, chrom_sizes, rng,
                     cluster_near=smarca5_summits, cluster_frac=args.cluster_frac_dmrt1, cluster_sd=300)

    print(f"Wrote synthetic test data to {outdir}/")
    print("Try:")
    print(f"  cd {outdir} && python ../../scripts/proximity_summit_analysis.py")


if __name__ == "__main__":
    main()
