# 2026_Smarca5

Analysis code for selected supplementary figure panels in the SMARCA5
spermatogenesis manuscript (Kitamura et al., submitted to *Nature
Communications*).

## Contents

| Folder | Figure panel(s) | What it does |
|---|---|---|
| [`Fig.S9b/`](./Fig.S9b) | Supplementary Fig. S9b, panels B & C | SMARCA5/DMRT1/RAR peak co-localization: signal-intensity-vs-overlap statistics (panel C) and summit-distance/enrichment analysis (panel B) |
| [`Fig.S9c/`](./Fig.S9c) | Supplementary Fig. S9c | Standalone copy of the panel C (signal-intensity-vs-overlap) analysis above |
| [`Fig.S7b/`](./Fig.S7b) | Supplementary Fig. S7b | Spike-in normalized differential H3K27me3 analysis (pyDESeq2) between Ctrl and cKO, with MA plot |

Each folder is self-contained: its own README (with system requirements,
installation guide, a runnable demo on synthetic data, and instructions for
running on real data), `requirements.txt`, and scripts.

## Sequencing data

Raw and processed RNA-seq, ATAC-seq and CUT&Tag data are deposited in the
Gene Expression Omnibus under accession
[GSE303063](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303063).

## License

Released under the MIT License — see [LICENSE](./LICENSE).
