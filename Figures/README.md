# Figures — Source Code

Code used to generate each main and supplementary figure panel in
*"MIA-Induced Metabolic Priming Increases Vulnerability to Neonatal Hypoxia"*.

Each panel folder contains the analysis/plotting script and its rendered output:

- `*.Rmd` — R Markdown source (R panels)
- `*.ipynb` — Jupyter notebook source (Python panels)
- `*.html` — rendered notebook with code, output, and the final figure inline

> **Data not included.** Scripts read from a local `data/` directory that is not
> distributed in this repository (see *Data availability* in the manuscript and the
> top-level [README](../README.md)). The committed `.html` files show every step and
> result inline, so the analysis is fully reviewable without re-running the code.

Some spatial-transcriptomics panels are produced from upstream analyses kept in the
top-level method directories — these are cross-referenced below rather than duplicated:
[`../CoGAPS`](../CoGAPS), [`../LIANA`](../LIANA), [`../Tricycle`](../Tricycle).

---

## Main figures

### Figure 1 — Long-term cerebellar motor deficits
| Panel | File | Description |
|---|---|---|
| 1D–E | [`Fig1/Fig1D-E`](Fig1/Fig1D-E) | ErasmusLadder cerebellar motor-learning analysis |

### Figure 2 — Cerebellar cytoarchitecture (mouse & human histomorphometry)
| Panel | File | Description |
|---|---|---|
| 2B | [`Fig2/Fig2B`](Fig2/Fig2B) | P11 Purkinje-cell count (cells/mm²) |
| 2C | [`Fig2/Fig2C`](Fig2/Fig2C) | P45 Purkinje-cell count (cells/mm²) |
| 2D | [`Fig2/Fig2D`](Fig2/Fig2D) | P11 EGL surface area (µm²) |
| 2E | [`Fig2/Fig2E`](Fig2/Fig2E) | P11 granule-cell count per HPF (IGL) |
| 2F | [`Fig2/Fig2F`](Fig2/Fig2F) | P45 granule-cell count per HPF (IGL) |
| 2H | [`Fig2/Fig2H`](Fig2/Fig2H) | Human EGL surface area (Term vs Preterm) |
| 2I | [`Fig2/Fig2I`](Fig2/Fig2I) | Human GCP count per HPF (EGL), Term vs Preterm |
| 2J | [`Fig2/Fig2J`](Fig2/Fig2J) | Human IGL surface area (Term vs Preterm) |
| 2K | [`Fig2/Fig2K`](Fig2/Fig2K) | Human GCP count per HPF (IGL), Term vs Preterm |

### Figure 3 — Transcriptome reprogramming
| Panel | File | Description |
|---|---|---|
| 3B | [`Fig3/Fig3B`](Fig3/Fig3B) | CoGAPS patterns vs treatment (Spearman correlation) — see [`../CoGAPS`](../CoGAPS) |

### Figure 4 — Metabolic priming & mitochondrial bioenergetics
| Panel | File | Description |
|---|---|---|
| 4A | [`Fig4/Fig4A`](Fig4/Fig4A) | Mitochondrial proton leak at P11 |
| 4B | [`Fig4/Fig4B`](Fig4/Fig4B) | Mitochondrial proton leak at P45 |
| 4C | [`Fig4/Fig4C`](Fig4/Fig4C) | Seahorse XF OCR trace (P45: Sal-N vs LPS-Hx) |
| 4E | [`Fig4/Fig4E`](Fig4/Fig4E) | P11 granule-cell mitochondrial density |
| 4F | [`Fig4/Fig4F`](Fig4/Fig4F) | P45 granule-cell mitochondrial density |
| 4G | [`Fig4/Fig4G`](Fig4/Fig4G) | P11 granule-cell mitochondrial surface area |
| 4H | [`Fig4/Fig4H`](Fig4/Fig4H) | P45 granule-cell mitochondrial surface area |
| 4I | [`Fig4/Fig4I`](Fig4/Fig4I) | P11 Purkinje-cell mitochondrial density |
| 4J | [`Fig4/Fig4J`](Fig4/Fig4J) | P45 Purkinje-cell mitochondrial density |
| 4K | [`Fig4/Fig4K`](Fig4/Fig4K) | P11 Purkinje-cell mitochondrial surface area |
| 4L | [`Fig4/Fig4L`](Fig4/Fig4L) | P45 Purkinje-cell mitochondrial surface area |
| 4N | [`Fig4/Fig4N`](Fig4/Fig4N) | Mitochondrial-health pie charts (granule cells) |
| 4O | [`Fig4/Fig4O`](Fig4/Fig4O) | Mitochondrial-health pie charts (Purkinje cells) |

### Figure 5 — Neuronal cell-cycle dynamics
| Panel | File | Description |
|---|---|---|
| 5A | [`Fig5/Fig5A`](Fig5/Fig5A) | Spatial Tricycle theta on tissue — see [`../Tricycle`](../Tricycle) |
| 5B | [`Fig5/Fig5B`](Fig5/Fig5B) | Schwabe cell-cycle stage proportions (GC*/PC* × treatment) |
| 5C | [`Fig5/Fig5C`](Fig5/Fig5C) | p21 immunoreactivity at P11 |
| 5D | [`Fig5/Fig5D`](Fig5/Fig5D) | p27 immunoreactivity at P11 |

### Figure 6 — Synaptic puncta & Purkinje dendritic morphometry
| Panel | File | Description |
|---|---|---|
| 6B | [`Fig6/Fig6B`](Fig6/Fig6B) | VGLUT2⁺ puncta density (total ML); incl. `R/build_data.R` |
| 6C | [`Fig6/Fig6C`](Fig6/Fig6C) | VGLUT2⁺ laminar distribution (Outer/Inner ML ratio); incl. `R/build_data.R` |
| 6E | [`Fig6/Fig6E`](Fig6/Fig6E) | Purkinje-cell dendritic length (P45) |
| 6F | [`Fig6/Fig6F`](Fig6/Fig6F) | Purkinje-cell dendritic complexity (P45) |

### Figure 7 — Social, exploratory & kinematic behavior
| Panel | File | Description |
|---|---|---|
| 7B | [`Fig7/Fig7B`](Fig7/Fig7B) | Three-chamber Social Interaction Test (SIT) |
| 7C | [`Fig7/Fig7C`](Fig7/Fig7C) | Open-field trajectory occupancy (KDE); incl. `lib/oft_movement_analysis.py` |
| 7D | [`Fig7/Fig7D`](Fig7/Fig7D) | Open-field Test (OFT) |
| 7E | [`Fig7/Fig7E`](Fig7/Fig7E) | B-SOiD 3D-UMAP classified behaviors |
| 7F | [`Fig7/Fig7F`](Fig7/Fig7F) | B-SOiD three-chamber social-behavior clustermap |

---

## Supplementary figures

| Figure | File | Description |
|---|---|---|
| S1 | [`FigS1`](FigS1) | Cross-species enrichment (CAMERA) of P11 MIA+Hx insult-unique signature in human preterm cerebellum |
| S3 | [`FigS3`](FigS3) | CoGAPS patterns vs cerebellar cell types (Spearman) — see [`../CoGAPS`](../CoGAPS) |
| S4 | [`FigS4`](FigS4) | GSEA GO enrichment of CoGAPS patterns M-1 and M-4 — see [`../CoGAPS`](../CoGAPS) |
| S8 | [`FigS8`](FigS8) | Chamber-resolved B-SOiD kinematic profiles (heatmaps) |
| S9 | [`FigS9`](FigS9) | Chamber-resolved B-SOiD bout-duration distributions (14-cluster CDF grid) |
