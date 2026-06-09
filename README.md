# Mouse-Prematurity-Insults

Code and analysis notebooks for **"MIA-Induced Metabolic Priming Increases Vulnerability to Neonatal Hypoxia"**. Mouse cerebellar phenotyping across maternal immune activation (MIA), neonatal hypoxia (Hx), and sequential insults at P11 and P45 — bioenergetics, ultrastructural and synaptic morphometry, spatial transcriptomics, and behavior.

## Repository structure

| Path | Contents |
|---|---|
| [`Figures/`](Figures) | Code and rendered notebooks (`.html`) for every main and supplementary figure panel. See [`Figures/README.md`](Figures/README.md) for the full panel-by-panel index. |
| [`CoGAPS/`](CoGAPS) | CoGAPS pattern factorization of the spatial transcriptome (feeds Fig 3B, Fig S3, Fig S4). |
| [`LIANA/`](LIANA) | LIANA ligand–receptor / cell–cell communication analysis. |
| [`Tricycle/`](Tricycle) | Tricycle cell-cycle analysis (feeds Fig 5). |

Each figure-panel folder contains the analysis source (`.Rmd` for R, `.ipynb` for Python) alongside its rendered `.html`, which shows the code, output, and final panel inline.

## Data availability

Raw and processed data are **not** included in this repository. Analysis scripts read
from a local `data/` directory that is not distributed; the committed `.html` outputs
embed all intermediate results and the final figures so the analysis is fully reviewable
without the underlying data. Spatial transcriptomics data analysed in this study were
retrieved from the NoCodeSeg repository. Additional data are available from the
corresponding author on reasonable request (see the manuscript's *Data availability*
statement).

## Citation

If you use this code, please cite the manuscript (citation details to follow).
