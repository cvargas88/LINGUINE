# LINGUINE 🍝
**LINkage GroUps INfErence for Ancestral Genomes**

[![R CMD Check](https://github.com/cvargas88/LINGUINE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cvargas88/LINGUINE/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Reconstructing the chromosomal architecture of ancestral genomes is a central challenge in comparative genomics. Existing methods struggle with rapid genome evolution, extensive chromosomal rearrangements, and whole-genome duplications (WGDs) followed by gene loss. 

**LINGUINE** is a phylogeny-aware R pipeline designed to solve these limitations. Instead of relying purely on single-copy orthologs, LINGUINE leverages full orthogroups to maximize marker density. It employs a Hidden Markov Model (HMM) framework to seamlessly bridge syntenic blocks across gene loss, micro-inversions, and assembly fragmentation.

## ✨ Key Features
* **Paralogy-Aware:** Automatically resolves paralogy and collapses duplicated linkage groups arising from WGDs prior to ancestral reconstruction.
* **Orthogroup Integration:** Maximizes evolutionary signals by using OrthoFinder Hierarchical Orthogroups (HOGs) rather than strict single-copy genes.
* **Iterative Post-Order Traversal:** Progressively integrates information from all descendant lineages up the phylogenetic tree, making inferences highly robust to missing data in individual taxa.
* **Automated Visualization:** Natively generates Oxford Grids, stacked distributions, and chromosomal ideograms for every reconstructed node.

## ⚠️ Critical Data Requirement
**Your GFF feature IDs must perfectly match your OrthoFinder sequence IDs!**
If your genome annotation pipeline appends suffixes (e.g., `-PA`, `_t1`) to transcripts, but your OrthoFinder dictionary lacks them, LINGUINE will fail to bridge the physical coordinates with the evolutionary orthology. Ensure your GFF `ID=` attributes and OrthoFinder `.tsv` gene names are strictly harmonized before running the pipeline.

## 🛠 Installation

You can install the development version of LINGUINE from GitHub using `devtools`:

```r
# install.packages("devtools")
devtools::install_github("cvargas88/LINGUINE")
```

## 🚀 Quick Start

LINGUINE is operated entirely through a centralized configuration object. 

```r
library(LINGUINE)

# 1. Generate the configuration object
config <- create_linguine_config(
  dataset = "Nematoda_Analysis",
  base_dir = "/path/to/your/working/directory",
  tree_filename = "species_tree.nwk",
  orthology_type = "HOGs",
  orthology_filename = "N0.tsv",
  min_chromosome_length_bp = 4500000, 
  resolve_multimapped = "drop" # Options: "drop", "keep", "random"
)

# 2. Execute the entire pipeline
run_linguine(config)
```

## 📊 Pipeline Outputs
All outputs are automatically sorted into the directory defined in your `base_dir`:
* `processed_data/`: Sanitized gene coordinates, filtered chromosome sizes, and centralized orthology dictionaries.
* `intermediate_data/`: HMM emission matrices and Viterbi state classifications.
* `results/`: The finalized `.rds` objects for the Ancestral Genomes at every internal node.
* `plots/`: High-resolution PNGs and SVGs of phylogenetic trees, Oxford Grids, and physical karyotypes.

## 📖 Citation
If you use LINGUINE in your research, please cite:
> Vargas-Chávez, C., & Fernández, R. (202X). *[Paper Title]*. [Journal Name]. [DOI]
