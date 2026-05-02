# GMI lcWGS Population Genomics Pipeline

Analysis code for: **Population genomics and speciation of a Baja California 
disjunct species, the longjaw mudsucker (*Gillichthys mirabilis*)**

Allyson C. Salazar Sawkins, Giacomo Bernardi, 2026 
Ecology and Evolution — in review

---

## Overview
This repository contains code for population genomic analyses of *Gillichthys 
mirabilis* using low-coverage whole-genome sequencing (lcWGS) data from 113 
individuals across 16 localities.

Primary lcWGS processing and population genomic analyses — including PCA, 
admixture (sNMF), isolation by distance (IBD), FST estimation, and triangle 
plots — followed the pipeline of Toy et al. (2025) 
(https://github.com/Arthur-Oules/EJA_lcwgs_GenPop). Scripts in this 
repository cover Nanopore reference genome assembly and downstream analyses 
specific to this study, including annotation liftover, outlier locus 
detection, fixed difference detection between Pacific and Sea of Cortez 
populations, and synonymous/nonsynonymous substitution classification.

Data are deposited on Dryad: https://doi.org/10.5061/dryad.573n5tbq2

---

## Repository contents
- `dnds_pipeline.sh` — bash pipeline for annotation liftover, VCF filtering, 
   fixed difference detection, and syn/nonsyn classification
- `syn_nonsyn_from_gff.py` — Python script for codon-level classification of 
   fixed differences between Pacific and Sea of Cortez populations
- `GMI_lcwgs_outliers.qmd` — R/Quarto markdown for outlier locus detection 
   and functional annotation
- `environment.yml` — conda environment file specifying all dependencies
- `config.sh.example` — template for user-specific file paths

---

## Requirements
Recreate the conda environment:
```bash
conda env create -f environment.yml
conda activate liftoff
```

---

## Usage
1. Copy `config.sh` and fill in your local paths
2. Run the bash pipeline: `bash dnds_pipeline.sh`
3. Run the R analysis: open `GMI_lcwgs_outliers.qmd` in RStudio

---

## Citation
If you use this code, please cite:
[full citation once published]

