# GMI lcWGS Population Genomics Pipeline

Analysis code for: **Population genomics and speciation of a Baja California
disjunct species, the longjaw mudsucker (*Gillichthys mirabilis*)**

Allyson C. Salazar Sawkins, Giacomo Bernardi, 2026
Ecology and Evolution — in review

Archived code: Zenodo, https://doi.org/[ZENODO DOI]
Data: Dryad, https://doi.org/10.5061/dryad.573n5tbq2

---

## Overview
This repository contains code for population genomic analyses of *Gillichthys
mirabilis* using low-coverage whole-genome sequencing (lcWGS) data from 113
individuals across 16 localities.

Primary lcWGS processing and population genomic analyses, including PCA,
admixture (sNMF), isolation by distance (IBD, Mantel test), FST estimation, and
triangle plots, followed the pipeline of Toy et al. (2025)
(https://github.com/Arthur-Oules/EJA_lcwgs_GenPop); those scripts are not
duplicated here. Scripts in this repository cover the analyses specific to this
study: annotation liftover, outlier locus detection (pcadapt), windowed FST,
fixed differences between Pacific and Gulf of California populations,
synonymous/nonsynonymous (dN/dS) tests, and inbreeding/heterozygosity.

---

## Repository contents
| File | Purpose |
|---|---|
| `dnds_pipeline.sh` | Bash pipeline: annotation liftover, VCF filtering, fixed-difference detection, syn/nonsyn classification |
| `syn_nonsyn_from_gff.py` | Codon-level classification of fixed differences between Pacific and Gulf of California populations |
| `expected_nonsyn.py` | Expected proportion of nonsynonymous changes across all coding sites (transition-weighted, κ = 2 × ts/tv) |
| `expected_nonsyn_by_gene.py` | Same, per gene |
| `GMI_lcwgs_outliers.qmd` | R/Quarto: pcadapt outlier detection and functional annotation |
| `00_setup.R` | Shared settings and helper functions (sourced by the other R scripts) |
| `windowed_fst.R` | Hudson FST in 50 kb windows for six pairwise comparisons and pooled Pacific vs Gulf of California (Fig. S2); flags top 1% windows and pcadapt outliers |
| `dnds_test.R` | Fixed Pacific–Gulf of California SNPs, transition/transversion ratio, observed vs expected nonsynonymous proportion (binomial test, dN/dS), per-gene dN/dS |
| `inbreeding_heterozygosity.R` | Per-site and per-individual inbreeding coefficient (F) and observed heterozygosity with PLINK2 (Table S4) |
| `environment.yml` | Conda environment specifying dependencies |
| `config.sh.example` | Template for user-specific file paths |

---

## Requirements
Recreate the conda environment:
```bash
conda env create -f environment.yml
conda activate liftoff
```
Main software: PLINK2, Liftoff, bcftools, Python 3, R (tidyverse, here,
patchwork, pcadapt, vegan, geosphere, rtracklayer, GenomicRanges).

Variant calling was done with snpArcher (commit 765d83a, Snakemake 9.3.3;
GATK 4.1.8.0). Parameters are given in Table S1.

---

## Usage
Run all R scripts from the repository root. They read the filtered PLINK files
(`plink/GMIall.filtered.*`) and `data/individuals_coordinates.csv`.

1. Copy `config.sh.example` to `config.sh` and fill in your local paths.
2. Run the bash pipeline: `bash dnds_pipeline.sh`
3. Run outlier detection: open `GMI_lcwgs_outliers.qmd` in RStudio (writes the
   `output/pcadapt_lcwgs_<comparison>_K2.rds` files).
4. Windowed FST (needs the pcadapt files from step 3):
   `Rscript scripts/windowed_fst.R`
5. dN/dS: run `Rscript scripts/dnds_test.R` once to write the fixed-difference
   files and κ; run the Python classification
   (`syn_nonsyn_from_gff.py`, `expected_nonsyn.py`, `expected_nonsyn_by_gene.py`;
   see `dnds_pipeline.sh`); then run `dnds_test.R` again for the tests.
6. Inbreeding/heterozygosity: `Rscript scripts/inbreeding_heterozygosity.R`

Outputs are written to `output/`.

---

## Citation
If you use this code, please cite:
[full citation will become available once published]

Toy et al. (2025) [full reference] for the lcWGS population genomics pipeline.
