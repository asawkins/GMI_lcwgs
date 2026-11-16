# 00_setup.R -- shared settings and helpers. Sourced by the other R scripts.
# Run all scripts from the repository root (here() finds it automatically).
#
# Expected inputs (see README):
#   data/individuals_coordinates.csv   individual metadata (";" delimited)
#   plink/GMIall.filtered.{bed,bim,fam} filtered SNP set (1,665,506 biallelic SNPs before MAF)
#   output/pcadapt_lcwgs_<comparison>_K2.rds   pcadapt results (GMI_lcwgs_outliers.qmd)
# Optional: set PLINK2 in the environment, otherwise "plink2" on the PATH is used.
library(tidyverse)
library(here)

plink2 <- Sys.getenv("PLINK2", unset = "plink2")

# Directory for dN/dS intermediate files (written by dnds_pipeline.sh).
# Override with the environment variable DN_DIR.
dn_dir <- Sys.getenv("DN_DIR", unset = here("output", "dnds"))

win_size   <- 50000    # window size (bp)
min_contig <- 250000   # only contigs >= 250 kb are windowed
min_snps   <- 25       # minimum SNPs per window
groups     <- c("north", "central", "south", "gulf")

comparison_names <- c(
  north_south = "n_vs_s", north_central = "n_vs_c", central_south = "c_vs_s",
  north_gulf = "n_vs_gulf", central_gulf = "c_vs_gulf", south_gulf = "s_vs_gulf"
)

dir.create(here("output", "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(dn_dir, recursive = TRUE, showWarnings = FALSE)

# Hudson FST numerator and denominator per SNP (Bhatia et al. 2013)
hudson <- function(alt1, n1, alt2, n2) {
  p1 <- alt1 / n1
  p2 <- alt2 / n2
  tibble(
    num = (p1 - p2)^2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1),
    den = p1 * (1 - p2) + p2 * (1 - p1)
  )
}

# individuals, in .fam order
read_individuals <- function() {
  indiv <- here("data", "individuals_coordinates.csv") |>
    read_delim(delim = ";", locale = locale(decimal_mark = "."), show_col_types = FALSE)
  fam <- here("plink", "GMIall.filtered.fam") |>
    read_table(col_names = FALSE, show_col_types = FALSE)
  slice(indiv, match(pull(fam, 2), Indiv))
}

read_bim <- function() {
  here("plink", "GMIall.filtered.bim") |>
    read_tsv(col_names = c("chrom", "id", "cm", "pos", "a1", "a2"),
             col_types = "cccdcc")
}
