# inbreeding_heterozygosity.R -- per-site / per-individual inbreeding (Table S4)
#
# PLINK2 --het with allele frequencies estimated within each site
# (--read-freq), so F is relative to that site. Sites with very few individuals
# (e.g. Alviso n = 2, Bodega Bay n = 4) give unreliable F.
# Observed heterozygosity = 1 - O(HOM) / OBS_CT.
# Outputs: output/Table_S_inbreeding_by_individual.csv, ..._by_site.csv
source(here::here("scripts", "00_setup.R"))

indiv <- read_individuals() |>
  mutate(site = str_extract(Indiv, "(?<=GMI_)[A-Z]+"))
het_dir <- here("output", "het")
dir.create(het_dir, recursive = TRUE, showWarnings = FALSE)
bfile <- here("plink", "GMIall.filtered")

het_one_site <- function(s) {
  keep <- indiv |> filter(site == s) |> transmute(`#IID` = Indiv)
  keep_file <- file.path(het_dir, paste0(s, ".keep"))
  write_tsv(keep, keep_file)
  out <- file.path(het_dir, s)
  # --set-all-var-ids gives unique variant IDs, which --read-freq needs
  common <- c("--bfile", bfile, "--allow-extra-chr", "--keep", keep_file,
              "--set-all-var-ids", "@:#")
  system2(plink2, c(common, "--freq", "--out", out))
  system2(plink2, c(common, "--read-freq", paste0(out, ".afreq"), "--het", "--out", out))
  read_tsv(paste0(out, ".het"), show_col_types = FALSE) |>
    rename_with(\(x) str_remove(x, "^#")) |>
    mutate(site = s)
}

het_ind <- indiv |>
  distinct(site) |>
  pull(site) |>
  map(het_one_site) |>
  list_rbind() |>
  left_join(select(indiv, IID = Indiv, group = Localisation_b), by = "IID") |>
  transmute(IID, site, group, n_snps = OBS_CT, obs_het = 1 - `O(HOM)` / OBS_CT, F)

het_site <- het_ind |>
  group_by(site, group) |>
  summarise(n = n(), mean_obs_het = mean(obs_het), mean_F = mean(F), .groups = "drop")

het_site   # expect Elkhorn Slough (ELK): mean F about -0.055, observed het about 0.363

write_csv(het_ind,  here("output", "Table_S_inbreeding_by_individual.csv"))
write_csv(het_site, here("output", "Table_S_inbreeding_by_site.csv"))
