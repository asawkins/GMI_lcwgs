# dnds_test.R -- dN/dS test for fixed Pacific-Gulf differences (+ per-gene test)
#
# Order of operations (see README):
#   1. Rscript scripts/dnds_test.R   (part 1: writes fixed-difference files and kappa)
#   2. python3 scripts/syn_nonsyn_from_gff.py ...   (classifies fixed SNPs, see dnds_pipeline.sh)
#      python3 scripts/expected_nonsyn.py <gff> <fasta> <kappa>            -> expected proportion
#      python3 scripts/expected_nonsyn_by_gene.py <gff> <fasta> <kappa> <out.tsv>
#   3. Rscript scripts/dnds_test.R   (part 2 runs automatically once the classified file exists)
# Set EXP_P to the value printed by expected_nonsyn.py (0.7388 in the paper).
source(here::here("scripts", "00_setup.R"))

# ---- part 1: fixed differences (expect 26,458 SNPs) -------------------------
indiv <- read_individuals()
bim   <- read_bim()

count_group <- function(g) {
  keep <- here("plink", paste0("keep_", g, ".txt"))
  indiv |> filter(Localisation_b %in% g) |> transmute(`#IID` = Indiv) |> write_tsv(keep)
  out <- here("plink", paste0("cnt_", paste(g, collapse = "_")))
  system2(plink2, c("--bfile", here("plink", "GMIall.filtered"), "--keep", keep,
                    "--freq", "counts", "--out", out))
  read_tsv(paste0(out, ".acount"), show_col_types = FALSE) |>
    transmute(alt = ALT_CTS, n = OBS_CT) |>
    bind_cols(select(bim, chrom, pos, a1, a2))
}
pac  <- count_group(c("north", "central", "south"))
gulf <- count_group("gulf")

fixed_alleles <- inner_join(pac, gulf, by = c("chrom", "pos", "a1", "a2"),
                            suffix = c("_pac", "_gulf")) |>
  filter(n_pac > 0, n_gulf > 0) |>
  mutate(f_pac = alt_pac / n_pac, f_gulf = alt_gulf / n_gulf) |>
  filter((f_pac == 0 & f_gulf == 1) | (f_pac == 1 & f_gulf == 0)) |>
  mutate(REF = a2, ALT = a1,
         PAC = if_else(f_pac == 0, REF, ALT),
         GULF = if_else(f_gulf == 0, REF, ALT)) |>
  transmute(CHROM = chrom, POS = pos, REF, ALT, PAC, GULF)

nrow(fixed_alleles)   # expect 26,458
write_tsv(fixed_alleles, file.path(dn_dir, "PACvsGULF.fixed_alleles_ALL.tsv"))
fixed_alleles |>
  mutate(id = paste0(CHROM, ":", POS)) |>
  pull(id) |>
  write_lines(file.path(dn_dir, "fixed_ALL.ids.txt"))

# transition/transversion ratio -> kappa = 2 x ts/tv (~1.66), used by expected_nonsyn.py
tstv <- fixed_alleles |>
  mutate(ts = (REF %in% c("A", "G") & ALT %in% c("A", "G")) |
              (REF %in% c("C", "T") & ALT %in% c("C", "T"))) |>
  summarise(ts_tv = sum(ts) / sum(!ts))
kappa <- 2 * pull(tstv, ts_tv)
message("ts/tv = ", round(pull(tstv, ts_tv), 3), "; kappa = ", round(kappa, 3))

# ---- part 2: observed vs expected nonsynonymous proportion ------------------
classified <- file.path(dn_dir, "fixed_ALL.syn_nonsyn.tsv")
if (!file.exists(classified)) {
  message("Run the Python classification steps first, then re-run this script.")
  quit(save = "no")
}
exp_p <- as.numeric(Sys.getenv("EXP_P", unset = "0.7388"))

fixed_cds <- read_tsv(classified, show_col_types = FALSE) |>
  distinct(sid, .keep_all = TRUE) |>
  filter(class %in% c("nonsynonymous", "synonymous"))

n_nonsyn <- sum(pull(fixed_cds, class) == "nonsynonymous")
n_total  <- nrow(fixed_cds)       # expect 433 (246 nonsynonymous, 187 synonymous)
bt <- binom.test(n_nonsyn, n_total, p = exp_p)   # expect p = 1.9e-14
bt

odds <- \(p) p / (1 - p)
omega_tbl <- tibble(
  n_nonsyn = n_nonsyn, n_syn = n_total - n_nonsyn, observed_prop = n_nonsyn / n_total,
  expected_prop = exp_p, p_value = pluck(bt, "p.value"),
  omega = odds(n_nonsyn / n_total) / odds(exp_p),    # expect 0.47
  lower = odds(pluck(bt, "conf.int")[1]) / odds(exp_p),  # 0.38
  upper = odds(pluck(bt, "conf.int")[2]) / odds(exp_p)   # 0.57
)
omega_tbl
write_csv(omega_tbl, here("output", "dNdS_overall.csv"))

# ---- part 3: per-gene test (expect 388 genes, none significant after BH) ----
expected_gene <- file.path(dn_dir, "expected_nonsyn_by_gene.tsv")
if (file.exists(expected_gene)) {
  expected_gene <- read_tsv(expected_gene, show_col_types = FALSE) |>
    mutate(exp_p = nonsyn_w / (nonsyn_w + syn_w))

  gene_dnds <- fixed_cds |>
    group_by(gene) |>
    summarise(n_nonsyn = sum(class == "nonsynonymous"),
              n_syn = sum(class == "synonymous"), .groups = "drop") |>
    mutate(n = n_nonsyn + n_syn) |>
    left_join(select(expected_gene, gene, nonsyn_w, syn_w, exp_p), by = "gene") |>
    filter(!is.na(exp_p)) |>
    mutate(
      # 0.5 added to each count so genes with zero synonymous SNPs do not give infinity
      omega = ((n_nonsyn + 0.5) / (n_syn + 0.5)) / (nonsyn_w / syn_w),
      p_greater = pmap_dbl(list(n_nonsyn, n, exp_p),
                           \(k, n, p) binom.test(k, n, p = p, alternative = "greater") |>
                             pluck("p.value")),
      fdr = p.adjust(p_greater, method = "BH")
    ) |>
    arrange(p_greater)

  gene_dnds |>
    summarise(genes = n(), snps = sum(n), median_snps = median(n),
              genes_omega_gt1 = sum(omega > 1), genes_fdr_05 = sum(fdr < 0.05))
  write_csv(gene_dnds, here("output", "dNdS_by_gene.csv"))
}
