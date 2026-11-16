# windowed_fst.R -- windowed Hudson FST (Fig. S2) and pcadapt overlap
#
# Hudson FST (Bhatia et al. 2013) per SNP from plink2 allele counts, summed as a
# ratio of sums in 50 kb non-overlapping windows (contigs >= 250 kb, >= 25 SNPs).
# Top 1% of windows per comparison are flagged; pcadapt outliers
# (Bonferroni-adjusted p < 0.01) are flagged per window.
# Six pairwise comparisons plus pooled Pacific vs Gulf.
# Outputs: output/windows_<comparison>.csv, output/FST_genome_wide.csv,
#          output/plots/Figure_S2_windowed_FST.png
source(here::here("scripts", "00_setup.R"))
library(patchwork)

indiv <- read_individuals()
indiv |> count(Localisation_b)   # expect central 53, gulf 15, north 21, south 24

bim <- read_bim()
nrow(bim)                        # expect 1,665,506
bim_pos <- select(bim, chrom, pos)

# allele counts per group (.acount rows are in .bim order)
walk(groups, \(g) {
  indiv |>
    filter(Localisation_b == g) |>
    transmute(`#IID` = Indiv) |>
    write_tsv(here("plink", paste0("keep_", g, ".txt")))
  system2(plink2,
          c("--bfile", here("plink", "GMIall.filtered"),
            "--keep", here("plink", paste0("keep_", g, ".txt")),
            "--freq", "counts",
            "--out", here("plink", paste0("cnt_", g))))
})

counts <- groups |>
  set_names() |>
  map(\(g) {
    here("plink", paste0("cnt_", g, ".acount")) |>
      read_tsv(show_col_types = FALSE) |>
      transmute(alt = ALT_CTS, n = OBS_CT) |>
      bind_cols(bim_pos)
  })

pac <- counts[["north"]] |>
  mutate(
    alt = alt + pull(counts[["central"]], alt) + pull(counts[["south"]], alt),
    n   = n   + pull(counts[["central"]], n)   + pull(counts[["south"]], n)
  )
gulf <- counts[["gulf"]]

pairs <- tribble(
  ~g1,       ~g2,
  "north",   "south",
  "north",   "central",
  "central", "south",
  "north",   "gulf",
  "central", "gulf",
  "south",   "gulf"
)

long_contigs <- bim_pos |>
  group_by(chrom) |>
  summarise(len = max(pos)) |>
  filter(len >= min_contig) |>
  pull(chrom)

window_fst <- function(s) {
  s |>
    filter(chrom %in% long_contigs) |>
    mutate(win = (pos - 1) %/% win_size) |>
    group_by(chrom, win) |>
    summarise(n_snps = n(), fst = sum(num) / sum(den), .groups = "drop") |>
    filter(n_snps >= min_snps, is.finite(fst)) |>
    mutate(start = win * win_size + 1, end = (win + 1) * win_size,
           top1pct = fst >= quantile(fst, 0.99))
}

fst_windows <- pairs |>
  mutate(
    name = unname(comparison_names[paste(g1, g2, sep = "_")]),
    snps = map2(g1, g2, \(x, y) {
      a <- counts[[x]]; b <- counts[[y]]
      bind_cols(select(a, chrom, pos),
                hudson(pull(a, alt), pull(a, n), pull(b, alt), pull(b, n))) |>
        filter(is.finite(num), is.finite(den))
    }),
    genome_fst = map_dbl(snps, \(s) summarise(s, f = sum(num) / sum(den)) |> pull(f)),
    windows = map(snps, window_fst)
  )

# pooled Pacific vs Gulf
pac_gulf_snps <- bind_cols(select(pac, chrom, pos),
                           hudson(pull(pac, alt), pull(pac, n), pull(gulf, alt), pull(gulf, n))) |>
  filter(is.finite(num), is.finite(den))

genome_fst <- bind_rows(
  select(fst_windows, name, genome_fst),
  tibble(name = "pac_vs_gulf",
         genome_fst = summarise(pac_gulf_snps, f = sum(num) / sum(den)) |> pull(f))
)
genome_fst   # expect n_vs_s 0.160, n_vs_c 0.096, c_vs_s 0.078, Gulf comparisons 0.46-0.54
write_csv(genome_fst, here("output", "FST_genome_wide.csv"))

# pcadapt outliers (Bonferroni < 0.01) -> windows containing >= 1 outlier
get_outlier_snps <- function(name) {
  pvals <- here("output", paste0("pcadapt_lcwgs_", name, "_K2.rds")) |>
    readRDS() |>
    pluck("pvalues")
  stopifnot(length(pvals) == nrow(bim))
  bim_pos |>
    mutate(outlier = p.adjust(pvals, method = "bonferroni") < 0.01) |>
    filter(outlier)
}

pcadapt_win <- c("n_vs_s", "n_vs_c", "c_vs_s", "n_vs_gulf", "c_vs_gulf", "s_vs_gulf", "all") |>
  set_names() |>
  map(\(nm) get_outlier_snps(nm) |>
        mutate(win = (pos - 1) %/% win_size) |>
        distinct(chrom, win)) |>
  list_rbind(names_to = "name") |>
  mutate(pcadapt = TRUE,
         name = if_else(name == "all", "pac_vs_gulf", name))

# cumulative x-position along the long contigs
contig_x <- bim_pos |>
  group_by(chrom) |>
  summarise(len = max(pos), .groups = "drop") |>
  filter(chrom %in% long_contigs) |>
  arrange(desc(len)) |>
  mutate(offset = dplyr::lag(cumsum(len), default = 0),
         shade  = row_number() %% 2)

win_tbl <- bind_rows(
  fst_windows |> select(name, windows) |> unnest(windows),
  window_fst(pac_gulf_snps) |> mutate(name = "pac_vs_gulf")
) |>
  left_join(contig_x, by = "chrom") |>
  left_join(pcadapt_win, by = c("name", "chrom", "win")) |>
  mutate(
    pcadapt = replace_na(pcadapt, FALSE),
    x = offset + (start + end) / 2,
    col = case_when(
      top1pct & pcadapt ~ "#D7191C",   # top 1% and pcadapt outlier
      top1pct           ~ "blue",      # top 1%
      shade == 1        ~ "black",
      TRUE              ~ "grey70"
    )
  ) |>
  arrange(top1pct)

# one csv per comparison
win_tbl |>
  select(name, chrom, start, end, n_snps, fst, top1pct, pcadapt) |>
  group_split(name) |>
  walk(\(d) write_csv(d, here("output", paste0("windows_", pull(d, name)[1], ".csv"))))

# summary: median window FST and the top-1% threshold
win_tbl |>
  group_by(name) |>
  summarise(median_fst = median(fst), top1pct_threshold = min(fst[top1pct]),
            n_windows = n(), n_top1pct_pcadapt = sum(top1pct & pcadapt))

# Figure S2: A = six pairwise comparisons, B = all individuals (Pacific vs Gulf)
fst_panel <- function(nm, y = FALSE, x = FALSE, ymax = NULL) {
  p <- ggplot(filter(win_tbl, name == nm), aes(x = x, y = fst, colour = col)) +
    geom_point(size = 0.5) +
    scale_colour_identity() +
    labs(x = if (x) "Genome position" else NULL,
         y = if (y) expression(F[ST]) else NULL) +
    theme_classic()
  if (!is.null(ymax)) {
    p <- p + scale_y_continuous(breaks = seq(0, ymax, by = 0.25)) +
      coord_cartesian(ylim = c(0, ymax))
  }
  p
}

# panel order below = A..G in patchwork; A (n_vs_c, C, D) = Pacific pairs, E-G = Gulf
# comparisons, B = pooled. Adjust this layout to match the published figure.
layout <- "
A#B
CDB
EFG
"
Figure_S2 <- (fst_panel("n_vs_c", y = TRUE) +
              fst_panel("pac_vs_gulf", y = TRUE, x = TRUE, ymax = 1) +
              fst_panel("n_vs_s", y = TRUE) +
              fst_panel("c_vs_s") +
              fst_panel("n_vs_gulf", y = TRUE, x = TRUE, ymax = 1) +
              fst_panel("c_vs_gulf", x = TRUE, ymax = 1) +
              fst_panel("s_vs_gulf", x = TRUE, ymax = 1)) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = list(c("A", "B", "", "", "", "", ""))) &
  theme(text = element_text(size = 8), axis.text = element_text(size = 6.5),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.tag = element_text(face = "bold", size = 12))

ggsave(here("output", "plots", "Figure_S2_windowed_FST.png"), Figure_S2,
       width = 180, height = 90, units = "mm", scale = 1.5)
