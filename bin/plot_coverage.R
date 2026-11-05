#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
species_map_str <- "c:cerevisiae,e:eubayanus,p:paradoxus,m:mikatae,k:kudriavzevii,u:uvarum,a:arboricola,j:jurei"
ploidy_mult <- 2.0
ploidy_tsv <- ""
files <- c()

i <- 1
while (i <= length(args)) {
  if (args[i] == "--species_map") {
    species_map_str <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--ploidy_mult") {
    ploidy_mult <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--ploidy_tsv") {
    ploidy_tsv <- args[i + 1]
    i <- i + 2
  } else {
    files <- c(files, args[i])
    i <- i + 1
  }
}

curated_g22 <- list()
curated_g24 <- list()
if (nchar(ploidy_tsv) > 0 && file.exists(ploidy_tsv)) {
  p_df <- read_tsv(ploidy_tsv, show_col_types = FALSE)
  if (nrow(p_df) > 0) {
    for (r_i in seq_len(nrow(p_df))) {
      s_id <- as.character(p_df$Sample_ID[r_i])
      curated_g22[[s_id]] <- suppressWarnings(as.numeric(p_df$G22_Baseline[r_i]))
      curated_g24[[s_id]] <- suppressWarnings(as.numeric(p_df$G24_Ploidy[r_i]))
    }
  }
}

species_pairs <- str_split(species_map_str, ",")[[1]]
species_map <- setNames(
  sapply(species_pairs, function(x) str_split(x, ":")[[1]][2]),
  sapply(species_pairs, function(x) str_split(x, ":")[[1]][1])
)

species_colors <- c(
  "c" = "#4F81BD", "e" = "#00B050", "p" = "#E07A5F", "m" = "#81B29A",
  "k" = "#9B5DE5", "u" = "#2A9D8F", "a" = "#F4A261", "j" = "#E76F51"
)

roman_levels <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")

dir.create("karyograms", showWarnings = FALSE, recursive = TRUE)

build_strain_plot <- function(df_nuc, sample_id, g22, g24_val) {
  y_max <- max(8.0, 2.5 * g24_val)
  y_step <- ifelse(y_max >= 8, 2, 1)

  # Lock subgenome facet order to canonical species_map order (c, e, p, m, k, u, a, j)
  ordered_labels <- df_nuc %>%
    distinct(prefix_order, species_label) %>%
    arrange(prefix_order) %>%
    pull(species_label)

  df_plot <- df_nuc %>%
    mutate(
      species_label = factor(species_label, levels = ordered_labels),
      norm_val = pmin((cov / g22) * g24_val, y_max)
    ) %>%
    arrange(prefix_order, chrom_order, start) %>%
    group_by(prefix, species_label) %>%
    mutate(win_idx = row_number()) %>%
    ungroup()

  bounds <- df_plot %>%
    group_by(prefix, species_label, chrom_short, chrom_order) %>%
    summarise(
      start_idx = min(win_idx),
      end_idx = max(win_idx),
      mid_idx = (min(win_idx) + max(win_idx)) / 2,
      .groups = "drop"
    )

  n_subgenomes <- n_distinct(df_plot$prefix)

  p <- ggplot(df_plot, aes(x = win_idx, y = norm_val, fill = prefix)) +
    geom_col(width = 1.0, show.legend = FALSE) +
    geom_vline(
      data = bounds,
      aes(xintercept = end_idx),
      linetype = "dashed", color = "#333333", linewidth = 0.35
    ) +
    geom_text(
      data = bounds,
      aes(x = mid_idx, y = -0.06 * y_max, label = chrom_short),
      inherit.aes = FALSE, size = 2.6, color = "#222222"
    ) +
    facet_wrap(~ species_label, ncol = 1, scales = "free_x") +
    scale_fill_manual(values = species_colors) +
    scale_y_continuous(
      limits = c(-0.10 * y_max, y_max),
      breaks = seq(0, y_max, by = y_step),
      expand = c(0, 0)
    ) +
    scale_x_continuous(expand = c(0.002, 0.002)) +
    labs(
      title = sample_id,
      x = "Chromosomes (Concatenated Sliding Windows)",
      y = "Copy Number"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      strip.text = element_text(face = "bold.italic", size = 10, hjust = 0),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.8, "lines")
    )

  return(list(plot = p, n_sub = n_subgenomes))
}

cohort_plots <- list()
max_cohort_subgenomes <- 2

for (fpath in sort(files)) {
  sample_id <- str_replace(basename(fpath), "\\.slidingwindow\\.tab$", "")
  df <- read_tsv(fpath, col_names = c("chrom", "start", "end", "cov"), show_col_types = FALSE) %>%
    mutate(cov = replace_na(suppressWarnings(as.numeric(cov)), 0.0))

  df_nuc <- df %>%
    filter(!str_ends(chrom, "_M") & !str_ends(chrom, "_P")) %>%
    mutate(
      prefix = ifelse(str_detect(chrom, "_"), str_extract(chrom, "^[^_]+"), "c"),
      chrom_short = ifelse(str_detect(chrom, "_"), str_replace(chrom, "^[^_]+_", ""), str_replace(chrom, "^chr", "")),
      chrom_order = match(chrom_short, roman_levels, nomatch = 50),
      prefix_order = match(prefix, names(species_map), nomatch = 99),
      species_name = ifelse(prefix %in% names(species_map), species_map[prefix], prefix),
      species_label = paste0("S. ", species_name)
    )

  chrom_meds <- df_nuc %>%
    group_by(prefix, chrom) %>%
    summarise(med_cov = median(cov, na.rm = TRUE), .groups = "drop")

  non_zero_meds <- chrom_meds %>% filter(med_cov > 0) %>% pull(med_cov)
  q75 <- ifelse(length(non_zero_meds) > 0, quantile(non_zero_meds, 0.75), 10.0)
  noise_cutoff <- max(5.0, 0.20 * q75)

  # Keep subgenomes with >=1 active chromosome OR >=8 active windows (captures segmental introgressions)
  chr_active_pref <- chrom_meds %>%
    filter(med_cov >= noise_cutoff) %>%
    pull(prefix) %>%
    unique()

  win_active_pref <- df_nuc %>%
    group_by(prefix) %>%
    summarise(n_active_win = sum(cov >= noise_cutoff, na.rm = TRUE), .groups = "drop") %>%
    filter(n_active_win >= 8) %>%
    pull(prefix)

  active_prefixes <- union(chr_active_pref, win_active_pref)
  if (length(active_prefixes) == 0) {
    active_prefixes <- unique(df_nuc$prefix)[1]
  }

  df_nuc_active <- df_nuc %>% filter(prefix %in% active_prefixes)

  sc_active  <- chrom_meds %>% filter(prefix == "c", med_cov >= noise_cutoff) %>% pull(med_cov)
  all_active <- chrom_meds %>% filter(med_cov >= noise_cutoff) %>% pull(med_cov)

  default_g22 <- if (length(sc_active) > 0) median(sc_active) else if (length(all_active) > 0) median(all_active) else 1.0

  g22_use <- default_g22
  if (sample_id %in% names(curated_g22)) {
    val_g22 <- curated_g22[[sample_id]]
    if (!is.na(val_g22) && val_g22 > 0) g22_use <- val_g22
  }

  g24_use <- ploidy_mult
  if (sample_id %in% names(curated_g24)) {
    val_g24 <- curated_g24[[sample_id]]
    if (!is.na(val_g24) && val_g24 > 0) g24_use <- val_g24
  }

  res <- build_strain_plot(df_nuc_active, sample_id, g22_use, g24_use)
  if (res$n_sub > max_cohort_subgenomes) {
    max_cohort_subgenomes <- res$n_sub
  }

  ggsave(
    filename = file.path("karyograms", paste0(sample_id, "_karyogram.png")),
    plot = res$plot, width = 14, height = max(3.0, 2.2 * res$n_sub), dpi = 200, bg = "white"
  )
  cohort_plots[[sample_id]] <- res$plot
}

pdf("cohort_karyograms.pdf", width = 14, height = max(5.5, 2.2 * max_cohort_subgenomes))
for (p in cohort_plots) print(p)
dev.off()