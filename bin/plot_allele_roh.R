#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
cohort_id <- "cohort"
species_map_str <- "c:cerevisiae,e:eubayanus,p:paradoxus,m:mikatae,k:kudriavzevii,u:uvarum,a:arboricola,j:jurei"
bin_dir <- "."

i <- 1
while (i <= length(args)) {
  if (args[i] == "--cohort_id") {
    cohort_id <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--species_map") {
    species_map_str <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--bin_dir") {
    bin_dir <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

sp_pairs <- str_split(species_map_str, ",")[[1]]
species_map <- setNames(
  sapply(sp_pairs, function(x) str_split(x, ":")[[1]][2]),
  sapply(sp_pairs, function(x) str_split(x, ":")[[1]][1])
)

png_dpi <- 300
png_width_px <- 4250
png_height_px <- 750

dir.create("allele_plots", showWarnings = FALSE, recursive = TRUE)
dir.create("roh_plots", showWarnings = FALSE, recursive = TRUE)
dir.create("combined_overlay", showWarnings = FALSE, recursive = TRUE)

chr_len_paths <- list.files("chr_lens", pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
subgenome_tables <- list()

for (cp in chr_len_paths) {
  df_c <- suppressWarnings(read_csv(cp, show_col_types = FALSE, progress = FALSE)) %>%
    mutate(
      length = as.numeric(length),
      prefix = ifelse(str_detect(chr, "_"), str_extract(chr, "^[^_]+"), "c"),
      offset = c(0, cumsum(head(length, -1))),
      end    = offset + length,
      mid    = offset + length / 2
    )
  if (nrow(df_c) > 0) {
    pref <- df_c$prefix[1]
    subgenome_tables[[pref]] <- df_c
  }
}

roh_paths <- list.files("roh_files", pattern = "interest\\.txt$", full.names = TRUE, recursive = TRUE)
allele_paths <- sort(list.files("allele_csvs", pattern = "\\.alleleplot\\.csv$", full.names = TRUE, recursive = TRUE))

allele_pdf_list <- list()
roh_pdf_list <- list()
combined_pdf_list <- list()
summary_rows <- list()

for (ap in allele_paths) {
  fname <- basename(ap)
  # Strip either .<pref>.alleleplot.csv or .alleleplot.csv cleanly
  sample_base <- str_replace(fname, "(\\.[A-Za-z0-9]+)?\\.alleleplot\\.csv$", "")

  depth_df <- suppressWarnings(read_delim(
    ap, delim = ";",
    col_types = cols(chr = col_character(), pos = col_double(), ad = col_character()),
    na = c(".", ""), progress = FALSE
  ))

  if (nrow(depth_df) == 0) next

  pref <- ifelse(str_detect(depth_df$chr[1], "_"), str_extract(depth_df$chr[1], "^[^_]+"), "c")
  if (!(pref %in% names(subgenome_tables))) next

  chr_len <- subgenome_tables[[pref]]
  sp_full <- ifelse(pref %in% names(species_map), paste0("S. ", species_map[[pref]]), pref)
  genome_end <- max(chr_len$end)
  offset_map <- setNames(chr_len$offset, chr_len$chr)

  depth_df <- depth_df %>%
    filter(chr %in% chr_len$chr, !is.na(pos), !is.na(ad)) %>%
    mutate(genomepos = offset_map[chr] + (pos - 1))

  if (nrow(depth_df) == 0) next

  # Fast vectorized AD parsing
  ad_split <- strsplit(depth_df$ad, ",", fixed = TRUE)
  lens <- lengths(ad_split)
  ad_num <- suppressWarnings(as.numeric(unlist(ad_split, use.names = FALSE)))
  site_idx <- rep.int(seq_len(nrow(depth_df)), lens)

  valid <- !is.na(ad_num)
  ad_num <- ad_num[valid]
  site_idx <- site_idx[valid]

  if (length(ad_num) == 0) next

  sums_mat <- rowsum(ad_num, site_idx, reorder = FALSE)
  totals <- numeric(nrow(depth_df))
  totals[as.integer(rownames(sums_mat))] <- sums_mat[, 1]

  site_tot <- totals[site_idx]
  keep <- site_tot > 0 & is.finite(site_tot)

  melted_df <- tibble(
    site_id   = site_idx[keep],
    genomepos = depth_df$genomepos[site_idx[keep]],
    value     = ad_num[keep] / site_tot[keep]
  )

  if (nrow(melted_df) == 0) next

  n_vars <- sum(totals > 0)
  het_vars <- length(unique(melted_df$site_id[melted_df$value >= 0.15 & melted_df$value <= 0.85]))
  het_pct <- ifelse(n_vars > 0, round(100.0 * het_vars / n_vars, 1), 0.0)

  # Match corresponding ROH file for this sample and subgenome
  roh_df <- tibble(Chromosome = character(), Start = numeric(), End = numeric(), genomepos_start = numeric(), genomepos_end = numeric())
  expected_roh_names <- c(
    paste0(sample_base, ".", pref, ".roh.interest.txt"),
    paste0(sample_base, ".roh.interest.txt")
  )
  cand_roh <- roh_paths[basename(roh_paths) %in% expected_roh_names]
  for (rp in cand_roh) {
    tmp_r <- suppressWarnings(read_tsv(rp, show_col_types = FALSE, progress = FALSE))
    if (nrow(tmp_r) > 0 && all(c("Chromosome", "Start", "End") %in% names(tmp_r))) {
      tmp_r <- tmp_r %>%
        filter(Chromosome %in% chr_len$chr) %>%
        mutate(
          Start = as.numeric(Start),
          End = as.numeric(End),
          genomepos_start = offset_map[Chromosome] + pmax(0, Start - 1),
          genomepos_end   = offset_map[Chromosome] + pmax(0, End - 1)
        )
      if (nrow(tmp_r) > 0) {
        roh_df <- tmp_r
        break
      }
    }
  }

  chr_roh_stats <- chr_len %>%
    left_join(
      roh_df %>%
        mutate(seg_len = pmax(0, End - Start + 1)) %>%
        group_by(Chromosome) %>%
        summarise(roh_bp = sum(seg_len, na.rm = TRUE), .groups = "drop"),
      by = c("chr" = "Chromosome")
    ) %>%
    mutate(
      roh_bp = pmin(length, replace_na(roh_bp, 0)),
      roh_pct = round(100.0 * roh_bp / length, 1),
      roh_label = paste0(round(roh_pct, 0), "%")
    )

  total_roh_bp <- sum(chr_roh_stats$roh_bp, na.rm = TRUE)
  genome_roh_pct <- round(100.0 * total_roh_bp / genome_end, 1)

  row_entry <- list(
    Sample_ID = sample_base,
    Subgenome = sp_full,
    Informative_Variants = n_vars,
    Heterozygous_Variants = het_vars,
    Heterozygosity_Pct = het_pct,
    Genome_Wide_ROH_Pct = genome_roh_pct
  )
  for (r_i in seq_len(nrow(chr_roh_stats))) {
    row_entry[[paste0("ROH_Pct_", chr_roh_stats$label[r_i])]] <- chr_roh_stats$roh_pct[r_i]
  }
  summary_rows[[paste0(sample_base, "_", pref)]] <- as_tibble(row_entry)

  tag_suffix <- ifelse(length(subgenome_tables) > 1, paste0(".", pref), "")

  # 1. Enhanced Allele Plot
  p_allele <- ggplot(melted_df, aes(x = genomepos, y = value)) +
    geom_hline(yintercept = c(0.25, 0.75), linetype = "dotted", color = "#6C757D", linewidth = 0.3) +
    geom_hline(yintercept = c(1/3, 2/3), linetype = "dashed", color = "#2A9D8F", linewidth = 0.35) +
    geom_hline(yintercept = 0.50, linetype = "dashed", color = "#1D3557", linewidth = 0.4) +
    geom_point(size = 0.1, color = "black") +
    geom_vline(xintercept = c(chr_len$offset, genome_end), color = "red", linewidth = 0.4) +
    annotate("text", x = chr_len$mid, y = 0.0, label = chr_len$label, vjust = 1.5, size = 5.5, fontface = "bold") +
    scale_x_continuous(limits = c(0, genome_end), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0(
        sample_base, "  |  ", sp_full,
        "  |  Informative Variants: ", format(n_vars, big.mark = ","),
        "  |  Het Variants: ", het_pct, "%  (Guides: 0.50 = 1:1, 0.33/0.67 = 1:2, 0.25/0.75 = 1:3)"
      ),
      x = "", y = "Allele Ratio"
    ) +
    theme_gray(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.01),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 6, r = 10, b = 18, l = 6, unit = "pt")
    )

  ggsave(
    filename = paste0(sample_base, tag_suffix, ".alleleplot.png"),
    plot = p_allele, path = "allele_plots",
    dpi = png_dpi, width = png_width_px, height = png_height_px, units = "px", bg = "white"
  )
  allele_pdf_list[[paste0(sample_base, "_", pref)]] <- p_allele

  # 2. Enhanced ROH Plot
  p_roh <- ggplot() +
    {
      if (nrow(roh_df) > 0) {
        geom_rect(data = roh_df, aes(xmin = genomepos_start, xmax = genomepos_end, ymin = 0, ymax = 1), fill = "red", alpha = 0.15)
      }
    } +
    geom_vline(xintercept = c(chr_len$offset, genome_end), color = "red", linewidth = 0.3) +
    annotate("text", x = chr_roh_stats$mid, y = 0.93, label = chr_roh_stats$roh_label, size = 3.3, fontface = "bold", color = "#9B2226") +
    annotate("text", x = chr_len$mid, y = 0, label = chr_len$label, vjust = 1.5, size = 5.5, fontface = "bold") +
    scale_x_continuous(limits = c(0, genome_end), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0(
        sample_base, "  |  ", sp_full,
        "  |  Genome-Wide ROH: ", genome_roh_pct, "%  (",
        format(round(total_roh_bp / 1e6, 2), nsmall = 2), " / ",
        format(round(genome_end / 1e6, 2), nsmall = 2), " Mb)"
      ),
      x = "", y = "ROH"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.01),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), plot.background  = element_blank(),
      legend.position  = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      plot.margin      = margin(t = 6, r = 10, b = 18, l = 6, unit = "pt")
    )

  ggsave(
    filename = paste0(sample_base, tag_suffix, ".roh.interest.png"),
    plot = p_roh, path = "roh_plots",
    dpi = png_dpi, width = png_width_px, height = png_height_px, units = "px", bg = "transparent"
  )
  roh_pdf_list[[paste0(sample_base, "_", pref)]] <- p_roh

  # 3. Combined Overlay Plot
  p_comb <- ggplot() +
    {
      if (nrow(roh_df) > 0) {
        geom_rect(data = roh_df, aes(xmin = genomepos_start, xmax = genomepos_end, ymin = 0, ymax = 1), fill = "red", alpha = 0.16)
      }
    } +
    geom_hline(yintercept = c(0.25, 0.75), linetype = "dotted", color = "#6C757D", linewidth = 0.3) +
    geom_hline(yintercept = c(1/3, 2/3), linetype = "dashed", color = "#2A9D8F", linewidth = 0.35) +
    geom_hline(yintercept = 0.50, linetype = "dashed", color = "#1D3557", linewidth = 0.4) +
    geom_point(data = melted_df, aes(x = genomepos, y = value), size = 0.1, color = "black") +
    geom_vline(xintercept = c(chr_len$offset, genome_end), color = "red", linewidth = 0.4) +
    annotate("text", x = chr_roh_stats$mid, y = 1.04, label = chr_roh_stats$roh_label, size = 3.1, fontface = "bold", color = "#9B2226") +
    annotate("text", x = chr_len$mid, y = 0, label = chr_len$label, vjust = 1.5, size = 5.5, fontface = "bold") +
    scale_x_continuous(limits = c(0, genome_end), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0(
        sample_base, "  |  ", sp_full,
        "  |  Variants: ", format(n_vars, big.mark = ","),
        " (Het: ", het_pct, "%)  |  Genome-Wide ROH: ", genome_roh_pct, "%"
      ),
      x = "", y = "Allele Ratio & ROH"
    ) +
    theme_gray(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.01),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      plot.margin = margin(t = 12, r = 10, b = 18, l = 6, unit = "pt")
    )

  ggsave(
    filename = paste0(sample_base, tag_suffix, ".allele_roh_combined.png"),
    plot = p_comb, path = "combined_overlay",
    dpi = png_dpi, width = png_width_px, height = png_height_px, units = "px", bg = "white"
  )
  combined_pdf_list[[paste0(sample_base, "_", pref)]] <- p_comb
}

if (length(allele_pdf_list) > 0) {
  pdf("cohort_allele_plots.pdf", width = 14.17, height = 2.8)
  for (p in allele_pdf_list) print(p)
  dev.off()
}
if (length(roh_pdf_list) > 0) {
  pdf("cohort_roh_plots.pdf", width = 14.17, height = 2.8)
  for (p in roh_pdf_list) print(p)
  dev.off()
}
if (length(combined_pdf_list) > 0) {
  pdf("cohort_allele_roh_combined.pdf", width = 14.17, height = 2.8)
  for (p in combined_pdf_list) print(p)
  dev.off()
}

summary_df <- bind_rows(summary_rows)
if (nrow(summary_df) == 0) {
  summary_df <- tibble(
    Sample_ID = character(), Subgenome = character(),
    Informative_Variants = numeric(), Heterozygous_Variants = numeric(),
    Heterozygosity_Pct = numeric(), Genome_Wide_ROH_Pct = numeric()
  )
}
write_tsv(summary_df, "cohort_allele_roh_summary.tsv")

# Safe, non-fatal Excel workbook export via pure Python 3 helper
tryCatch({
  py_cmd <- sprintf(
    "import csv, sys; sys.path.insert(0, '%s'); from summarize_coverage import write_multi_sheet_xlsx; rows = list(csv.reader(open('cohort_allele_roh_summary.tsv'), delimiter='\\t')); rows_parsed = [[(float(x) if i>=2 and x.replace('.','',1).isdigit() else x) for i,x in enumerate(r)] for r in rows]; write_multi_sheet_xlsx('cohort_allele_roh_summary.xlsx', {'Cohort_Allele_ROH_Summary': rows_parsed}, color_scale_sheets={'Cohort_Allele_ROH_Summary': (2, 0.0, 50.0, 100.0)})",
    bin_dir
  )
  system2("python3", args = c("-c", shQuote(py_cmd)), stdout = FALSE, stderr = FALSE)
}, error = function(e) {
  message("Note: Optional Excel export skipped (", e$message, ")")
})