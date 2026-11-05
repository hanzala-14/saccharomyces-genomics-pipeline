#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
cohort_id   <- "Cohort"
species_map <- "c:cerevisiae,e:eubayanus,p:paradoxus,m:mikatae,k:kudriavzevii,u:uvarum,a:arboricola,j:jurei"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--cohort_id" && i < length(args)) {
    cohort_id <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--species_map" && i < length(args)) {
    species_map <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

sp_lookup <- setNames(
  sapply(strsplit(strsplit(species_map, ",")[[1]], ":"), `[`, 2),
  sapply(strsplit(strsplit(species_map, ",")[[1]], ":"), `[`, 1)
)

dir.create("plots", showWarnings = FALSE, recursive = TRUE)

win_files <- sort(list.files("window_tsvs", pattern = "\\.windows\\.tsv$", full.names = TRUE))

all_chrom_summaries <- list()
all_win_summaries   <- list()

for (wf in win_files) {
  base_tag <- sub("\\.windows\\.tsv$", "", basename(wf))
  parts    <- strsplit(base_tag, "__")[[1]]
  pref     <- parts[1]
  c_slug   <- ifelse(length(parts) > 1, parts[2], "All_Strains")

  sp_raw   <- ifelse(pref %in% names(sp_lookup), sp_lookup[[pref]], pref)
  sp_cap   <- paste0(toupper(substring(sp_raw, 1, 1)), substring(sp_raw, 2))
  sp_ital  <- paste0("S. ", sp_raw)

  cf <- file.path("chrom_tsvs", paste0(base_tag, ".chrom_summary.tsv"))
  lf <- file.path("chr_lens",   paste0(pref, ".chr_length.csv"))
  if (!file.exists(cf) || !file.exists(lf)) next

  chr_lens <- read_csv(lf, show_col_types = FALSE) %>%
    mutate(
      cum_end   = cumsum(as.numeric(length)),
      cum_start = cum_end - as.numeric(length),
      cum_mid   = cum_start + as.numeric(length) / 2,
      shade     = row_number() %% 2 == 0
    )

  win_df <- read_tsv(wf, show_col_types = FALSE) %>%
    mutate(Tajima_D = suppressWarnings(as.numeric(Tajima_D))) %>%
    inner_join(chr_lens %>% select(chr, cum_start), by = c("Chromosome" = "chr")) %>%
    mutate(
      Genome_Pos = cum_start + Midpoint,
      Subgenome  = sp_cap
    )

  chr_df <- read_tsv(cf, show_col_types = FALSE) %>%
    mutate(
      Tajima_D  = suppressWarnings(as.numeric(Tajima_D)),
      Subgenome = sp_cap,
      Prefix    = pref
    )

  clade_label <- unique(win_df$Clade)[1]
  if (is.na(clade_label) || clade_label == "") clade_label <- "All Strains"

  # Format clean header: e.g., "Cerevisiae (Saaz Clade)" or "Cerevisiae"
  if (tolower(clade_label) %in% c("all strains", "all", "cohort")) {
    header_prefix <- sp_cap
    file_tag      <- sp_raw
  } else {
    header_prefix <- paste0(sp_cap, " (", clade_label, ")")
    file_tag      <- paste0(sp_raw, "_", gsub("[^A-Za-z0-9._-]+", "_", clade_label))
  }

  all_chrom_summaries[[base_tag]] <- chr_df
  all_win_summaries[[base_tag]]   <- win_df

  # 1. Sliding-window Pi & Theta_W Track
  p_pi <- ggplot() +
    geom_rect(
      data = chr_lens %>% filter(shade),
      aes(xmin = cum_start, xmax = cum_end, ymin = -Inf, ymax = Inf),
      fill = "#F2F4F7", alpha = 0.8
    ) +
    geom_vline(data = chr_lens, aes(xintercept = cum_end), color = "#D0D5DD", linetype = "dashed", linewidth = 0.3) +
    geom_line(data = win_df, aes(x = Genome_Pos, y = Pi_per_bp * 1000, color = "Pi (Nucleotide Diversity)"), linewidth = 0.55) +
    geom_line(data = win_df, aes(x = Genome_Pos, y = Theta_W_per_bp * 1000, color = "Watterson's Theta"), linewidth = 0.4, alpha = 0.75) +
    scale_x_continuous(breaks = chr_lens$cum_mid, labels = chr_lens$label, expand = c(0.005, 0.005)) +
    scale_color_manual(values = c("Pi (Nucleotide Diversity)" = "#1D4ED8", "Watterson's Theta" = "#D97706")) +
    labs(
      title = paste0(header_prefix, " — Sliding-Window Nucleotide Diversity (π & θw)"),
      x = paste0(sp_cap, " Chromosome"),
      y = "Diversity per kb (× 10⁻³ / bp)",
      color = "Metric"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "top"
    )

  ggsave(file.path("plots", paste0(file_tag, "_sliding_window_pi.png")), p_pi, width = 14, height = 4.5, dpi = 300, bg = "white")
  ggsave(file.path("plots", paste0(file_tag, "_sliding_window_pi.pdf")), p_pi, width = 14, height = 4.5, bg = "white")

  # 2. Sliding-window Tajima's D Track
  if (any(!is.na(win_df$Tajima_D))) {
    p_td <- ggplot() +
      geom_rect(
        data = chr_lens %>% filter(shade),
        aes(xmin = cum_start, xmax = cum_end, ymin = -Inf, ymax = Inf),
        fill = "#F2F4F7", alpha = 0.8
      ) +
      geom_hline(yintercept = 0, color = "#475467", linetype = "solid", linewidth = 0.4) +
      geom_hline(yintercept = c(-2, 2), color = "#DC2626", linetype = "dotted", linewidth = 0.4) +
      geom_line(data = win_df %>% filter(!is.na(Tajima_D)), aes(x = Genome_Pos, y = Tajima_D), color = "#0F766E", linewidth = 0.55) +
      scale_x_continuous(breaks = chr_lens$cum_mid, labels = chr_lens$label, expand = c(0.005, 0.005)) +
      labs(
        title = paste0(header_prefix, " — Sliding-Window Tajima's D"),
        x = paste0(sp_cap, " Chromosome"),
        y = "Tajima's D"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
      )

    ggsave(file.path("plots", paste0(file_tag, "_tajima_d.png")), p_td, width = 14, height = 4.2, dpi = 300, bg = "white")
    ggsave(file.path("plots", paste0(file_tag, "_tajima_d.pdf")), p_td, width = 14, height = 4.2, bg = "white")
  }
}

master_chr <- bind_rows(all_chrom_summaries)
master_win <- bind_rows(all_win_summaries)

if (nrow(master_chr) > 0) {
  subgenome_overview <- master_chr %>%
    group_by(Clade, Subgenome, Prefix) %>%
    summarise(
      Chromosomes       = n(),
      Total_Length_bp   = sum(Length_bp, na.rm = TRUE),
      Segregating_Sites = sum(Segregating_Sites, na.rm = TRUE),
      Genome_Wide_Pi    = sum(Pi_Total_Sum, na.rm = TRUE) / sum(Length_bp, na.rm = TRUE),
      Mean_Chr_Pi       = mean(Pi_per_bp, na.rm = TRUE),
      Mean_Chr_Theta_W  = mean(Theta_W_per_bp, na.rm = TRUE),
      Mean_Chr_Tajima_D = mean(Tajima_D, na.rm = TRUE),
      .groups = "drop"
    )

  write_tsv(master_chr, "cohort_nucleotide_diversity_summary.tsv")
  write_tsv(subgenome_overview, "subgenome_overview.tsv")
  write_tsv(master_win, "all_sliding_windows.tsv")

  # Cross-clade / cross-subgenome comparison plot
  n_clades <- length(unique(master_chr$Clade))
  if (n_clades > 1) {
    p_comp <- ggplot(master_chr, aes(x = Clade, y = Pi_per_bp * 1000, fill = Subgenome)) +
      geom_boxplot(alpha = 0.75, width = 0.55, outlier.shape = NA, position = position_dodge(width = 0.65)) +
      geom_jitter(aes(group = Subgenome), position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.65), size = 1.6, alpha = 0.8, color = "#1F2937") +
      labs(
        title = "Chromosome-Level Nucleotide Diversity (π) Across Clades & Subgenomes",
        x = "Clade",
        y = "Nucleotide Diversity π (× 10⁻³ / bp)",
        fill = "Subgenome"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), legend.position = "top")
  } else {
    cl_lbl <- unique(master_chr$Clade)[1]
    comp_title <- ifelse(tolower(cl_lbl) %in% c("all strains", "all"),
                         "Chromosome-Level Nucleotide Diversity (π) Across Subgenomes",
                         paste0("Chromosome-Level Nucleotide Diversity (π) — ", cl_lbl))
    p_comp <- ggplot(master_chr, aes(x = Subgenome, y = Pi_per_bp * 1000, fill = Subgenome)) +
      geom_boxplot(alpha = 0.65, width = 0.45, outlier.shape = NA) +
      geom_jitter(width = 0.12, size = 2.0, alpha = 0.85, color = "#1F2937") +
      labs(
        title = comp_title,
        x = "Subgenome",
        y = "Nucleotide Diversity π (× 10⁻³ / bp)"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), legend.position = "none")
  }

  ggsave(file.path("plots", "subgenome_pi_comparison.png"), p_comp, width = 9, height = 5.5, dpi = 300, bg = "white")
  ggsave(file.path("plots", "subgenome_pi_comparison.pdf"), p_comp, width = 9, height = 5.5, bg = "white")

  # Build multi-sheet Excel summary via Python
  py_cmd <- "
import csv, zipfile, xml.sax.saxutils as saxutils

def read_tsv(p):
    with open(p, 'r', encoding='utf-8') as f:
        return list(csv.reader(f, delimiter='\\t'))

sheets = {
    'Subgenome_Overview': read_tsv('subgenome_overview.tsv'),
    'Chromosome_Summary': read_tsv('cohort_nucleotide_diversity_summary.tsv'),
    'Sliding_Windows': read_tsv('all_sliding_windows.tsv')
}

with zipfile.ZipFile('cohort_nucleotide_diversity_summary.xlsx', 'w', zipfile.ZIP_DEFLATED) as zf:
    names = list(sheets.keys())
    ovr = '\\n'.join(f'<Override PartName=\"/xl/worksheets/sheet{i+1}.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>' for i in range(len(names)))
    zf.writestr('[Content_Types].xml', f'<?xml version=\"1.0\" encoding=\"UTF-8\"?><Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\"><Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/><Default Extension=\"xml\" ContentType=\"application/xml\"/><Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>{ovr}</Types>')
    zf.writestr('_rels/.rels', '<?xml version=\"1.0\" encoding=\"UTF-8\"?><Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\"><Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"xl/workbook.xml\"/></Relationships>')
    s_xml = ''.join(f'<sheet name=\"{name}\" sheetId=\"{i+1}\" r:id=\"rId{i+1}\"/>' for i, name in enumerate(names))
    zf.writestr('xl/workbook.xml', f'<?xml version=\"1.0\" encoding=\"UTF-8\"?><workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\"><sheets>{s_xml}</sheets></workbook>')
    r_xml = ''.join(f'<Relationship Id=\"rId{i+1}\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet{i+1}.xml\"/>' for i in range(len(names)))
    zf.writestr('xl/_rels/workbook.xml.rels', f'<?xml version=\"1.0\" encoding=\"UTF-8\"?><Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">{r_xml}</Relationships>')
    for idx, name in enumerate(names):
        r_list = []
        for r_i, row in enumerate(sheets[name], start=1):
            c_list = []
            for val in row:
                try:
                    fv = float(val)
                    c_list.append(f'<c><v>{fv}</v></c>')
                except ValueError:
                    c_list.append(f'<c t=\"inlineStr\"><is><t>{saxutils.escape(str(val))}</t></is></c>')
            r_list.append(f'<row r=\"{r_i}\">' + ''.join(c_list) + '</row>')
        zf.writestr(f'xl/worksheets/sheet{idx+1}.xml', '<?xml version=\"1.0\" encoding=\"UTF-8\"?><worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\"><sheetData>' + ''.join(r_list) + '</sheetData></worksheet>')
"
  system2("python3", args = c("-c", shQuote(py_cmd)))
} else {
  writeLines("Clade\tChromosome\tLabel\tLength_bp\tCalled_Haplotypes\tSegregating_Sites\tPi_Total_Sum\tPi_per_bp\tTheta_W_per_bp\tTajima_D", "cohort_nucleotide_diversity_summary.tsv")
}