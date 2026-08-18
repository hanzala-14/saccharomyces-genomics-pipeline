#!/usr/bin/env Rscript

# ==============================================================================
# CLI Argument Parsing
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: run_ibs.R <input_dir> <samples_txt> <threads> <block_rows> <report_blocks> <bin_dir>")
}

INPUT_DIR     <- args[1]
SAMPLES_TXT   <- args[2]
THREADS       <- as.integer(args[3])
BLOCK_ROWS    <- as.integer(args[4])
REPORT_BLOCKS <- as.integer(args[5])
BIN_DIR       <- args[6]

# ==============================================================================
# Load Dependencies from the Nextflow bin/ Directory
# ==============================================================================
source(file.path(BIN_DIR, "ibs_utils.R"))
source(file.path(BIN_DIR, "ibs_pipeline.R"))

# ==============================================================================
# Internal Config Construction
# ==============================================================================
CONFIG <- list(
  allele1_tsv = "",
  allele2_tsv = "",
  single_id = "",
  input_dir = INPUT_DIR,
  samples_txt = SAMPLES_TXT,
  code_dir = BIN_DIR,
  output_root = "ibs_outputs",
  write_final_combined = TRUE,
  log_file = "ibs_outputs/ibs_run.log",
  threads = THREADS,
  block_rows = BLOCK_ROWS,
  report_blocks = REPORT_BLOCKS,
  overwrite = TRUE,
  distance = list(
    set_den0_to = NA_real_
  )
)

# ==============================================================================
# Execution
# ==============================================================================
result <- run_ibs_pipeline(CONFIG)
print(paste("Results saved to:", result$output_root))