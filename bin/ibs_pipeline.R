run_ibs_pipeline <- function(config) {
  defaults <- list(
    allele1_tsv = NULL,
    allele2_tsv = NULL,
    single_id = NA_character_,
    input_dir = NULL,
    samples_txt = NULL,
    code_dir = getwd(),
    output_root = file.path(getwd(), "ibs_outputs"),
    write_final_combined = TRUE,
    log_file = file.path(getwd(), sprintf("ibs_run_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))),
    threads = 16L,
    block_rows = 8192L,
    report_blocks = 5L,
    overwrite = FALSE,
    distance = list(set_den0_to = NA_real_)
  )

  config <- modifyList(defaults, config)

  dir.create(dirname(config$log_file), recursive = TRUE, showWarnings = FALSE)
  log_con <- file(config$log_file, open = "wt")
  sink(log_con, type = "output")
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(log_con)
  }, add = TRUE)

  cat("Log file:", config$log_file, "\n")
  cat("Start time:", as.character(Sys.time()), "\n\n")

  cat("Compiling C++ code...\n")
  Sys.setenv("PKG_CXXFLAGS" = "-O3 -march=native -fopenmp")
  Sys.setenv("PKG_LIBS"     = "-fopenmp")
  load_ibs_cpp(config$code_dir)

  single_ready <- !is.na(config$allele1_tsv) && nzchar(config$allele1_tsv) &&
                  !is.na(config$allele2_tsv) && nzchar(config$allele2_tsv)

  if (single_ready) {
    id <- config$single_id
    if (is.na(id) || !nzchar(id)) id <- infer_id_from_pair(config$allele1_tsv, config$allele2_tsv)
    datasets <- list(list(id = id, allele1_tsv = config$allele1_tsv, allele2_tsv = config$allele2_tsv))
    cat("AUTO selected: SINGLE\n")
  } else {
    if (is.na(config$input_dir) || !nzchar(config$input_dir)) {
      stop("AUTO needs either (allele1_tsv + allele2_tsv) OR input_dir.")
    }
    datasets <- discover_pairs_in_dir(config$input_dir)
    cat("AUTO selected: BATCH\n")
  }

  cat("Datasets:", length(datasets), "\n")
  for (d in datasets) cat("  -", d$id, "\n")

  dir.create(config$output_root, recursive = TRUE, showWarnings = FALSE)

  dataset_results <- vector("list", length(datasets))
  names(dataset_results) <- vapply(datasets, `[[`, character(1), "id")

  for (idx in seq_along(datasets)) {
    d <- datasets[[idx]]
    id <- d$id

    cat("\n============================================================\n")
    cat("DATASET", idx, "/", length(datasets), ":", id, "\n")

    out_dir <- file.path(config$output_root, id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    bin_dir <- file.path(out_dir, "bin_cache")
    dir.create(bin_dir, recursive = TRUE, showWarnings = FALSE)

    bin_a1 <- file.path(bin_dir, paste0(id, "_allele1.u8bin"))
    bin_a2 <- file.path(bin_dir, paste0(id, "_allele2.u8bin"))

    ibs0_tsv <- file.path(out_dir, "IBS0.tsv")
    ibs1_tsv <- file.path(out_dir, "IBS1.tsv")
    ibs2_tsv <- file.path(out_dir, "IBS2.tsv")
    ibsx_tsv <- file.path(out_dir, "IBSx.tsv")
    dist_tsv <- file.path(out_dir, "distance_matrix.tsv")

    if (!config$overwrite &&
        file.exists(ibs0_tsv) && file.exists(ibs1_tsv) &&
        file.exists(ibs2_tsv) && file.exists(ibsx_tsv) &&
        file.exists(dist_tsv)) {
      cat("Skipping (outputs already exist):", id, "\n")
      dataset_results[[id]] <- list(
        id = id, skipped = TRUE, out_dir = out_dir,
        ibs0_tsv = ibs0_tsv, ibs1_tsv = ibs1_tsv,
        ibs2_tsv = ibs2_tsv, ibsx_tsv = ibsx_tsv,
        dist_tsv = dist_tsv
      )
      next
    }

    n_cols1 <- infer_n_cols(d$allele1_tsv)
    n_cols2 <- infer_n_cols(d$allele2_tsv)
    if (n_cols1 != n_cols2) stop("Column mismatch between TSVs in dataset: ", id)

    n_rows1 <- count_lines_fast(d$allele1_tsv)
    n_rows2 <- count_lines_fast(d$allele2_tsv)
    if (n_rows1 != n_rows2) stop("Row mismatch between TSVs in dataset: ", id)

    n_cols <- n_cols1
    n_rows <- n_rows1

    cat("allele1:", d$allele1_tsv, "\n")
    cat("allele2:", d$allele2_tsv, "\n")
    cat("n_cols (samples):", n_cols, "\n")
    cat("n_rows (sites):  ", n_rows, "\n")

    sample_names <- read_sample_names_optional(config$samples_txt, n_cols)

    if (check_bin_ok(bin_a1, n_rows, n_cols)) {
      cat("BIN OK (skip):", bin_a1, "\n")
    } else {
      cat("Converting allele1 TSV -> BIN...\n")
      tsv_to_u8bin(d$allele1_tsv, bin_a1, n_cols = n_cols, report_every = 50000)
    }

    if (check_bin_ok(bin_a2, n_rows, n_cols)) {
      cat("BIN OK (skip):", bin_a2, "\n")
    } else {
      cat("Converting allele2 TSV -> BIN...\n")
      tsv_to_u8bin(d$allele2_tsv, bin_a2, n_cols = n_cols, report_every = 50000)
    }

    cat("Computing IBS with threads=", config$threads,
        ", block_rows=", config$block_rows, " ...\n", sep = "")

    t0 <- proc.time()[["elapsed"]]
    res <- computeIBS_from_u8bin_single(
      bin_a1, bin_a2,
      n_cols = n_cols,
      n_rows = n_rows,
      block_rows = config$block_rows,
      n_threads = config$threads,
      report_blocks = config$report_blocks
    )
    t1 <- proc.time()[["elapsed"]]
    cat("IBS finished. Elapsed:", fmt_time(t1 - t0), "\n")

    dimnames(res$ibs0) <- list(sample_names, sample_names)
    dimnames(res$ibs1) <- list(sample_names, sample_names)
    dimnames(res$ibs2) <- list(sample_names, sample_names)
    dimnames(res$ibsx) <- list(sample_names, sample_names)

    write_matrix_tsv(res$ibs0, ibs0_tsv)
    write_matrix_tsv(res$ibs1, ibs1_tsv)
    write_matrix_tsv(res$ibs2, ibs2_tsv)
    write_matrix_tsv(res$ibsx, ibsx_tsv)

    dist_m <- compute_distance_from_ibs(
      res$ibs0, res$ibs1, res$ibs2,
      sample_names = sample_names,
      den0_value = config$distance$set_den0_to
    )

    dist_m <- force_all_ibsx_samples_to_distance_one(
      dist_m,
      res$ibs0, res$ibs1, res$ibs2, res$ibsx
    )

    write_matrix_tsv(dist_m, dist_tsv)

    dataset_results[[id]] <- list(
      id = id, skipped = FALSE, out_dir = out_dir,
      n_rows = n_rows, n_cols = n_cols,
      ibs0_tsv = ibs0_tsv, ibs1_tsv = ibs1_tsv,
      ibs2_tsv = ibs2_tsv, ibsx_tsv = ibsx_tsv,
      dist_tsv = dist_tsv
    )
  }

  combined_result <- NULL

  if (isTRUE(config$write_final_combined) && length(datasets) > 1) {
    cat("\n============================================================\n")
    cat("Combining IBS count tables by summation\n")

    valid_ids <- names(dataset_results)
    valid_ids <- valid_ids[vapply(valid_ids, function(id) file.exists(dataset_results[[id]]$ibs0_tsv), logical(1))]
    if (length(valid_ids) == 0) stop("No IBS tables found to combine.")

    first <- read_matrix_tsv(dataset_results[[valid_ids[1]]]$ibs0_tsv)
    combined_names <- rownames(first)
    if (is.null(combined_names)) stop("First IBS table has no rownames; cannot align.")

    n <- nrow(first)

    IBS0_all <- matrix(0, n, n, dimnames = list(combined_names, combined_names))
    IBS1_all <- matrix(0, n, n, dimnames = list(combined_names, combined_names))
    IBS2_all <- matrix(0, n, n, dimnames = list(combined_names, combined_names))
    IBSx_all <- matrix(0, n, n, dimnames = list(combined_names, combined_names))

    for (id in valid_ids) {
      info <- dataset_results[[id]]
      ibs0 <- read_matrix_tsv(info$ibs0_tsv)
      ibs1 <- read_matrix_tsv(info$ibs1_tsv)
      ibs2 <- read_matrix_tsv(info$ibs2_tsv)
      ibsx <- read_matrix_tsv(info$ibsx_tsv)

      if (is.null(rownames(ibs0)) || is.null(colnames(ibs0))) {
        stop("IBS matrix for dataset ", id, " has no dimnames; cannot align.")
      }
      if (!setequal(rownames(ibs0), combined_names)) {
        stop("Sample name set mismatch in dataset ", id, " (cannot combine IBS tables).")
      }

      ibs0 <- ibs0[combined_names, combined_names]
      ibs1 <- ibs1[combined_names, combined_names]
      ibs2 <- ibs2[combined_names, combined_names]
      ibsx <- ibsx[combined_names, combined_names]

      IBS0_all <- IBS0_all + ibs0
      IBS1_all <- IBS1_all + ibs1
      IBS2_all <- IBS2_all + ibs2
      IBSx_all <- IBSx_all + ibsx

      cat("  added:", id, "\n")
    }

    final_dist_m <- compute_distance_from_ibs(
      IBS0_all, IBS1_all, IBS2_all,
      sample_names = combined_names,
      den0_value = config$distance$set_den0_to
    )

    final_dist_m <- force_all_ibsx_samples_to_distance_one(
      final_dist_m,
      IBS0_all, IBS1_all, IBS2_all, IBSx_all
    )

    final_ibs0_tsv <- file.path(config$output_root, "IBS0_all.tsv")
    final_ibs1_tsv <- file.path(config$output_root, "IBS1_all.tsv")
    final_ibs2_tsv <- file.path(config$output_root, "IBS2_all.tsv")
    final_ibsx_tsv <- file.path(config$output_root, "IBSx_all.tsv")
    final_dist_tsv <- file.path(config$output_root, "FINAL_combined_distance_matrix.tsv")

    write_matrix_tsv(IBS0_all, final_ibs0_tsv)
    write_matrix_tsv(IBS1_all, final_ibs1_tsv)
    write_matrix_tsv(IBS2_all, final_ibs2_tsv)
    write_matrix_tsv(IBSx_all, final_ibsx_tsv)
    write_matrix_tsv(final_dist_m, final_dist_tsv)

    cat("Wrote combined IBS0:", final_ibs0_tsv, "\n")
    cat("Wrote combined IBS1:", final_ibs1_tsv, "\n")
    cat("Wrote combined IBS2:", final_ibs2_tsv, "\n")
    cat("Wrote combined IBSx:", final_ibsx_tsv, "\n")
    cat("Wrote combined distance:", final_dist_tsv, "\n")

    combined_result <- list(
      ibs0_tsv = final_ibs0_tsv,
      ibs1_tsv = final_ibs1_tsv,
      ibs2_tsv = final_ibs2_tsv,
      ibsx_tsv = final_ibsx_tsv,
      distance_tsv = final_dist_tsv
    )
  } else {
    cat("\nSkipping combined IBS output (only 1 dataset or disabled).\n")
  }

  cat("\nAll done.\nEnd time:", as.character(Sys.time()), "\n")
  cat("Outputs in:", config$output_root, "\n")
  cat("Log file:", config$log_file, "\n")

  list(
    config = config,
    datasets = dataset_results,
    combined = combined_result,
    output_root = config$output_root,
    log_file = config$log_file
  )
}
