suppressPackageStartupMessages({
  library(Rcpp)
})

fmt_time <- function(seconds) {
  seconds <- as.numeric(seconds)
  if (!is.finite(seconds) || seconds < 0) return("NA")
  h <- floor(seconds / 3600); seconds <- seconds - 3600 * h
  m <- floor(seconds / 60);   seconds <- seconds - 60 * m
  s <- floor(seconds)
  if (h > 0) sprintf("%dh %02dm %02ds", h, m, s) else sprintf("%dm %02ds", m, s)
}

infer_n_cols <- function(path) {
  first <- readLines(path, n = 1L)
  if (length(first) == 0L) stop("Empty file: ", path)
  length(strsplit(first, "\t", fixed = TRUE)[[1]])
}

count_lines_fast <- function(path, buf_size = 4 * 1024^2) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)

  n <- 0L
  last_byte <- as.raw(0)
  saw_any <- FALSE

  repeat {
    x <- readBin(con, what = "raw", n = buf_size)
    if (length(x) == 0L) break
    saw_any <- TRUE
    n <- n + sum(x == as.raw(0x0A))
    last_byte <- x[length(x)]
  }

  if (saw_any && last_byte != as.raw(0x0A)) n <- n + 1L
  n
}

read_sample_names_optional <- function(samples_path, n_cols) {
  if (is.na(samples_path) || !nzchar(samples_path)) {
    return(paste0("S", seq_len(n_cols)))
  }
  nm <- readLines(samples_path, warn = FALSE)
  nm <- trimws(nm)
  nm <- nm[nzchar(nm)]
  if (length(nm) != n_cols) {
    stop("Sample name count mismatch: ", samples_path, " has ", length(nm),
         " lines but TSV has n_cols=", n_cols)
  }
  if (anyDuplicated(nm)) {
    dup <- unique(nm[duplicated(nm)])
    stop("Duplicate sample names in ", samples_path, ", e.g.: ", paste(head(dup, 10), collapse = ", "))
  }
  nm
}

check_bin_ok <- function(bin_path, n_rows, n_cols) {
  expected <- as.numeric(n_rows) * as.numeric(n_cols)
  if (!file.exists(bin_path)) return(FALSE)
  info <- file.info(bin_path)
  if (is.na(info$size)) return(FALSE)
  info$size == expected
}

compute_distance_from_ibs <- function(ibs0, ibs1, ibs2, sample_names, den0_value = NA_real_) {
  den <- ibs0 + ibs1 + ibs2
  sim <- matrix(NA_real_, nrow = nrow(ibs0), ncol = ncol(ibs0),
                dimnames = list(sample_names, sample_names))
  ok <- den > 0
  sim[ok] <- (2 * ibs2[ok] + ibs1[ok]) / (2 * den[ok])
  dist_m <- 1 - sim
  if (!is.na(den0_value)) dist_m[!ok] <- den0_value
  diag(dist_m) <- 0
  dist_m
}

force_all_ibsx_samples_to_distance_one <- function(dist_m, ibs0_m, ibs1_m, ibs2_m, ibsx_m) {
  stopifnot(
    nrow(dist_m) == nrow(ibs0_m),
    nrow(dist_m) == nrow(ibs1_m),
    nrow(dist_m) == nrow(ibs2_m),
    nrow(dist_m) == nrow(ibsx_m),
    ncol(dist_m) == ncol(ibs0_m),
    ncol(dist_m) == ncol(ibs1_m),
    ncol(dist_m) == ncol(ibs2_m),
    ncol(dist_m) == ncol(ibsx_m),
    identical(rownames(dist_m), rownames(ibs0_m)),
    identical(rownames(dist_m), rownames(ibs1_m)),
    identical(rownames(dist_m), rownames(ibs2_m)),
    identical(rownames(dist_m), rownames(ibsx_m))
  )

  n <- nrow(dist_m)

  for (i in seq_len(n)) {
    offdiag <- setdiff(seq_len(n), i)
    if (length(offdiag) == 0) next

    all_ibsx_only <- all(
      ibs0_m[i, offdiag] == 0 &
      ibs1_m[i, offdiag] == 0 &
      ibs2_m[i, offdiag] == 0 &
      ibsx_m[i, offdiag] > 0,
      na.rm = TRUE
    )

    if (all_ibsx_only) {
      dist_m[i, offdiag] <- 1
      dist_m[offdiag, i] <- 1
      dist_m[i, i] <- 0
    }
  }

  dist_m
}

infer_id_from_pair <- function(a1, a2) {
  b1 <- sub("\\.tsv$", "", basename(a1))
  b2 <- sub("\\.tsv$", "", basename(a2))
  b1 <- sub("^allele1[_-]?", "", b1)
  b2 <- sub("^allele2[_-]?", "", b2)
  if (identical(b1, b2)) return(b1)
  b1
}

discover_pairs_in_dir <- function(input_dir) {
  files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (length(files) == 0) stop("No .tsv files found in: ", input_dir)

  base <- basename(files)

  m1 <- regexec("^(.+)\\.allele1\\.tsv$", base, perl = TRUE)
  p1 <- regmatches(base, m1)
  ok1 <- lengths(p1) == 2
  df1 <- data.frame(
    file = files[ok1],
    id = vapply(p1[ok1], `[[`, character(1), 2),
    stringsAsFactors = FALSE
  )

  m2 <- regexec("^(.+)\\.allele2\\.tsv$", base, perl = TRUE)
  p2 <- regmatches(base, m2)
  ok2 <- lengths(p2) == 2
  df2 <- data.frame(
    file = files[ok2],
    id = vapply(p2[ok2], `[[`, character(1), 2),
    stringsAsFactors = FALSE
  )

  if (nrow(df1) == 0) stop("No *.allele1.tsv files found in ", input_dir)
  if (nrow(df2) == 0) stop("No *.allele2.tsv files found in ", input_dir)

  ids <- sort(intersect(df1$id, df2$id))
  if (length(ids) == 0) stop("No matching allele1/allele2 pairs found.")

  datasets <- list()
  for (id in ids) {
    a1 <- df1$file[df1$id == id]
    a2 <- df2$file[df2$id == id]
    if (length(a1) != 1 || length(a2) != 1) {
      message("Skipping id '", id, "' (need exactly one allele1 and one allele2).")
      next
    }

    datasets[[length(datasets) + 1]] <- list(
      id = id,
      allele1_tsv = a1,
      allele2_tsv = a2
    )
  }

  if (length(datasets) == 0) stop("No valid dataset pairs after filtering.")
  datasets
}

load_ibs_cpp <- function(code_dir) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  setwd(code_dir)
  Rcpp::sourceCpp("tsv_to_u8bin.cpp")
  Rcpp::sourceCpp("ibs_from_bin_single_omp.cpp")
}

write_matrix_tsv <- function(mat, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    mat,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
  )
}

read_matrix_tsv <- function(path) {
  as.matrix(read.table(
    file = path,
    sep = "\t",
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  ))
}
