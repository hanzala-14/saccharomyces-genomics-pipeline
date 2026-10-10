#include <Rcpp.h>
#include <fstream>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

static const unsigned char NA_CODE = 255;

static inline int shared_count(unsigned char a1, unsigned char a2,
                               unsigned char b1, unsigned char b2,
                               bool &is_x) {
  const bool a1m = (a1 == NA_CODE);
  const bool a2m = (a2 == NA_CODE);
  const bool b1m = (b1 == NA_CODE);
  const bool b2m = (b2 == NA_CODE);

  if ((a1m && a2m) || (b1m && b2m)) { is_x = true; return -1; }

  unsigned char ak1 = 0, ak2 = 0, bk1 = 0, bk2 = 0;
  int an = 0, bn = 0;

  if (!a1m) { ak1 = a1; an = 1; }
  if (!a2m) { if (an == 0) ak1 = a2; else ak2 = a2; an++; }

  if (!b1m) { bk1 = b1; bn = 1; }
  if (!b2m) { if (bn == 0) bk1 = b2; else bk2 = b2; bn++; }

  int shared = 0;

  if (an == 1 && bn == 1) {
    shared = (ak1 == bk1) ? 1 : 0;
  } else if (an == 1 && bn == 2) {
    shared = (ak1 == bk1 || ak1 == bk2) ? 1 : 0;
  } else if (an == 2 && bn == 1) {
    shared = (bk1 == ak1 || bk1 == ak2) ? 1 : 0;
  } else {
    if (ak1 == bk1 || ak1 == bk2) shared++;
    if (ak2 == bk1 || ak2 == bk2) shared++;

    if (shared == 2) {
      const bool a_homo = (ak1 == ak2);
      const bool b_homo = (bk1 == bk2);
      if (a_homo != b_homo) shared = 1;
    }
  }

  is_x = false;
  return shared;
}

// [[Rcpp::export]]
List computeIBS_from_u8bin_single(std::string bin_allele1,
                                  std::string bin_allele2,
                                  int n_cols,
                                  int64_t n_rows,
                                  int block_rows = 8192,
                                  int n_threads = 1,
                                  int report_blocks = 25) {

  std::ifstream A1(bin_allele1.c_str(), std::ios::in | std::ios::binary);
  std::ifstream A2(bin_allele2.c_str(), std::ios::in | std::ios::binary);
  if (!A1) stop("Cannot open BIN allele1: " + bin_allele1);
  if (!A2) stop("Cannot open BIN allele2: " + bin_allele2);

  IntegerMatrix ibs0(n_cols, n_cols);
  IntegerMatrix ibs1(n_cols, n_cols);
  IntegerMatrix ibs2(n_cols, n_cols);
  IntegerMatrix ibsx(n_cols, n_cols);

  const size_t row_bytes = (size_t)n_cols;

  std::vector<unsigned char> buf1((size_t)block_rows * row_bytes);
  std::vector<unsigned char> buf2((size_t)block_rows * row_bytes);

#ifdef _OPENMP
  if (n_threads < 1) n_threads = 1;
  omp_set_num_threads(n_threads);
#else
  (void)n_threads;
#endif

  int64_t rows_done = 0;
  int block_id = 0;

  while (rows_done < n_rows) {
    const int cur = (int)std::min<int64_t>((int64_t)block_rows, n_rows - rows_done);
    const std::streamsize want = (std::streamsize)cur * (std::streamsize)row_bytes;

    A1.read(reinterpret_cast<char*>(buf1.data()), want);
    A2.read(reinterpret_cast<char*>(buf2.data()), want);

    if (A1.gcount() != want || A2.gcount() != want) {
      stop("Unexpected EOF while reading BIN files (check n_rows and n_cols)");
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < n_cols; i++) {
      // IMPORTANT: start j at i+1 (off-diagonal only).
      // We will set the diagonal after reading all blocks.
      for (int j = i + 1; j < n_cols; j++) {
        int c0 = 0, c1 = 0, c2 = 0, cx = 0;

        for (int k = 0; k < cur; k++) {
          const size_t base = (size_t)k * row_bytes;

          unsigned char a1 = buf1[base + (size_t)i];
          unsigned char a2 = buf2[base + (size_t)i];
          unsigned char b1 = buf1[base + (size_t)j];
          unsigned char b2 = buf2[base + (size_t)j];

          bool is_x = false;
          int sh = shared_count(a1, a2, b1, b2, is_x);
          if (is_x) { cx++; continue; }
          if (sh == 0) c0++;
          else if (sh == 1) c1++;
          else c2++;
        }

        // write symmetric off-diagonal
        ibs0(i,j) += c0; ibs0(j,i) += c0;
        ibs1(i,j) += c1; ibs1(j,i) += c1;
        ibs2(i,j) += c2; ibs2(j,i) += c2;
        ibsx(i,j) += cx; ibsx(j,i) += cx;
      }
    }

    rows_done += cur;
    block_id++;

    if (report_blocks > 0 && (block_id % report_blocks) == 0) {
      Rcpp::Rcout << "Processed rows: " << rows_done << " / " << n_rows << "\n";
    }
  }

  // Force diagonal to represent perfect self-match on ALL considered sites:
  // ibs2[i,i] = n_rows; ibs0/ibs1/ibsx diag = 0
  for (int i = 0; i < n_cols; i++) {
    ibs0(i,i) = 0;
    ibs1(i,i) = 0;
    ibsx(i,i) = 0;
    ibs2(i,i) = (int)n_rows; // ok unless n_rows > 2.1e9; if so, we should switch to numeric/double
  }

  return List::create(
    Named("ibs0") = ibs0,
    Named("ibs1") = ibs1,
    Named("ibs2") = ibs2,
    Named("ibsx") = ibsx
  );
}
