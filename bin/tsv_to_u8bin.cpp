#include <Rcpp.h>
#include <fstream>
#include <string>
#include <vector>
#include <cctype>
#include <cstdint>

using namespace Rcpp;

static const unsigned char NA_CODE = 255;

static inline unsigned char tok_to_code(const char* p, size_t len) {
  if (len == 0) return 254;

  if (len == 1 && p[0] == '.') return NA_CODE;

  unsigned int v = 0;
  for (size_t i = 0; i < len; ++i) {
    unsigned char c = (unsigned char)p[i];
    if (!std::isdigit(c)) return 254;
    v = v * 10u + (unsigned int)(c - '0');
    if (v > 254u) return 254;
  }

  return (unsigned char)v;
}

// [[Rcpp::export]]

List tsv_to_u8bin(std::string tsv_path,
                  std::string bin_path,
                  int n_cols,
                  int skip_lines = 0,
                  int report_every = 50000) {

  std::ifstream in(tsv_path.c_str(), std::ios::in | std::ios::binary);
  if (!in) stop("Cannot open input TSV: " + tsv_path);

  std::ofstream out(bin_path.c_str(), std::ios::out | std::ios::binary);
  if (!out) stop("Cannot open output BIN: " + bin_path);

  std::string line;
  for (int i = 0; i < skip_lines; ++i) {
    if (!std::getline(in, line)) stop("skip_lines exceeds file length");
  }

  std::vector<unsigned char> row((size_t)n_cols);

  long long rows_written = 0;
  while (std::getline(in, line)) {
    int col = 0;
    size_t start = 0;

    while (true) {
      size_t tab = line.find('\t', start);
      size_t end = (tab == std::string::npos) ? line.size() : tab;

      if (col >= n_cols) {
        stop("More columns than n_cols at row " + std::to_string(rows_written + 1));
      }

      unsigned char code = tok_to_code(line.data() + start, end - start);
      if (code == 254) {
        stop("Unexpected token at row " + std::to_string(rows_written + 1) +
             ", col " + std::to_string(col + 1) + " (token: \'" +
             line.substr(start, end - start) + "\')");
      }
      row[(size_t)col] = code;

      col++;
      if (tab == std::string::npos) break;
      start = tab + 1;
    }

    if (col != n_cols) {
      stop("Wrong number of columns at row " + std::to_string(rows_written + 1) +
           ": got " + std::to_string(col) + ", expected " + std::to_string(n_cols));
    }

    out.write(reinterpret_cast<const char*>(row.data()), (std::streamsize)row.size());
    rows_written++;

    if (report_every > 0 && (rows_written % report_every) == 0) {
      Rcpp::Rcout << "Converted rows: " << rows_written << "\n";
    }
  }

  out.flush();
  return List::create(
    Named("rows_written") = rows_written,
    Named("cols") = n_cols,
    Named("bin_path") = bin_path
  );
}
