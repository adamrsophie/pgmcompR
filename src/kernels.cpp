#include <Rcpp.h>
using namespace Rcpp;

// Zero-based index into a factor for every cell of the full grid over a
// union scope. `ucard` gives the union cardinalities, `pos` the zero-based
// positions of the factor's variables in the union, `fcard` the factor's
// cardinalities (same order as `pos`). Uses an odometer instead of division.
// [[Rcpp::export(rng = false)]]
NumericVector cpp_expand_index(IntegerVector ucard, IntegerVector pos, IntegerVector fcard) {
  int nu = ucard.size();
  double total = 1.0;
  for (int j = 0; j < nu; ++j) total *= ucard[j];
  R_xlen_t n = (R_xlen_t) total;
  std::vector<double> stride(nu, 0.0);
  double s = 1.0;
  for (int k = 0; k < pos.size(); ++k) {
    stride[pos[k]] = s;
    s *= fcard[k];
  }
  std::vector<int> digit(nu, 0);
  NumericVector out(n);
  double idx = 0.0;
  for (R_xlen_t c = 0; c < n; ++c) {
    out[c] = idx;
    for (int j = 0; j < nu; ++j) {
      digit[j]++;
      idx += stride[j];
      if (digit[j] < ucard[j]) break;
      idx -= stride[j] * digit[j];
      digit[j] = 0;
    }
  }
  return out;
}

// Sum (or log-sum-exp) over the trailing block of a factor already permuted
// so that kept variables come first: values form an nk x nd matrix.
// [[Rcpp::export(rng = false)]]
NumericVector cpp_row_reduce(NumericVector values, double nk, bool logscale) {
  R_xlen_t nrow = (R_xlen_t) nk;
  R_xlen_t ncol = values.size() / nrow;
  NumericVector out(nrow);
  if (!logscale) {
    for (R_xlen_t j = 0; j < ncol; ++j) {
      const double *col = &values[j * nrow];
      for (R_xlen_t i = 0; i < nrow; ++i) out[i] += col[i];
    }
    return out;
  }
  std::vector<double> mx(nrow, R_NegInf);
  for (R_xlen_t j = 0; j < ncol; ++j)
    for (R_xlen_t i = 0; i < nrow; ++i)
      if (values[j * nrow + i] > mx[i]) mx[i] = values[j * nrow + i];
  for (R_xlen_t i = 0; i < nrow; ++i) {
    if (mx[i] == R_NegInf) { out[i] = R_NegInf; continue; }
    double acc = 0.0;
    for (R_xlen_t j = 0; j < ncol; ++j) acc += std::exp(values[j * nrow + i] - mx[i]);
    out[i] = mx[i] + std::log(acc);
  }
  return out;
}

// Marginalize variables out of a factor without materializing a permutation:
// for each cell, add its value into the output cell given by `outidx`.
// [[Rcpp::export(rng = false)]]
NumericVector cpp_scatter_sum(NumericVector values, NumericVector outidx, double nout) {
  NumericVector out((R_xlen_t) nout);
  R_xlen_t n = values.size();
  for (R_xlen_t c = 0; c < n; ++c) out[(R_xlen_t) outidx[c]] += values[c];
  return out;
}
