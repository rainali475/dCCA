#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {
// R matrices are column-major. Shared matrices have a zero group stride.
R_xlen_t covariance_stride(const NumericVector& x, int rows, int cols,
                           int groups, const char* label) {
  if (!x.hasAttribute("dim")) stop("%s must be a matrix or array", label);
  IntegerVector dims = x.attr("dim");
  if ((dims.size() != 2 && dims.size() != 3) ||
      dims[0] != rows || dims[1] != cols ||
      (dims.size() == 3 && dims[2] != groups))
    stop("%s has incompatible dimensions", label);
  R_xlen_t stride = static_cast<R_xlen_t>(rows) * cols;
  if (x.size() != stride * (dims.size() == 3 ? groups : 1))
    stop("%s has incompatible length", label);
  return dims.size() == 3 ? stride : 0;
}

int group_count(const NumericVector& xy) {
  if (!xy.hasAttribute("dim")) stop("xy must be a matrix or array");
  IntegerVector dims = xy.attr("dim");
  if (dims.size() != 2 && dims.size() != 3)
    stop("xy must be a matrix or three-dimensional array");
  int groups = dims.size() == 3 ? dims[2] : 1;
  if (groups < 1) stop("At least one group is required");
  return groups;
}

// Keep the R tolerance exactly: eps * dimension * max(abs(S)) * ||a||^2.
// Computing sqrt(vx) * sqrt(vy) later avoids underflow at small scales.
double score(const NumericVector& weights, const NumericVector& covariance,
             R_xlen_t offset, std::vector<double>& product, bool& zero) {
  int n = weights.size();
  std::fill(product.begin(), product.end(), 0.0);
  double maximum = 0, norm = 0, variance = 0;
  for (int j = 0; j < n; ++j) {
    norm += weights[j] * weights[j];
    for (int i = 0; i < n; ++i) {
      double value = covariance[offset + i + static_cast<R_xlen_t>(n) * j];
      product[i] += value * weights[j];
      maximum = std::max(maximum, std::abs(value));
    }
  }
  for (int i = 0; i < n; ++i) variance += weights[i] * product[i];
  zero = variance <= std::numeric_limits<double>::epsilon() * n * maximum * norm;
  return zero ? 0 : variance;
}
}

// Objective, rho, derivative, and zero-variance flags share one group traversal.
// Flags are returned to R, which keeps per-fit diagnostic counts and view names.
// [[Rcpp::export(rng = false)]]
List covariance_kernel(NumericVector a, NumericVector b, NumericVector xy,
                       NumericVector xx, NumericVector yy, NumericVector z,
                       bool fisher = true, bool need_gradient = false) {
  int p = a.size(), q = b.size(), groups = group_count(xy);
  if (p < 1 || q < 1 || z.size() != groups)
    stop("Weights and group counts must match covariance dimensions");
  R_xlen_t xy_stride = covariance_stride(xy, p, q, groups, "xy");
  R_xlen_t xx_stride = covariance_stride(xx, p, p, groups, "xx");
  R_xlen_t yy_stride = covariance_stride(yy, q, q, groups, "yy");
  NumericVector rho(groups), gradient(p);
  LogicalVector zero_x(groups), zero_y(groups);
  std::vector<double> xa(p), yb(q), xyb(p);
  double objective = 0;
  const double limit = 1 - std::sqrt(std::numeric_limits<double>::epsilon());
  for (int g = 0; g < groups; ++g) {
    if (g % 256 == 0) checkUserInterrupt();
    bool zx, zy;
    double vx = score(a, xx, xx_stride * g, xa, zx);
    double vy = score(b, yy, yy_stride * g, yb, zy);
    zero_x[g] = zx; zero_y[g] = zy;
    if (zx || zy) continue;
    std::fill(xyb.begin(), xyb.end(), 0.0);
    for (int j = 0; j < q; ++j)
      for (int i = 0; i < p; ++i)
        xyb[i] += xy[xy_stride * g + i + static_cast<R_xlen_t>(p) * j] * b[j];
    double cross = 0;
    for (int i = 0; i < p; ++i) cross += a[i] * xyb[i];
    double denominator = std::sqrt(vx) * std::sqrt(vy);
    double r = cross / denominator;
    rho[g] = r;
    objective += z[g] * (fisher ? std::atanh(std::max(-limit, std::min(limit, r))) : r);
    if (need_gradient) {
      double weight = z[g];
      // The clipped Fisher objective is flat beyond the clipping boundary.
      if (fisher) weight = std::abs(r) >= limit ? 0 : weight / (1 - r * r);
      for (int i = 0; i < p; ++i)
        gradient[i] += (xyb[i] - cross / vx * xa[i]) * weight / denominator;
    }
  }
  return List::create(_["rho"] = rho, _["objective"] = objective,
                      _["gradient"] = gradient, _["zero_x"] = zero_x,
                      _["zero_y"] = zero_y);
}

// Variance-only path for final diagnostics and shared-variance solver updates.
// [[Rcpp::export(rng = false)]]
List score_variance_kernel(NumericVector weights, NumericVector covariance,
                           int groups) {
  int n = weights.size();
  if (n < 1 || groups < 1) stop("Weights and groups must be nonempty");
  R_xlen_t stride = covariance_stride(covariance, n, n, groups, "covariance");
  NumericVector variance(groups);
  LogicalVector zero(groups);
  std::vector<double> product(n);
  for (int g = 0; g < groups; ++g) {
    if (g % 256 == 0) checkUserInterrupt();
    bool is_zero;
    variance[g] = score(weights, covariance, stride * g, product, is_zero);
    zero[g] = is_zero;
  }
  return List::create(_["variance"] = variance, _["zero"] = zero);
}

// Deflate each valid group's cross-covariance; undefined groups stay unchanged.
// [[Rcpp::export(rng = false)]]
List deflate_covariance_kernel(NumericVector xy, NumericVector a, NumericVector b,
                               NumericVector rho, NumericVector xx, NumericVector yy) {
  int p = a.size(), q = b.size(), groups = group_count(xy);
  if (p < 1 || q < 1 || rho.size() != groups) stop("Invalid deflation dimensions");
  R_xlen_t xy_stride = covariance_stride(xy, p, q, groups, "xy");
  R_xlen_t xx_stride = covariance_stride(xx, p, p, groups, "xx");
  R_xlen_t yy_stride = covariance_stride(yy, q, q, groups, "yy");
  NumericVector output = clone(xy);
  LogicalVector zero_x(groups), zero_y(groups);
  std::vector<double> xa(p), yb(q);
  for (int g = 0; g < groups; ++g) {
    if (g % 256 == 0) checkUserInterrupt();
    bool zx, zy;
    double vx = score(a, xx, xx_stride * g, xa, zx);
    double vy = score(b, yy, yy_stride * g, yb, zy);
    zero_x[g] = zx; zero_y[g] = zy;
    if (zx || zy) continue;
    double denominator = std::sqrt(vx) * std::sqrt(vy);
    for (int j = 0; j < q; ++j)
      for (int i = 0; i < p; ++i) {
        R_xlen_t index = xy_stride * g + i + static_cast<R_xlen_t>(p) * j;
        output[index] -= rho[g] * (xa[i] * yb[j]) / denominator;
      }
  }
  return List::create(_["covariance"] = output, _["zero_x"] = zero_x,
                      _["zero_y"] = zero_y);
}
