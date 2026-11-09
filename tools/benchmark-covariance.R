# Run from the repository root after compiling the package, for example with
# pkgbuild::compile_dll(debug = FALSE). Timings are illustrative, not test gates.
pkgload::load_all(quiet = TRUE, compile = FALSE)
set.seed(904)
p <- 30L
q <- 20L
groups <- 40L
iterations <- 300L
xy <- array(rnorm(p * q * groups, sd = .01), c(p, q, groups))
xx <- array(rep(diag(p), groups), c(p, p, groups))
yy <- array(rep(diag(q), groups), c(q, q, groups))
a <- rnorm(p)
b <- rnorm(q)
z <- rnorm(groups)

measure <- function(use_cpp) {
  invisible(gc())
  unname(system.time(for (i in seq_len(iterations)) {
    dCCA:::obj_fn(a, b, z, xy, xx, yy, use_cpp = use_cpp)
    dCCA:::obj_gr(a, b, z, xy, xx, yy, use_cpp = use_cpp)
  })["elapsed"])
}

# Warm both paths before measuring repeated objective + gradient evaluations.
invisible(dCCA:::obj_fn(a, b, z, xy, xx, yy))
invisible(dCCA:::obj_gr(a, b, z, xy, xx, yy, use_cpp = FALSE))
r_time <- median(replicate(3, measure(FALSE)))
cpp_time <- median(replicate(3, measure(TRUE)))
print(data.frame(backend = c("R", "C++"), elapsed_seconds = c(r_time, cpp_time)))
cat(sprintf("%.1fx covariance evaluation speedup (%d pairs; p=%d, q=%d, groups=%d)\n",
            r_time / cpp_time, iterations, p, q, groups))
