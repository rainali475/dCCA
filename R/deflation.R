# Remove the previous component's rank-one contribution within each group.
# Marginal covariances stay fixed while cross-covariances are deflated.
deflate_cov <- function(last_var.xy, last_a, last_b, last_rho, array_var.x, array_var.y,
                        p, q, n_samps, diagnostics = NULL) {
  sapply(1:n_samps, function(i) {
    vx <- score_variance(last_a, mat_slice(array_var.x, i), diagnostics, i, "X", "deflation")
    vy <- score_variance(last_b, mat_slice(array_var.y, i), diagnostics, i, "Y", "deflation")
    if (vx == 0 || vy == 0) return(mat_slice(last_var.xy, i))
    var_Xa <- array_var.x[, , i] %*% last_a
    var_Yb <- array_var.y[, , i] %*% last_b
    last_var.xy[, , i] - last_rho[i] *
      tcrossprod(var_Xa, var_Yb) /
      (sqrt(vx) * sqrt(vy))
  }) %>% array(dim = c(p, q, n_samps))
}
