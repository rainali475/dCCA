# Remove the previous component's rank-one contribution within each group.
# Marginal covariances stay fixed while cross-covariances are deflated.
deflate_cov <- function(last_var.xy, last_a, last_b, last_rho, array_var.x, array_var.y,
                        p, q, n_samps) {
  sapply(1:n_samps, function(i) {
    var_Xa <- array_var.x[, , i] %*% last_a
    var_Yb <- array_var.y[, , i] %*% last_b
    last_var.xy[, , i] - last_rho[i] *
      tcrossprod(var_Xa, var_Yb) /
      (sqrt(sum(last_a * var_Xa) * sum(last_b * var_Yb)))
  }) %>% array(dim = c(p, q, n_samps))
}
