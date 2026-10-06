# Objective functions and gradients ---------------------------------------

# Retain matrix dimensions for one-feature views and accept shared matrices.
mat_slice <- function(A, i) {
  d <- dim(A)
  if (length(d) == 2L) return(A)
  matrix(A[, , i, drop = FALSE], nrow = d[1], ncol = d[2])
}

# Fisher's transform is finite only inside (-1, 1). Clipping also limits the
# influence of nearly perfect correlations on the optimization objective.
clip_rho <- function(rho, eps = sqrt(.Machine$double.eps)) {
  pmin(pmax(rho, -1 + eps), 1 - eps)
}

safe_fisher_transform <- function(rho, eps = sqrt(.Machine$double.eps)) {
  atanh(clip_rho(rho, eps))
}

# Compute a single rho for each sample group based on common canonical vectors
compute_rho <- function(a, b, var.xy, var.x, var.y, diagnostics = NULL) {
  sapply(1:(dim(var.y)[3]), function(samp_i) {
    vx <- score_variance(a, mat_slice(var.x, samp_i), diagnostics, samp_i, "X")
    vy <- score_variance(b, mat_slice(var.y, samp_i), diagnostics, samp_i, "Y")
    if (vx == 0 || vy == 0) return(0)
    sum(a * (mat_slice(var.xy, samp_i) %*% b)) /
      (sqrt(vx) * sqrt(vy))
  })
}

# Objective function for maximization
obj_fn <- function(a, b, z, var.xy, var.x, var.y, fisher.transform = TRUE,
                   diagnostics = NULL) {
  rho <- compute_rho(a, b, var.xy, var.x, var.y, diagnostics)
  if (fisher.transform) rho <- safe_fisher_transform(rho)
  sum(rho * z)
}

# Gradient with respect to a while b is fixed. Each group contributes the
# derivative of its covariance numerator and its X-variance denominator.
obj_gr <- function(a, b, z, var.xy, var.x, var.y, fisher.transform = TRUE,
                   diagnostics = NULL) {
  vapply(seq_along(z), function(samp_i) {
    cov_XYb <- mat_slice(var.xy, samp_i) %*% b
    cov_XaYb <- sum(a * cov_XYb)
    cov_XXa <- mat_slice(var.x, samp_i) %*% a
    var_Xa <- score_variance(a, mat_slice(var.x, samp_i), diagnostics,
                            samp_i, "X", "gradient")
    var_Yb <- score_variance(b, mat_slice(var.y, samp_i), diagnostics,
                            samp_i, "Y", "gradient")
    if (var_Xa == 0 || var_Yb == 0) return(rep(0, length(a)))
    grad_weight <- z[samp_i]

    if (fisher.transform) {
      rho <- cov_XaYb / (sqrt(var_Xa) * sqrt(var_Yb))
      limit <- 1 - sqrt(.Machine$double.eps)
      # The clipped objective is flat outside its bounds. Inside the bounds,
      # apply the chain rule: d atanh(rho) / d rho = 1 / (1 - rho^2).
      grad_weight <- if (abs(rho) >= limit) 0 else grad_weight / (1 - rho^2)
    }

    as.vector((cov_XYb - (cov_XaYb / var_Xa * cov_XXa)) *
      grad_weight / (sqrt(var_Xa) * sqrt(var_Yb)))
  }, numeric(length(a))) %>%
    matrix(nrow = length(a)) %>%
    rowSums()
}
