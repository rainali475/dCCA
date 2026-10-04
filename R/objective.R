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
compute_rho <- function(a, b, var.xy, var.x, var.y) {
  sapply(1:(dim(var.y)[3]), function(samp_i) {
    sum(a * (mat_slice(var.xy, samp_i) %*% b)) /
      sqrt(sum(a * (mat_slice(var.x, samp_i) %*% a)) *
        sum(b * (mat_slice(var.y, samp_i) %*% b)))
  })
}

# Objective function for maximization
obj_fn <- function(a, b, z, var.xy, var.x, var.y, fisher.transform = TRUE) {
  rho <- compute_rho(a, b, var.xy, var.x, var.y)
  if (fisher.transform) rho <- safe_fisher_transform(rho)
  sum(rho * z)
}

# Gradient with respect to a while b is fixed. Each group contributes the
# derivative of its covariance numerator and its X-variance denominator.
obj_gr <- function(a, b, z, var.xy, var.x, var.y, fisher.transform = TRUE) {
  vapply(seq_along(z), function(samp_i) {
    cov_XYb <- mat_slice(var.xy, samp_i) %*% b
    cov_XaYb <- sum(a * cov_XYb)
    cov_XXa <- mat_slice(var.x, samp_i) %*% a
    var_Xa <- sum(a * cov_XXa)
    var_Yb <- sum(b * (mat_slice(var.y, samp_i) %*% b))
    grad_weight <- z[samp_i]

    if (fisher.transform) {
      rho <- cov_XaYb / sqrt(var_Xa * var_Yb)
      limit <- 1 - sqrt(.Machine$double.eps)
      # The clipped objective is flat outside its bounds. Inside the bounds,
      # apply the chain rule: d atanh(rho) / d rho = 1 / (1 - rho^2).
      grad_weight <- if (abs(rho) >= limit) 0 else grad_weight / (1 - rho^2)
    }

    as.vector((cov_XYb - (cov_XaYb / var_Xa * cov_XXa)) *
      grad_weight / sqrt(var_Xa * var_Yb))
  }, numeric(length(a))) %>%
    matrix(nrow = length(a)) %>%
    rowSums()
}
