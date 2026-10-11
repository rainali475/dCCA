# CCA closed form solution
CCA <- function(var.xy, var.x, var.y, ncc) {
  # Whiten both feature sets, solve by SVD, then map weights back to the
  # original coordinates. Euclidean normalization preserves correlations.
  sqrt_varx <- expm::sqrtm(var.x)
  sqrt_vary <- expm::sqrtm(var.y)
  K <- solve(sqrt_varx, t(solve(sqrt_vary, t(var.xy))))
  svd_res <- svd(K, nu = ncc, nv = ncc)
  rho <- svd_res$d[1:ncc]

  a <- solve(sqrt_varx, svd_res$u)
  a <- sweep(a, 2, sqrt(colSums(a^2)), "/")

  b <- solve(sqrt_vary, svd_res$v)
  b <- sweep(b, 2, sqrt(colSums(b^2)), "/")

  list(a = a, b = b, rho = rho)
}

# Equal variance update, solve a with b fixed
eqvar_dCCA_update <- function(var.xy, var.x, var.y, z, curr_b, c1,
                              diagnostics = NULL) {
  sqrt_varx <- expm::sqrtm(var.x)
  A <- apply(var.xy, MARGIN = 3, function(var.xyi) {
    solve(sqrt_varx, var.xyi %*% curr_b)
  }, simplify = FALSE)

  B <- sapply(1:dim(var.xy)[3], function(samp_i) {
    var.yi <- var.y[, , samp_i]
    if (length(dim(var.y)) != 3) var.yi <- var.y
    variance <- score_variance(curr_b, var.yi, diagnostics, samp_i, "Y", "update")
    if (variance == 0) 0 else z[samp_i] / sqrt(variance)
  })

  # Combine the whitened cross-covariance directions using group weights
  # adjusted for the variance of the fixed Y score.
  u <- do.call(cbind, A) %*% B
  curr_a <- solve(sqrt_varx, u)

  # A zero cross-covariance objective has no preferred direction.
  if (sum(curr_a^2) == 0) curr_a <- rep(1, nrow(var.x))
  curr_a <- curr_a / norm(curr_a, type = "2")

  if (!is.null(c1)) {
    search_res <- search_delta(curr_a, c1)
    curr_a <- search_res$output
  }

  return(as.vector(curr_a))
}

# Unequal variance update, solve a with b fixed
dCCA_update <- function(init_a, var.xy, array_var.x, array_var.y, z, curr_b, c1,
                        solver, solver_eps, fisher.transform = TRUE,
                            diagnostics = NULL) {
  pos_fn <- function(a) obj_fn(
    a, curr_b, z, var.xy, array_var.x, array_var.y,
    fisher.transform = fisher.transform, diagnostics = diagnostics
  )
  pos_gr <- function(a) obj_gr(
    a, curr_b, z, var.xy, array_var.x, array_var.y,
    fisher.transform = fisher.transform, diagnostics = diagnostics
  )

  # optim minimizes by default; fnscale = -1 maximizes the correlation sum.
  # Conjugate gradients use the analytic derivative; Nelder-Mead does not.
  if (solver == "CG") {
    pos_res <- optim(
      par = init_a,
      fn = pos_fn,
      gr = pos_gr,
      method = "CG",
      control = list(
        abstol = solver_eps,
        reltol = solver_eps,
        fnscale = -1
      )
    )
  } else if (solver == "Nelder-Mead") {
    pos_res <- optim(
      par = init_a,
      fn = pos_fn,
      method = "Nelder-Mead",
      control = list(
        abstol = solver_eps,
        reltol = solver_eps,
        fnscale = -1
      )
    )
  } else {
    stop("solver must be one of 'CG' or 'Nelder-Mead'")
  }

  # Correlations are scale invariant, so normalize before imposing sparsity.
  if (sum(pos_res$par^2) == 0) pos_res$par <- init_a
  curr_a <- pos_res$par / norm(pos_res$par, type = "2")

  if (!is.null(c1)) {
    search_res <- search_delta(curr_a, c1)
    curr_a <- search_res$output
  }

  return(curr_a)
}

