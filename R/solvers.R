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
eqvar_dCCA_update <- function(var.xy, var.x, var.y, z, curr_b, c1) {
  sqrt_varx <- expm::sqrtm(var.x)
  A <- apply(var.xy, MARGIN = 3, function(var.xyi) {
    solve(sqrt_varx, var.xyi %*% curr_b)
  }, simplify = FALSE)

  B <- sapply(1:dim(var.xy)[3], function(samp_i) {
    var.yi <- var.y[, , samp_i]
    if (length(dim(var.y)) != 3) var.yi <- var.y
    z[samp_i] / sqrt(sum(curr_b * (var.yi %*% curr_b)))
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
                        solver, solver_eps, fisher.transform = TRUE) {
  pos_fn <- function(a) obj_fn(
    a, curr_b, z, var.xy, array_var.x, array_var.y,
    fisher.transform = fisher.transform
  )
  pos_gr <- function(a) obj_gr(
    a, curr_b, z, var.xy, array_var.x, array_var.y,
    fisher.transform = fisher.transform
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
  curr_a <- pos_res$par / norm(pos_res$par, type = "2")

  if (!is.null(c1)) {
    search_res <- search_delta(curr_a, c1)
    curr_a <- search_res$output
  }

  return(curr_a)
}

# Solve dCCA iteratively
iteratively_solve_CCA <- function(p, q, n_samps, z,
                                  var.x, var.y, var.xy,
                                  array_var.x, array_var.y,
                                  equal_var.x, equal_var.y,
                                  c1, c2,
                                  a_init, b_init,
                                  epsilon, patience, max_iter,
                                  solver, solver_eps, fisher.transform = TRUE) {
  curr_a <- a_init / norm(a_init, type = "2")
  curr_b <- b_init / norm(b_init, type = "2")

  if (!is.null(c1)) {
    search_res <- search_delta(curr_a, c1)
    curr_a <- search_res$output
  }
  if (!is.null(c2)) {
    search_res <- search_delta(curr_b, c2)
    curr_b <- search_res$output
  }

  # Transposing each group lets the same update routine solve for b.
  var.yx <- apply(var.xy, MARGIN = 3, t) %>% array(dim = c(q, p, n_samps))

  a_record <- curr_a
  b_record <- curr_b
  convergence <- FALSE
  n_iter <- 0
  wait <- 0
  if (is.null(max_iter)) max_iter <- Inf

  obj_val <- obj_fn(
    a = curr_a, b = curr_b, z = z, var.xy = var.xy,
    var.x = array_var.x, var.y = array_var.y,
    fisher.transform = fisher.transform
  )

  # Sparse thresholding can lower the objective. Keep the best weights
  # separately from the complete iteration history.
  last_best_res <- list(
    a = curr_a,
    b = curr_b,
    rho = compute_rho(curr_a, curr_b, var.xy, array_var.x, array_var.y),
    obj_val = obj_val[length(obj_val)],
    n_iter = n_iter,
    convergence = convergence,
    a_record = a_record,
    b_record = b_record,
    obj_val_record = obj_val
  )


  while ((!convergence) && (wait <= patience) && (n_iter < max_iter)) {
    # Fix b, update a
    if (equal_var.x && !fisher.transform) {
      curr_a <- eqvar_dCCA_update(var.xy, var.x, array_var.y, z, curr_b, c1)
    } else {
      curr_a <- dCCA_update(
        init_a = curr_a,
        var.xy = var.xy,
        array_var.x = array_var.x,
        array_var.y = array_var.y,
        z = z,
        curr_b = curr_b,
        c1 = c1,
        solver = solver,
        solver_eps = solver_eps,
        fisher.transform = fisher.transform
      )
    }
    a_record <- rbind(a_record, curr_a)

    # Fix a, update b
    if (equal_var.y && !fisher.transform) {
      curr_b <- eqvar_dCCA_update(var.yx, var.y, array_var.x, z, curr_a, c2)
    } else {
      curr_b <- dCCA_update(
        init_a = curr_b,
        var.xy = var.yx,
        array_var.x = array_var.y,
        array_var.y = array_var.x,
        z = z,
        curr_b = curr_a,
        c1 = c2,
        solver = solver,
        solver_eps = solver_eps,
        fisher.transform = fisher.transform
      )
    }
    b_record <- rbind(b_record, curr_b)

    curr_obj_val <- obj_fn(
      curr_a, curr_b, z, var.xy, array_var.x, array_var.y,
      fisher.transform = fisher.transform
    )
    obj_val <- c(obj_val, curr_obj_val)


    # Stop after sustained deterioration, convergence, or the iteration limit.
    if (curr_obj_val < last_best_res$obj_val) {
      wait <- wait + 1
    } else {
      wait <- 0
    }

    convergence <- converge(obj_val[length(obj_val) - 1], curr_obj_val, epsilon)
    n_iter <- n_iter + 1

    if (obj_val[length(obj_val)] > last_best_res$obj_val) {
      last_best_res <- list(
        a = curr_a,
        b = curr_b,
        rho = compute_rho(curr_a, curr_b, var.xy, array_var.x, array_var.y),
        obj_val = obj_val[length(obj_val)],
        n_iter = n_iter,
        convergence = convergence,
        a_record = a_record,
        b_record = b_record,
        obj_val_record = obj_val
      )
    }
  }

  ret_res <- list(
    a = last_best_res$a,
    b = last_best_res$b,
    rho = last_best_res$rho,
    obj_val = last_best_res$obj_val,
    n_iter = n_iter,
    convergence = convergence,
    a_record = a_record,
    b_record = b_record,
    obj_val_record = obj_val
  )

  return(ret_res)
}
