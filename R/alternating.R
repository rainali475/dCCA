# Solve dCCA iteratively
iteratively_solve_CCA <- function(p, q, n_samps, z,
                                  var.x, var.y, var.xy,
                                  array_var.x, array_var.y,
                                  equal_var.x, equal_var.y,
                                  c1, c2,
                                  a_init, b_init,
                                  epsilon, patience, max_iter,
                                  solver, solver_eps, fisher.transform = TRUE,
                            diagnostics = NULL, covariates = NULL) {
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

  lambda_state <- NULL
  lambda_record <- NULL
  if (!is.null(covariates)) {
    lambda_state <- update_lambda(
      compute_rho(curr_a, curr_b, var.xy, array_var.x, array_var.y, diagnostics),
      covariates, fisher.transform
    )
    z <- lambda_state$weights
    lambda_record <- matrix(lambda_state$lambda, nrow = 1)
  }

  obj_val <- obj_fn(
    a = curr_a, b = curr_b, z = z, var.xy = var.xy,
    var.x = array_var.x, var.y = array_var.y,
    fisher.transform = fisher.transform, diagnostics = diagnostics
  )

  # Sparse thresholding can lower the objective. Keep the best weights
  # separately from the complete iteration history.
  last_best_res <- list(
    a = curr_a,
    b = curr_b,
    lambda_state = lambda_state,
    rho = compute_rho(curr_a, curr_b, var.xy, array_var.x, array_var.y, diagnostics),
    obj_val = obj_val[length(obj_val)],
    n_iter = n_iter,
    convergence = convergence,
    a_record = a_record,
    b_record = b_record,
    obj_val_record = obj_val
  )


  while ((!convergence) && (wait <= patience) && (n_iter < max_iter)) {
    # Fix b, update a
    if (!is.null(diagnostics)) diagnostics$swapped <- FALSE
    if (equal_var.x && !fisher.transform) {
      curr_a <- eqvar_dCCA_update(var.xy, var.x, array_var.y, z, curr_b, c1, diagnostics)
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
        fisher.transform = fisher.transform, diagnostics = diagnostics
      )
    }
    a_record <- rbind(a_record, curr_a)

    # Fix a, update b
    if (!is.null(diagnostics)) diagnostics$swapped <- TRUE
    if (equal_var.y && !fisher.transform) {
      curr_b <- eqvar_dCCA_update(var.yx, var.y, array_var.x, z, curr_a, c2, diagnostics)
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
        fisher.transform = fisher.transform, diagnostics = diagnostics
      )
    }
    b_record <- rbind(b_record, curr_b)
    if (!is.null(diagnostics)) diagnostics$swapped <- FALSE

    # Update lambda after a and b; the next sweep holds these group weights fixed.
    if (!is.null(covariates)) {
      lambda_state <- update_lambda(
        compute_rho(curr_a, curr_b, var.xy, array_var.x, array_var.y, diagnostics),
        covariates, fisher.transform, lambda_state
      )
      z <- lambda_state$weights
      lambda_record <- rbind(lambda_record, lambda_state$lambda)
    }

    curr_obj_val <- obj_fn(
      curr_a, curr_b, z, var.xy, array_var.x, array_var.y,
      fisher.transform = fisher.transform, diagnostics = diagnostics
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
        lambda_state = lambda_state,
        rho = compute_rho(curr_a, curr_b, var.xy, array_var.x, array_var.y, diagnostics),
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

  if (!is.null(covariates)) {
    ret_res$lambda <- last_best_res$lambda_state$lambda
    ret_res$group_weights <- last_best_res$lambda_state$weights
    ret_res$lambda_record <- lambda_record
    colnames(ret_res$lambda_record) <- colnames(covariates$Z)
  }

  return(ret_res)
}
