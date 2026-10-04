# Solve single canonical layer
solve_single_cc <- function(p, q, n_samps, z,
                            var.x, var.y, var.xy,
                            array_var.x, array_var.y,
                            equal_var.x, equal_var.y,
                            c1, c2,
                            a_inits, b_inits,
                            L2_sol_init,
                            epsilon, patience, max_iter,
                            solver, solver_eps, fisher.transform = TRUE) {
  if (L2_sol_init && !(is.null(c1) && is.null(c2))) {
    # Select the best unconstrained fit across starts, then use its weights
    # as the starting point for the sparse fit.


    all_dcca_res <- purrr::map2(a_inits, b_inits, function(a_init, b_init) {
      iteratively_solve_CCA(
        p = p, q = q, n_samps = n_samps, z = z,
        var.x = var.x, var.y = var.y, var.xy = var.xy,
        array_var.x = array_var.x, array_var.y = array_var.y,
        equal_var.x = equal_var.x, equal_var.y = equal_var.y,
        c1 = NULL, c2 = NULL,
        a_init = a_init, b_init = b_init,
        epsilon = epsilon, patience = patience, max_iter = max_iter,
        solver = solver, solver_eps = solver_eps,
        fisher.transform = fisher.transform
      )
    })

    obj_vals <- sapply(all_dcca_res, function(res) res$obj_val)
    dcca_res <- all_dcca_res[[which.max(obj_vals)]]
    l2_a_init <- dcca_res$a
    l2_b_init <- dcca_res$b


    dcca_res <- iteratively_solve_CCA(
      p = p, q = q, n_samps = n_samps, z = z,
      var.x = var.x, var.y = var.y, var.xy = var.xy,
      array_var.x = array_var.x, array_var.y = array_var.y,
      equal_var.x = equal_var.x, equal_var.y = equal_var.y,
      c1 = c1, c2 = c2,
      a_init = l2_a_init, b_init = l2_b_init,
      epsilon = epsilon, patience = patience, max_iter = max_iter,
      solver = solver, solver_eps = solver_eps,
      fisher.transform = fisher.transform
    )
  } else {
    all_dcca_res <- purrr::map2(a_inits, b_inits, function(a_init, b_init) {
      iteratively_solve_CCA(
        p, q, n_samps, z,
        var.x, var.y, var.xy,
        array_var.x, array_var.y,
        equal_var.x, equal_var.y,
        c1, c2,
        a_init, b_init,
        epsilon, patience, max_iter,
        solver, solver_eps, fisher.transform = fisher.transform
      )
    })

    # Alternating optimization can reach different local optima across starts.
    obj_vals <- sapply(all_dcca_res, function(res) res$obj_val)
    dcca_res <- all_dcca_res[[which.max(obj_vals)]]
  }

  dcca_res
}
