#' Differential canonical correlation analysis
#'
#' Finds common canonical weights maximizing a weighted sum of Fisher-transformed
#' group correlations by default. Set fisher.transform = FALSE to optimize raw
#' correlations. Sparse fits use alternating optimization and soft thresholding;
#' they are heuristic local solutions, so multiple starts are recommended.
#' @export
#' @importFrom dplyr %>%
#' @importFrom stats optim sd
#' @return A list with `result` (weights `a`, `b`, group correlations `rho`,
#' objective `obj_val`, and iterative diagnostics), `parameters`, and
#' `diagnostics`. For multiple
#' iterative components, `result` is a list of component results. Closed-form
#' weights and correlations are matrices with components in columns.
#' The top-level `diagnostics$zero_variance` table reports component, evaluation
#' source, group, view, and counts of zero-variance divisions avoided across all
#' starts and iterations. Undefined group correlations contribute zero during
#' fitting and are returned as NA for the selected weights.
#' @details Group marginal covariances must be symmetric positive semidefinite,
#' but may be singular. Positive definiteness is checked for the whole raw
#' dataset before within-group centering, or for the mean marginal covariance
#' across groups for covariance-only input. Singular shared covariances use
#' numerical optimization instead of closed-form whitening. Zero variance is
#' detected with a scale-dependent floating-point tolerance. No ridge is added.
#' @examples
#' set.seed(1)
#' X <- list(matrix(rnorm(120), 40, 3), matrix(rnorm(120), 40, 3))
#' Y <- lapply(X, function(x) x + matrix(rnorm(120), 40, 3))
#' fit <- dCCA(X = X, Y = Y, z = c(-1, 1), max_iter = 20)
#' @param z a vector of length n_samps indicating sample groups features.
#' NULL gives equal weights; otherwise supply one finite value per group.
#' Nonconstant weights are standardized; constant weights are used unchanged.
#' @param X Numeric matrix, data frame, or list of paired group matrices with rows being observations and columns
#' being first set of features. Either X and Y or var.xy need to be provided
#' @param Y Numeric matrix, data frame, or list matching X with rows being observations and columns
#' being second set of features. Either X and Y or var.xy need to be provided
#' @param samps Row group labels for matrix inputs. Factors use level order;
#' other labels use first appearance. Omit for list inputs, which use list order.
#' @param equal_var.x Assume shared X covariance. Raw data use pooled within-group
#' covariance; supplied arrays must have identical slices. Default FALSE.
#' @param equal_var.y Assume shared Y covariance, analogous to equal_var.x.
#' @param var.x Positive semidefinite X covariance: a shared p-by-p matrix or
#' a p-by-p-by-group array. Defaults to identity for covariance-only input.
#' @param var.y Positive semidefinite Y covariance: a shared q-by-q matrix or
#' a q-by-q-by-group array. Defaults to identity for covariance-only input.
#' @param var.xy covariance between first and second set of features in each sample.
#' Computed from within-group centered X and Y if not given. Supply a
#' p-by-q matrix for one group or a p-by-q-by-group array.
#' @param c1 L1 bound on unit L2 X weights, at least 1. NULL means unconstrained.
#' @param c2 L1 bound on unit L2 Y weights, at least 1. NULL means unconstrained.
#' @param a_inits List of initial values for a in solving CC if iterative method is used
#' @param b_inits List of initial values for b in solving CC if iterative method is used
#' @param L2_sol_init use non-sparse dCCA solution as initial value for a and b in solving
#' dCCA problem for each CC. Default is FALSE
#' @param epsilon convergence criterion for iterative update of canonical weights
#' @param patience number of iterations to wait before ending iterations when recent
#' iterations do not improve objective value compared to last best result.
#' @param max_iter Maximum iterations, default 1000; NULL also uses 1000.
#' @param solver nested iterative solver for updating one canonical vector given the other.
#' One of "Nelder-Mead" and "CG". Default is "CG"
#' @param solver_eps convergence criterion for nested iterative solver
#' @param ncc number of canonical vectors to find
#' @param fisher.transform Logical scalar. If TRUE, maximize the weighted sum
#' of Fisher-transformed correlations, atanh(rho), using iterative optimization.
#' Correlations are clipped to +/- (1 - sqrt(.Machine$double.eps)) before
#' transformation. Returned rho values and covariance deflation remain on the
#' raw correlation scale. Default TRUE uses the Fisher objective; FALSE selects
#' the original raw-correlation objective.
dCCA <- function(z = NULL, X = NULL, Y = NULL, samps = NULL,
                 equal_var.x = FALSE, equal_var.y = FALSE,
                 var.x = NULL, var.y = NULL, var.xy = NULL,
                 c1 = NULL, c2 = NULL,
                 a_inits = NULL, b_inits = NULL,
                 L2_sol_init = FALSE,
                 epsilon = 1e-6, patience = 10, max_iter = 1000,
                 solver = c("CG", "Nelder-Mead"),
                 solver_eps = 1e-6, ncc = 1, fisher.transform = TRUE) {
  if (!is.logical(fisher.transform) || length(fisher.transform) != 1L ||
      is.na(fisher.transform)) {
    stop("fisher.transform must be a single logical value")
  }
  solver <- match.arg(solver)
  # Normalize user input and validate controls before selecting a solution path.
  prepared <- prepare_input(
    X, Y, samps, z, var.x, var.y, var.xy,
    equal_var.x, equal_var.y
  )
  X <- prepared$X
  Y <- prepared$Y
  samps <- prepared$samps
  z <- prepared$z
  var.x <- prepared$var.x
  var.y <- prepared$var.y
  var.xy <- prepared$var.xy
  p <- dim(var.xy)[1]
  q <- dim(var.xy)[2]
  check_controls(
    p, q, ncc, c1, c2, epsilon, patience, max_iter,
    solver_eps, a_inits, b_inits
  )
  if (is.null(max_iter)) max_iter <- 1000

  p <- ifelse(is.null(X), nrow(var.xy), ncol(X))
  q <- ifelse(is.null(Y), ncol(var.xy), ncol(Y))

  # Group weights are validated and standardized by prepare_input().
  n_samps <- length(z)

  var_args <- get_var(
    X, Y, samps, p, q, n_samps, var.xy, var.x, var.y,
    equal_var.x, equal_var.y
  )
  var.xy <- var_args$var.xy
  var.x <- var_args$var.x
  var.y <- var_args$var.y
  array_var.x <- var_args$array_var.x
  array_var.y <- var_args$array_var.y

  if (is.null(a_inits)) a_inits <- rep(list(rep(1, p)), max(1, length(b_inits)))
  if (is.null(b_inits)) b_inits <- rep(list(rep(1, q)), length(a_inits))

  params <- list(
    var.x = var.x,
    var.y = var.y,
    var.xy = var.xy,
    z = z,
    equal_var.x = equal_var.x,
    equal_var.y = equal_var.y,
    c1 = c1,
    c2 = c2,
    fisher.transform = fisher.transform
  )

  diagnostics <- new_diagnostics(prepared$group_names)
  # A shared but singular within-group covariance can occur even when the
  # whole dataset is full rank. Use numerical updates instead of whitening it.
  if (equal_var.x) equal_var.x <- covariance_is_pd(var.x)
  if (equal_var.y) equal_var.y <- covariance_is_pd(var.y)

  # Shared marginal covariance and no sparsity reduce the weighted objective
  # to a single whitened SVD of the weighted cross-covariance matrix.
  if (!fisher.transform && equal_var.x && equal_var.y && is.null(c1) && is.null(c2)) {
    A <- apply(var.xy, MARGIN = c(1, 2), function(x) sum(x * z))
    CCA_res <- CCA(A, var.x, var.y, ncc)

    params[["method"]] <- "CCA closed form solution"
    return(finish_fit(
      result = list(
        a = CCA_res$a,
        b = CCA_res$b,
        obj_val = CCA_res$rho,
        rho = vapply(seq_len(ncc), function(k) {
          compute_rho(
            CCA_res$a[, k], CCA_res$b[, k], var.xy,
            array_var.x, array_var.y
          )
        }, numeric(n_samps))
      ),
      parameters = params, d = diagnostics,
      var.x = array_var.x, var.y = array_var.y
    ))
  }

  a_method <- dplyr::case_when(
    equal_var.x && !fisher.transform ~ "closed form solution in each iteration",
    !is.null(c1) ~ paste(solver, "with soft threshold"),
    is.null(c1) ~ solver
  )

  b_method <- dplyr::case_when(
    equal_var.y && !fisher.transform ~ "closed form solution in each iteration",
    !is.null(c2) ~ paste(solver, "with soft threshold"),
    is.null(c2) ~ solver
  )

  method_param <- c(a_method, b_method)
  names(method_param) <- c("a", "b")

  params[["epsilon"]] <- epsilon
  if (fisher.transform || !equal_var.x || !equal_var.y) {
    params[["method"]] <- method_param
    params[["solver_eps"]] <- solver_eps
  }

  dcca_res <- solve_single_cc(
    p, q, n_samps, z,
    var.x, var.y, var.xy,
    array_var.x, array_var.y,
    equal_var.x, equal_var.y,
    c1, c2,
    a_inits, b_inits,
    L2_sol_init,
    epsilon, patience, max_iter,
    solver, solver_eps, fisher.transform = fisher.transform,
    diagnostics = diagnostics
  )


  if (ncc == 1) {
    return(finish_fit(dcca_res, params, diagnostics, array_var.x, array_var.y))
  }

  # Extract subsequent components from successively deflated cross-covariances.
  res_list <- list(dcca_res)
  def.var.xy <- var.xy

  for (i in 2:ncc) {
    last_res <- res_list[[length(res_list)]]
    def.var.xy <- deflate_cov(
      last_var.xy = def.var.xy,
      last_a = last_res$a,
      last_b = last_res$b,
      last_rho = last_res$rho,
      array_var.x = array_var.x,
      array_var.y = array_var.y,
      p = p, q = q,
      n_samps = n_samps, diagnostics = diagnostics
    )

    diagnostics$component <- i

    dcca_res <- solve_single_cc(
      p, q, n_samps, z,
      var.x, var.y, def.var.xy,
      array_var.x, array_var.y,
      equal_var.x, equal_var.y,
      c1, c2,
      a_inits, b_inits,
      L2_sol_init,
      epsilon, patience, max_iter,
      solver, solver_eps, fisher.transform = fisher.transform,
      diagnostics = diagnostics
    )

    res_list[[i]] <- dcca_res
  }

  finish_fit(res_list, params, diagnostics, array_var.x, array_var.y)
}
