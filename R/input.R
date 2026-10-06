# Validate inputs once and center observations within their sample groups.
prepare_input <- function(X, Y, samps, z, var.x, var.y, var.xy,
                          equal_var.x, equal_var.y) {
  global.x <- global.y <- NULL
  group_names <- if (is.list(X) && !is.data.frame(X)) names(X) else NULL
  for (flag in list(equal_var.x, equal_var.y)) {
    if (!is.logical(flag) || length(flag) != 1 || is.na(flag)) stop("equal_var flags must be logical scalars")
  }

  data_matrix <- function(x) {
    x <- as.matrix(x)
    if (!is.numeric(x) || any(!is.finite(x)) || nrow(x) < 2 || ncol(x) < 1) {
      stop("Data matrices require finite numeric values and at least two rows")
    }
    x
  }

  # Stack paired list elements into a common row layout before computing
  # covariances. Group labels retain the original list order.
  if (!is.null(X) || !is.null(Y)) {
    if (is.null(X) || is.null(Y)) stop("Supply both X and Y")
    if (is.list(X) && !is.data.frame(X)) {
      if (!is.list(Y) || is.data.frame(Y) || !length(X) || length(X) != length(Y)) {
        stop("X and Y must be paired lists of equal nonzero length")
      }
      if (!is.null(samps)) stop("Do not supply samps with list inputs")
      if (!identical(names(X), names(Y))) stop("List names must match in order")
      X <- lapply(X, data_matrix)
      Y <- lapply(Y, data_matrix)

      if (length(unique(vapply(X, ncol, integer(1)))) != 1 ||
        length(unique(vapply(Y, ncol, integer(1)))) != 1) {
        stop("Feature counts must match across groups")
      }
      if (!identical(vapply(X, nrow, integer(1)), vapply(Y, nrow, integer(1)))) {
        stop("Paired matrices must have equal row counts")
      }
      samps <- rep(seq_along(X), vapply(X, nrow, integer(1)))
      X <- do.call(rbind, X)
      Y <- do.call(rbind, Y)
    } else {
      X <- data_matrix(X)
      Y <- data_matrix(Y)
    }
    if (nrow(X) != nrow(Y)) stop("X and Y must have equal row counts")
    if (is.null(samps)) samps <- rep(1, nrow(X))
    if (length(samps) != nrow(X) || anyNA(samps)) stop("samps must label every row")

    # Preserve factor level order, otherwise use order of first appearance.
    samps <- if (is.factor(samps)) droplevels(samps) else factor(samps, levels = unique(samps))
    if (any(table(samps) < 2)) stop("Each group needs at least two observations")

    # Check the whole dataset before within-group centering. Between-group
    # differences can make the combined data full rank despite singular groups.
    global.x <- stats::cov(X)
    global.y <- stats::cov(Y)

    # Center each group separately so cross-products estimate within-group
    # covariance rather than differences between group means.
    for (group in levels(samps)) {
      rows <- samps == group
      X[rows, ] <- scale(X[rows, , drop = FALSE], scale = FALSE)
      Y[rows, ] <- scale(Y[rows, , drop = FALSE], scale = FALSE)
    }
    empirical <- in_samp_cov(X, Y, samps)
    if (is.null(var.xy)) var.xy <- empirical$var.xy

    # Pooled estimates lose one degree of freedom per centered group.
    if (is.null(var.x)) {
      var.x <- if (equal_var.x) {
        crossprod(X) / (nrow(X) - nlevels(samps))
      } else {
        empirical$var.x
      }
    }
    if (is.null(var.y)) {
      var.y <- if (equal_var.y) {
        crossprod(Y) / (nrow(Y) - nlevels(samps))
      } else {
        empirical$var.y
      }
    }
  }

  # Internally, cross-covariances always have a third (group) dimension,
  # including when the user supplies a matrix for a single group.
  if (is.null(var.xy)) stop("Supply X and Y or var.xy")
  if (length(dim(var.xy)) == 2) var.xy <- array(var.xy, c(dim(var.xy), 1))
  if (!is.numeric(var.xy) || length(dim(var.xy)) != 3 || any(dim(var.xy) < 1) || any(!is.finite(var.xy))) {
    stop("var.xy must be a finite numeric matrix or three-dimensional array")
  }
  p <- dim(var.xy)[1]
  q <- dim(var.xy)[2]
  groups <- dim(var.xy)[3]
  if (!is.null(samps) && nlevels(samps) != groups) stop("Covariance group counts must match samps")

  # Shared marginal matrices may be supplied directly or as identical slices.
  # Group matrices may be singular, but not genuinely indefinite.
  # For covariance-only input, the mean of the slices represents the whole data.
  validate_cov <- function(v, size, equal, global) {
    if (is.null(v)) v <- diag(size)
    dims <- dim(v)
    if (!is.numeric(v) || !(length(dims) %in% c(2, 3)) ||
      !identical(as.integer(dims[1:2]), c(size, size)) ||
      (length(dims) == 3 && dims[3] != groups) || any(!is.finite(v))) {
      stop("Variance matrices require matching feature and group dimensions")
    }
    if (length(dims) == 3 && equal) {
      first <- matrix(v[, , 1], size, size)
      for (i in seq_len(groups)) {
        if (!isTRUE(all.equal(matrix(v[, , i], size, size), first))) {
          stop("Equal variance arrays must contain identical matrices")
        }
      }
      v <- first
    }
    matrices <- lapply(seq_len(if (length(dim(v)) == 3) groups else 1), function(i) {
      m <- if (length(dim(v)) == 3) matrix(v[, , i], size, size) else v
      if (!isTRUE(all.equal(m, t(m), check.attributes = FALSE))) {
        stop("Variance matrices must be symmetric")
      }
      values <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
      if (min(values) < -sqrt(.Machine$double.eps) * max(abs(values))) {
        stop("Group variance matrices must be positive semidefinite")
      }
      m
    })
    if (is.null(global)) global <- Reduce(`+`, matrices) / length(matrices)
    if (!covariance_is_pd(global)) {
      stop("Whole-dataset covariance must be positive definite; remove globally redundant features or supply regularized covariances")
    }
    v
  }

  var.x <- validate_cov(var.x, p, equal_var.x, global.x)
  var.y <- validate_cov(var.y, q, equal_var.y, global.y)
  if (is.null(group_names)) group_names <- if (!is.null(samps)) levels(samps) else
    dimnames(var.xy)[[3]]
  if (is.null(group_names)) group_names <- as.character(seq_len(groups))

  # Standardize varying weights only: constant weights have zero standard
  # deviation and instead represent an equally weighted correlation objective.
  if (is.null(z)) z <- rep(1, groups)
  if (!is.numeric(z) || length(z) != groups || any(!is.finite(z))) stop("z requires one finite numeric weight per group")
  if (groups > 1 && stats::sd(z) > 0) z <- (z - mean(z)) / stats::sd(z)
  if (all(z == 0)) stop("z must contain a nonzero weight")

  list(
    X = X, Y = Y, samps = samps, z = z,
    var.x = var.x, var.y = var.y, var.xy = var.xy, group_names = group_names
  )
}

# Reject invalid constraints and starting vectors before entering a solver.
check_controls <- function(p, q, ncc, c1, c2, epsilon, patience,
                           max_iter, solver_eps, a_inits, b_inits) {
  scalar <- function(x, lower, integer = FALSE) {
    is.numeric(x) && length(x) == 1 && is.finite(x) && x >= lower && (!integer || x == floor(x))
  }
  if (!scalar(ncc, 1, TRUE) || ncc > min(p, q)) stop("ncc must be between 1 and min(p, q)")
  # Any unit L2 vector has an L1 norm of at least one.
  for (bound in list(c1, c2)) {
    if (!is.null(bound) && !scalar(bound, 1)) stop("L1 bounds must be finite and at least 1")
  }
  if (!scalar(epsilon, .Machine$double.eps) || !scalar(solver_eps, .Machine$double.eps) ||
    !scalar(patience, 0, TRUE) || (!is.null(max_iter) && !scalar(max_iter, 1, TRUE))) {
    stop("Invalid iteration controls")
  }
  check_starts <- function(starts, size) {
    if (is.null(starts)) {
      return(invisible(NULL))
    }
    if (!is.list(starts) || !length(starts) || any(!vapply(
      starts, function(x) {
        is.numeric(x) && length(x) == size && all(is.finite(x)) && sum(x^2) > 0
      },
      logical(1)
    ))) {
      stop("Initial values must be lists of finite nonzero vectors of matching length")
    }
  }
  check_starts(a_inits, p)
  check_starts(b_inits, q)
  if (!is.null(a_inits) && !is.null(b_inits) && length(a_inits) != length(b_inits)) stop("Initial lists must have equal lengths")
}
