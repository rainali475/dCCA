# Compute sample covariances from data already centered within each group.
# The third dimension of each returned array follows the levels of samps.
in_samp_cov <- function(X, Y, samps) {
  p <- ncol(X)
  q <- ncol(Y)
  n_samps <- length(levels(samps))

  # Compute empirical covariances
  var.xy <- sapply(levels(samps), function(samp) {
    crossprod(
      X[samps == samp, , drop = FALSE],
      Y[samps == samp, , drop = FALSE]
    ) / (sum(samps == samp) - 1)
  }) %>% array(dim = c(p, q, n_samps))

  # Compute empirical variances
  var.x <- sapply(levels(samps), function(samp) {
    crossprod(X[samps == samp, , drop = FALSE]) / (sum(samps == samp) - 1)
  }) %>% array(dim = c(p, p, n_samps))

  var.y <- sapply(levels(samps), function(samp) {
    crossprod(Y[samps == samp, , drop = FALSE]) / (sum(samps == samp) - 1)
  }) %>% array(dim = c(q, q, n_samps))

  return(list(
    var.xy = var.xy,
    var.x = var.x,
    var.y = var.y
  ))
}

# Organize input variance matrices or compute them from data
get_var <- function(X, Y, samps, p, q, n_samps, var.xy, var.x, var.y,
                    equal_var.x, equal_var.y) {
  if (equal_var.x) {
    if (is.null(var.x)) {
      var.x <- diag(p)
      if (!is.null(X)) var.x <- crossprod(X) / nrow(X)
    } else if (length(dim(var.x)) > 2) {
      var.x <- var.x[, , 1]
    }
  }

  if (equal_var.y) {
    if (is.null(var.y)) {
      var.y <- diag(q)
      if (!is.null(Y)) var.y <- crossprod(Y) / nrow(Y)
    } else if (length(dim(var.y)) > 2) {
      var.y <- var.y[, , 1]
    }
  }

  if (is.null(var.xy) || is.null(var.x) || is.null(var.y)) {
    if (is.null(samps)) samps <- rep(1, nrow(X))
    samps <- factor(samps)
    var_res <- in_samp_cov(X, Y, samps)

    if (is.null(var.xy)) {
      var.xy <- var_res[["var.xy"]]
    }
    if (!equal_var.x && is.null(var.x)) var.x <- var_res[["var.x"]]
    if (!equal_var.y && is.null(var.y)) var.y <- var_res[["var.y"]]
  }

  # Repeat shared covariance matrices so group-wise calculations can use the
  # same array indexing as calculations with distinct group covariances.
  if (length(dim(var.x)) == 2) {
    array_var.x <- rep(var.x, n_samps) %>% array(dim = c(p, p, n_samps))
  } else {
    array_var.x <- var.x
  }

  if (length(dim(var.y)) == 2) {
    array_var.y <- rep(var.y, n_samps) %>% array(dim = c(q, q, n_samps))
  } else {
    array_var.y <- var.y
  }

  return(list(
    var.xy = var.xy,
    var.x = var.x,
    var.y = var.y,
    array_var.x = array_var.x,
    array_var.y = array_var.y
  ))
}
