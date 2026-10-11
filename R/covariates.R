# Standardize group-level covariates and retain their identifiable column space.
prepare_covariates <- function(z, groups) {
  z <- as.matrix(z)
  if (!is.numeric(z) || nrow(z) != groups || ncol(z) < 1L ||
      any(!is.finite(z))) {
    stop("Covariates require a finite numeric matrix with one row per group")
  }
  center <- colMeans(z)
  scale <- apply(z, 2, stats::sd)
  keep <- is.finite(scale) & scale > 0
  if (!any(keep)) stop("Covariates require at least one nonconstant column")
  Z <- sweep(sweep(z[, keep, drop = FALSE], 2, center[keep], "-"),
             2, scale[keep], "/")
  decomposition <- svd(Z)
  # A relative singular-value tolerance avoids inverting redundant directions.
  retained <- decomposition$d > max(dim(Z)) * .Machine$double.eps *
    max(decomposition$d)
  U <- decomposition$u[, retained, drop = FALSE]
  V <- decomposition$v[, retained, drop = FALSE]
  singular <- decomposition$d[retained]
  # (Z'Z)^+ Z' = V D^-1 U'. Applying this factorization avoids squaring
  # the condition number by explicitly constructing the Gram matrix.
  inverse <- sweep(V, 2, singular, "/") %*% t(U)
  list(Z = Z, inverse = inverse, basis = U, rank = length(singular),
       center = center, scale = scale, retained_columns = which(keep),
       dropped_columns = which(!keep), column_names = colnames(z))
}

# Closed-form maximizer subject to lambda' Z' Z lambda = 1. The minimum-norm
# representative is used when coefficients are not uniquely identifiable.
update_lambda <- function(rho, covariates, fisher.transform, previous = NULL) {
  target <- if (fisher.transform) safe_fisher_transform(rho) else rho
  coordinates <- as.vector(crossprod(covariates$basis, target))
  size <- sqrt(sum(coordinates^2))
  tolerance <- max(dim(covariates$Z)) * .Machine$double.eps *
    sqrt(sum(target^2))
  if (size <= tolerance) {
    # Every feasible lambda has objective zero for an orthogonal target.
    if (!is.null(previous)) return(previous)
    weights <- covariates$basis[, 1]
    lambda <- as.vector(covariates$inverse %*% weights)
  } else {
    weights <- as.vector(covariates$basis %*% (coordinates / size))
    lambda <- as.vector(covariates$inverse %*% (target / size))
  }
  names(lambda) <- colnames(covariates$Z)
  list(lambda = lambda, weights = weights)
}
