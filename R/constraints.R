# Sparsity constraints and convergence ------------------------------------

soft_threshold <- function(x, lambda) {
  sign(x) * pmax(0, abs(x) - lambda)
}

# Relative-change tolerance with an epsilon-squared floor near zero.
converge <- function(x_curr, x_next, epsilon) {
  abs(x_next - x_curr) <= epsilon * (abs(x_curr) + epsilon)
}

# Find a normalized soft-thresholded vector satisfying the L1 bound.
# Tied largest entries use a deterministic sparse fallback, preserving signs.
search_delta <- function(x, c) {
  x <- as.vector(x)
  x <- x / sqrt(sum(x^2))
  if (sum(abs(x)) <= c) {
    return(list(delta = 0, output = x))
  }

  lo <- 0
  hi <- max(abs(x))
  best <- rep(0, length(x))
  best[which.max(abs(x))] <- sign(x[which.max(abs(x))])

  for (i in seq_len(100)) {
    delta <- (lo + hi) / 2
    sx <- soft_threshold(x, delta)
    size <- sqrt(sum(sx^2))
    if (size == 0) {
      hi <- delta
      next
    }

    u <- sx / size
    if (sum(abs(u)) > c) {
      lo <- delta
    } else {
      hi <- delta
      best <- u
    }
    if (abs(sum(abs(u)) - c) < 1e-6) break
  }

  list(delta = hi, output = best)
}
