# Backend choice belongs to one fit, never to global options or mutable state.
covariance_use_cpp <- function(diagnostics) {
  is.null(diagnostics) || !identical(diagnostics$use_cpp, FALSE)
}

# Record only exceptional groups in R. Keep X/Y interleaving and swapped-view
# semantics identical to score_variance(), including counts for repeated calls.
record_kernel_zeros <- function(result, diagnostics, source = "objective") {
  if (is.null(diagnostics)) return(invisible(NULL))
  for (g in which(result$zero_x | result$zero_y)) {
    if (result$zero_x[g]) record_zero_variance(diagnostics, g, "X", source)
    if (result$zero_y[g]) record_zero_variance(diagnostics, g, "Y", source)
  }
  invisible(NULL)
}

score_variance_groups <- function(weights, covariance, groups, use_cpp = TRUE) {
  if (use_cpp) return(score_variance_kernel(weights, covariance, groups)$variance)
  vapply(seq_len(groups), function(g) {
    score_variance(weights, mat_slice(covariance, g))
  }, numeric(1))
}
