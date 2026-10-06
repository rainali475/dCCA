# Keep bounded, per-fit counts rather than storing every optimizer evaluation.
new_diagnostics <- function(group_names) {
  d <- new.env(parent = emptyenv())
  d$group_names <- group_names
  d$component <- 1L
  d$swapped <- FALSE
  d$events <- list()
  d
}

record_zero_variance <- function(d, group, view, source) {
  if (is.null(d)) return(invisible(NULL))
  if (d$swapped) view <- if (view == "X") "Y" else "X"
  key <- paste(d$component, source, group, view, sep = ":")
  row <- d$events[[key]]
  if (is.null(row)) {
    row <- data.frame(component = d$component, source = source,
                      group_index = group, group = d$group_names[group],
                      view = view, count = 0L, stringsAsFactors = FALSE)
  }
  row$count <- row$count + 1L
  d$events[[key]] <- row
  invisible(NULL)
}

# Roundoff can leave a tiny positive or negative variance in a null direction.
# The tolerance scales with both the covariance and the weight vector.
score_variance <- function(weights, covariance, d = NULL, group = 1L,
                           view = "X", source = "objective") {
  value <- sum(weights * (covariance %*% weights))
  tolerance <- .Machine$double.eps * nrow(covariance) *
    max(abs(covariance)) * sum(weights^2)
  if (value <= tolerance) {
    record_zero_variance(d, group, view, source)
    return(0)
  }
  value
}

covariance_is_pd <- function(m) {
  values <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  min(values) > .Machine$double.eps * nrow(m) * max(abs(values))
}

# Solver internals use zero contributions for undefined group correlations.
# Expose those correlations as NA in the returned scientific results.
finish_fit <- function(result, parameters, d, var.x, var.y) {
  finalize_component <- function(res, component) {
    a <- as.matrix(res$a)
    b <- as.matrix(res$b)
    for (k in seq_len(ncol(a))) {
      d$component <- component + k - 1L
      for (g in seq_along(d$group_names)) {
        vx <- score_variance(a[, k], mat_slice(var.x, g), d, g, "X", "result")
        vy <- score_variance(b[, k], mat_slice(var.y, g), d, g, "Y", "result")
        if (vx == 0 || vy == 0) {
          if (is.matrix(res$rho)) res$rho[g, k] <- NA_real_ else res$rho[g] <- NA_real_
        }
      }
    }
    res
  }
  d$swapped <- FALSE
  if (!is.null(result$a)) {
    result <- finalize_component(result, 1L)
  } else {
    result <- lapply(seq_along(result), function(i) finalize_component(result[[i]], i))
  }
  events <- if (length(d$events)) do.call(rbind, unname(d$events)) else
    data.frame(component = integer(), source = character(), group_index = integer(),
               group = character(), view = character(), count = integer())
  rownames(events) <- NULL
  list(result = result, parameters = parameters,
       diagnostics = list(zero_variance_encountered = nrow(events) > 0,
                          zero_variance = events,
                          zero_variance_policy = "Zero objective/gradient contribution; returned rho is NA"))
}
