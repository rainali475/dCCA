test_that("compiled objectives, correlations and gradients match R", {
  set.seed(904)
  p <- 4; q <- 3; groups <- 5
  xx <- array(0, c(p, p, groups)); yy <- array(0, c(q, q, groups))
  for (g in seq_len(groups)) {
    xx[, , g] <- crossprod(matrix(rnorm(p * p), p)) + diag(p)
    yy[, , g] <- crossprod(matrix(rnorm(q * q), q)) + diag(q)
  }
  xy <- array(rnorm(p * q * groups, sd = .1), c(p, q, groups))
  a <- rnorm(p); b <- rnorm(q); z <- rnorm(groups)
  for (fisher in c(FALSE, TRUE)) {
    for (shared in c(FALSE, TRUE)) {
      vx <- if (shared) xx[, , 1] else xx
      vy <- if (shared) yy[, , 1] else yy
      kernel <- dCCA:::covariance_kernel(a, b, xy, vx, vy, z, fisher, TRUE)
      expect_equal(kernel$rho, dCCA:::compute_rho(a, b, xy, vx, vy, use_cpp = FALSE), tolerance = 1e-12)
      expect_equal(kernel$objective, dCCA:::obj_fn(a, b, z, xy, vx, vy, fisher, use_cpp = FALSE), tolerance = 1e-12)
      expect_equal(kernel$gradient, dCCA:::obj_gr(a, b, z, xy, vx, vy, fisher, use_cpp = FALSE), tolerance = 1e-12)
      numeric_gradient <- vapply(seq_along(a), function(i) {
        step <- numeric(p); step[i] <- 1e-6
        (dCCA:::obj_fn(a + step, b, z, xy, vx, vy, fisher) -
           dCCA:::obj_fn(a - step, b, z, xy, vx, vy, fisher)) / 2e-6
      }, numeric(1))
      expect_equal(kernel$gradient, numeric_gradient, tolerance = 1e-7)
      expect_false(any(kernel$zero_x | kernel$zero_y))
    }
  }
})

test_that("compiled clipping, tiny scales and single-feature views match R", {
  for (r in c(-1, -1 + 1e-9, -.5, 0, .5, 1 - 1e-9, 1)) {
    for (fisher in c(FALSE, TRUE)) {
      xx <- yy <- array(1e-200, c(1, 1, 1))
      xy <- array(r * 1e-200, c(1, 1, 1))
      expect_equal(dCCA:::obj_fn(1, 1, 1, xy, xx, yy, fisher),
                   dCCA:::obj_fn(1, 1, 1, xy, xx, yy, fisher, use_cpp = FALSE))
      expect_equal(dCCA:::obj_gr(1, 1, 1, xy, xx, yy, fisher),
                   dCCA:::obj_gr(1, 1, 1, xy, xx, yy, fisher, use_cpp = FALSE))
    }
  }
  xy <- array(c(.9, .1, .1, .2), c(2, 2, 1))
  a <- b <- c(1, 0)
  for (r in c(1, -1, 1 - 1e-9)) {
    xy[1, 1, 1] <- r
    expect_equal(dCCA:::obj_gr(a, b, 1, xy, diag(2), diag(2)), c(0, 0))
  }
})

test_that("zero flags and diagnostic counts retain component and view semantics", {
  xx <- array(c(diag(c(1, 0)), diag(c(0, 1))), c(2, 2, 2))
  xy <- .4 * xx
  a <- b <- c(1, 0)
  diagnostics <- lapply(c(TRUE, FALSE), function(cpp) {
    d <- dCCA:::new_diagnostics(c("first", "second"), cpp)
    d$component <- 2L
    d$swapped <- TRUE
    for (i in 1:2) {
      dCCA:::compute_rho(a, b, xy, xx, xx, d)
      dCCA:::obj_fn(a, b, c(1, -1), xy, xx, xx, diagnostics = d)
      dCCA:::obj_gr(a, b, c(1, -1), xy, xx, xx, diagnostics = d)
    }
    d
  })
  expect_equal(diagnostics[[1]]$events, diagnostics[[2]]$events)
  expect_equal(diagnostics[[1]]$events[[1]]$view, "Y")
  kernel <- dCCA:::covariance_kernel(a, b, xy, xx, xx, c(1, -1), TRUE, TRUE)
  expect_identical(kernel$zero_x, c(FALSE, TRUE))
  expect_identical(kernel$zero_y, c(FALSE, TRUE))
  expect_equal(kernel$rho, c(.4, 0))
})

test_that("compiled variances and deflation match R including null groups", {
  xx <- array(c(diag(c(1, 0)), diag(c(0, 1))), c(2, 2, 2))
  xy <- .4 * xx
  a <- b <- c(1, 0)
  expect_equal(dCCA:::score_variance_groups(a, xx, 2),
               dCCA:::score_variance_groups(a, xx, 2, FALSE))
  cpp_d <- dCCA:::new_diagnostics(c("first", "second"))
  r_d <- dCCA:::new_diagnostics(c("first", "second"), FALSE)
  cpp <- dCCA:::deflate_cov(xy, a, b, c(.4, 0), xx, xx, 2, 2, 2, cpp_d)
  ref <- dCCA:::deflate_cov(xy, a, b, c(.4, 0), xx, xx, 2, 2, 2, r_d)
  expect_equal(cpp, ref)
  expect_equal(cpp_d$events, r_d$events)
  expect_equal(xy[, , 1], diag(c(.4, 0))) # Native code must not mutate inputs.
  expect_equal(cpp[, , 2], xy[, , 2])
})

test_that("both backends fit raw and Fisher covariate models", {
  xx <- array(rep(diag(2), 3), c(2, 2, 3))
  xy <- array(c(diag(c(.6, .2)), diag(c(.3, .1)), diag(c(.1, .4))), c(2, 2, 3))
  for (fisher in c(FALSE, TRUE)) for (solver in c("CG", "Nelder-Mead")) {
    for (z in list(NULL, c(-1, 0, 1), cbind(c(-1, 0, 1), c(1, -2, 1)))) {
      args <- list(var.xy = xy, var.x = xx, var.y = xx, z = z,
                   fisher.transform = fisher, solver = solver, max_iter = 3,
                   ncc = 2, c1 = 1, c2 = 1, L2_sol_init = TRUE,
                   a_inits = list(c(1, 0)), b_inits = list(c(1, 0)))
      cpp <- do.call(dCCA, args)
      ref <- do.call(dCCA, c(args, list(use_cpp = FALSE)))
      expect_true(cpp$parameters$use_cpp)
      expect_false(ref$parameters$use_cpp)
      expect_equal(cpp$result, ref$result, tolerance = 1e-7)
      expect_equal(cpp$diagnostics, ref$diagnostics)
    }
  }
})

test_that("public backend flag is validated and native dimensions are checked", {
  expect_identical(formals(dCCA)$use_cpp, TRUE)
  for (bad in list(NA, NULL, 1, c(TRUE, FALSE)))
    expect_error(dCCA(var.xy = diag(2), use_cpp = bad), "single logical")
  expect_error(dCCA:::covariance_kernel(1, 1, array(0, c(2, 1, 1)), matrix(1), matrix(1), 1), "dimensions")
  expect_error(dCCA:::covariance_kernel(1, 1, matrix(1), matrix(1), matrix(1), c(1, 1)), "group counts")
  expect_error(dCCA:::score_variance_kernel(1, matrix(1), 0), "nonempty")
})

test_that("dense raw-data fits and shared closed-form fits support either backend", {
  set.seed(173)
  X <- lapply(seq_len(3), function(g) matrix(rnorm(180), 60, 3))
  Y <- lapply(X, function(x) .3 * x[, 1:2] + matrix(rnorm(120), 60, 2))
  for (fisher in c(FALSE, TRUE)) {
    cpp <- dCCA(X = X, Y = Y, fisher.transform = fisher, max_iter = 5)
    ref <- dCCA(X = X, Y = Y, fisher.transform = fisher, max_iter = 5, use_cpp = FALSE)
    expect_equal(cpp$result$obj_val, ref$result$obj_val, tolerance = 1e-7)
    expect_equal(cpp$result$a, ref$result$a, tolerance = 1e-6)
    expect_equal(cpp$result$b, ref$result$b, tolerance = 1e-6)
    expect_equal(cpp$result$rho, ref$result$rho, tolerance = 1e-6)
  }
  args <- list(var.xy = diag(c(.8, .3)), equal_var.x = TRUE,
               equal_var.y = TRUE, fisher.transform = FALSE, ncc = 2)
  expect_equal(do.call(dCCA, args)$result,
               do.call(dCCA, c(args, list(use_cpp = FALSE)))$result)
})
