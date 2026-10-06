test_that("constant group features are allowed when the whole data are full rank", {
  t <- seq(-2, 2, length.out = 12)
  X <- list(first = cbind(t, 0), second = cbind(0, t))
  Y <- list(first = cbind(t + sin(t), 0), second = cbind(0, t + cos(t)))
  for (fisher in c(FALSE, TRUE)) for (solver in c("CG", "Nelder-Mead")) {
    fit <- dCCA(X = X, Y = Y, z = c(1, 1), c1 = 1, c2 = 1,
                a_inits = list(c(1, 0)), b_inits = list(c(1, 0)),
                fisher.transform = fisher, solver = solver, max_iter = 5)
    expect_true(fit$diagnostics$zero_variance_encountered)
    expect_true(all(is.finite(fit$result$obj_val_record)))
    expect_true(anyNA(fit$result$rho))
    events <- fit$diagnostics$zero_variance
    expect_true(all(events$count > 0))
    expect_true("second" %in% events$group)
    expect_true("result" %in% events$source)
    if (solver == "CG") expect_true("gradient" %in% events$source)
  }
})

test_that("complementary collinear group covariances have positive definite mean", {
  vx <- array(c(matrix(1, 2, 2), matrix(c(1, -1, -1, 1), 2)), c(2, 2, 2))
  xy <- .4 * vx
  fit <- dCCA(var.xy = xy, var.x = vx, var.y = vx, z = c(1, 1), max_iter = 10)
  expect_true(is.finite(fit$result$obj_val))
  expect_true(fit$diagnostics$zero_variance_encountered)
  expect_error(dCCA(var.xy = xy, var.x = array(rep(matrix(1, 2, 2), 2), c(2, 2, 2)),
                    var.y = vx), "Whole-dataset covariance")
})

test_that("global raw covariance is checked before within-group centering", {
  X <- list(cbind(1:5, 0), cbind(1:5, 1))
  fit <- dCCA(X = X, Y = X, a_inits = list(c(0, 1)),
              b_inits = list(c(0, 1)), fisher.transform = FALSE,
              equal_var.x = TRUE, equal_var.y = TRUE, max_iter = 3)
  expect_true(fit$diagnostics$zero_variance_encountered)
  expect_true(all(is.na(fit$result$rho)))
  expect_true(is.finite(fit$result$obj_val))
  expect_equal(unname(fit$parameters$method), c("CG", "CG"))
  expect_error(dCCA(X = list(cbind(1:5, 1:5)), Y = list(cbind(1:5, 1:5))),
               "Whole-dataset covariance")
})

test_that("indefinite slices are rejected even when their mean is positive definite", {
  vx <- array(c(diag(c(3, -1)), diag(c(1, 3))), c(2, 2, 2))
  expect_error(dCCA(var.xy = array(0, c(2, 2, 2)), var.x = vx), "positive semidefinite")
})

test_that("zero-variance deflation leaves that group's cross-covariance unchanged", {
  vx <- array(c(diag(c(1, 0)), diag(c(0, 1))), c(2, 2, 2))
  xy <- .5 * vx
  d <- dCCA:::new_diagnostics(c("one", "two"))
  deflated <- dCCA:::deflate_cov(xy, c(1, 0), c(1, 0), c(.5, 0), vx, vx,
                                2, 2, 2, diagnostics = d)
  expect_equal(deflated[, , 2], xy[, , 2])
  expect_true(any(vapply(d$events, function(e) e$source == "deflation", logical(1))))
  fit <- dCCA(var.xy = xy, var.x = vx, var.y = vx, ncc = 2,
              c1 = 1, c2 = 1, max_iter = 5)
  expect_length(fit$result, 2)
  expect_true(2 %in% fit$diagnostics$zero_variance$component)
})

test_that("fully positive definite fits return empty diagnostics", {
  fit <- dCCA(var.xy = diag(c(.8, .3)))
  expect_false(fit$diagnostics$zero_variance_encountered)
  expect_equal(nrow(fit$diagnostics$zero_variance), 0L)
  expect_true(all(is.finite(fit$result$rho)))
})

test_that("shared-variance updates skip zero score variances in the other view", {
  vy <- array(c(diag(c(1, 0)), diag(c(0, 1))), c(2, 2, 2))
  fit <- dCCA(var.xy = .3 * vy, var.x = diag(2), var.y = vy,
              equal_var.x = TRUE, fisher.transform = FALSE,
              a_inits = list(c(1, 0)), b_inits = list(c(1, 0)),
              c1 = 1, c2 = 1, max_iter = 5)
  expect_true(is.finite(fit$result$obj_val))
  expect_true("update" %in% fit$diagnostics$zero_variance$source)
  expect_true(anyNA(fit$result$rho))
})

test_that("small covariance scales do not create artificial zero variances", {
  vx <- vy <- array(1e-200, c(1, 1, 1))
  xy <- array(5e-201, c(1, 1, 1))
  expect_equal(dCCA:::obj_fn(1, 1, 1, xy, vx, vy), atanh(.5))
  expect_equal(as.vector(dCCA:::obj_gr(1, 1, 1, xy, vx, vy)), 0)
})
