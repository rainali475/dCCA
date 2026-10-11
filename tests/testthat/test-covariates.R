test_that("closed-form lambda is feasible and maximizes the projected objective", {
  Z <- cbind(age = 1:5, duplicate = 2 * (1:5), other = c(0, 1, 0, 1, 0), constant = 1)
  design <- dCCA:::prepare_covariates(Z, 5)
  expect_equal(design$rank, 2)
  expect_equal(unname(design$dropped_columns), 4L)
  expect_equal(colMeans(design$Z), rep(0, 3), ignore_attr = TRUE)
  expect_equal(apply(design$Z, 2, sd), rep(1, 3), ignore_attr = TRUE)
  rho <- c(-.2, .1, .4, .3, .6)
  for (fisher in c(FALSE, TRUE)) {
    solution <- dCCA:::update_lambda(rho, design, fisher)
    f <- if (fisher) atanh(rho) else rho
    expect_equal(sum(solution$weights^2), 1)
    expect_equal(as.vector(design$Z %*% solution$lambda), solution$weights)
    expect_equal(sum(solution$weights * f), sqrt(sum(crossprod(design$basis, f)^2)))
    expect_equal(unname(solution$lambda[2]), unname(solution$lambda[1]), tolerance = 1e-7)
    gram <- crossprod(design$Z)
    eig <- eigen(gram, symmetric = TRUE)
    keep <- eig$values > 1e-10
    gram_inverse <- sweep(eig$vectors[, keep, drop = FALSE], 2,
                          eig$values[keep], "/") %*% t(eig$vectors[, keep, drop = FALSE])
    numerator <- gram_inverse %*% crossprod(design$Z, f)
    expected <- as.vector(numerator / sqrt(sum(f * (design$Z %*% numerator))))
    expect_equal(unname(solution$lambda), expected, tolerance = 1e-7)
  }
  initial <- dCCA:::update_lambda(rep(0, 5), design, FALSE)
  expect_equal(sum(initial$weights^2), 1)
  expect_identical(dCCA:::update_lambda(rep(0, 5), design, FALSE, initial), initial)
})

test_that("raw group inputs support covariates and zero-variance diagnostics", {
  t <- seq(-2, 2, length.out = 12)
  X <- list(first = cbind(t, 0), second = cbind(0, t))
  Y <- list(first = cbind(t + sin(t), 0), second = cbind(0, t + cos(t)))
  fit <- dCCA(X = X, Y = Y, z = c(-1, 1), c1 = 1, c2 = 1,
              a_inits = list(c(1, 0)), b_inits = list(c(1, 0)), max_iter = 5)
  expect_true(fit$diagnostics$zero_variance_encountered)
  expect_true(anyNA(fit$result$rho))
  expect_equal(sum(fit$result$group_weights^2), 1)
  expect_true(all(is.finite(fit$result$obj_val_record)))
})

test_that("vector and one-column matrix fits are equivalent", {
  xy <- array(c(diag(c(.2, .1)), diag(c(.5, .2)), diag(c(.8, .3))), c(2, 2, 3))
  for (fisher in c(FALSE, TRUE)) {
    args <- list(var.xy = xy, max_iter = 10, fisher.transform = fisher)
    vector <- do.call(dCCA, c(args, list(z = c(-1, 0, 1))))
    matrix <- do.call(dCCA, c(args, list(z = matrix(c(-1, 0, 1), ncol = 1))))
    expect_equal(vector, matrix)
    expect_equal(sum(vector$result$group_weights^2), 1)
    expect_equal(vector$result$obj_val,
                 sum(vector$result$group_weights * if (fisher) atanh(vector$result$rho) else vector$result$rho))
  }
})

test_that("covariates propagate through sparse starts and existing deflation", {
  xy <- array(c(diag(c(.2, .1)), diag(c(.5, .2)), diag(c(.8, .3))), c(2, 2, 3))
  for (solver in c("CG", "Nelder-Mead")) {
    fit <- dCCA(var.xy = xy, z = cbind(1:3, c(0, 1, 0)),
                c1 = 1.2, c2 = 1.2, L2_sol_init = TRUE, ncc = 2,
                solver = solver, max_iter = 5,
                a_inits = list(c(1, .1), c(.1, 1)), b_inits = list(c(1, .1), c(.1, 1)))
    expect_length(fit$result, 2)
    for (res in fit$result) {
      expect_equal(sum(res$group_weights^2), 1)
      expect_equal(res$obj_val, sum(res$group_weights * atanh(res$rho)))
      expect_equal(nrow(res$lambda_record), length(res$obj_val_record))
      expect_lte(sum(abs(res$a)), 1.200001)
      expect_lte(sum(abs(res$b)), 1.200001)
    }
  }
})

test_that("invalid covariates are rejected", {
  xy <- array(0, c(2, 2, 3))
  expect_error(dCCA(var.xy = xy, z = c(1, 1, 1)), "nonconstant")
  expect_error(dCCA(var.xy = xy, z = matrix(1:4, 2)), "one row")
  expect_error(dCCA(var.xy = xy, z = c(1, NA, 3)), "finite numeric")
})
