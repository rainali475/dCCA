test_that("closed form agrees with known canonical correlations", {
  fit <- dCCA(var.xy = diag(c(.8, .3)), var.x = diag(2), var.y = diag(2),
              equal_var.x = TRUE, equal_var.y = TRUE, ncc = 2,
              fisher.transform = FALSE)
  expect_equal(fit$result$obj_val, c(.8, .3))
  expect_equal(dim(fit$result$a), c(2L, 2L))
  expect_equal(as.vector(fit$result$rho), c(.8, .3))
})

test_that("list and matrix inputs agree and centering removes offsets", {
  set.seed(12)
  X <- list(matrix(rnorm(160), 80, 2), matrix(rnorm(160), 80, 2))
  Y <- lapply(X, function(x) x + matrix(rnorm(160), 80, 2))
  fit <- dCCA(X = X, Y = Y, z = c(-1, 1), max_iter = 20)
  stacked <- dCCA(X = do.call(rbind, X), Y = do.call(rbind, Y),
                  samps = rep(c("b", "a"), each = 80), z = c(-1, 1), max_iter = 20)
  shifted <- dCCA(X = lapply(X, function(x) x + 100), Y = Y,
                  z = c(-1, 1), max_iter = 20)
  expect_equal(fit$result$obj_val, stacked$result$obj_val)
  expect_equal(fit$parameters$var.xy, shifted$parameters$var.xy, tolerance = 1e-10)
  expect_true(all(is.finite(fit$result$rho)))
  expect_lte(fit$result$n_iter, 20)
})

test_that("constant weights and sparse tied vectors stay finite", {
  fit <- dCCA(var.xy = array(c(.8, .4), c(1, 1, 2)), z = c(1, 1),
              equal_var.x = TRUE, equal_var.y = TRUE)
  expect_equal(as.vector(fit$result$rho), c(.8, .4))
  sparse <- dCCA:::search_delta(c(-1, 1), 1)$output
  expect_equal(sum(abs(sparse)), 1)
  expect_equal(sum(sparse^2), 1)
  expect_lte(sparse[1], 0)
})

test_that("invalid inputs fail before optimization", {
  expect_error(dCCA(), "Supply")
  expect_error(dCCA(var.xy = diag(2), z = c(1, 2)), "one finite")
  expect_error(dCCA(var.xy = diag(2), c1 = .5), "L1")
  expect_error(dCCA(var.xy = diag(2), ncc = 3), "ncc")
  expect_error(dCCA(var.xy = diag(2), var.x = matrix(0, 2, 2)), "positive definite")
  expect_error(dCCA(X = list(matrix(1, 2, 2)), Y = list(matrix(1, 3, 2))), "row counts")
})

test_that("analytic gradient matches finite differences", {
  vx <- array(diag(c(2, 1)), c(2, 2, 1))
  vy <- array(diag(2), c(2, 2, 1))
  xy <- array(diag(c(.6, .2)), c(2, 2, 1))
  a <- c(.7, -.2); b <- c(.4, .8)
  numeric_gradient <- vapply(seq_along(a), function(i) {
    step <- c(0, 0); step[i] <- 1e-6
    (dCCA:::obj_fn(a + step, b, 1, xy, vx, vy, FALSE) -
       dCCA:::obj_fn(a - step, b, 1, xy, vx, vy, FALSE)) / 2e-6
  }, numeric(1))
  expect_equal(as.vector(dCCA:::obj_gr(a, b, 1, xy, vx, vy, FALSE)), numeric_gradient,
               tolerance = 1e-6)
})

test_that("sparse fits respect bounds and deflation extracts another component", {
  xy <- array(diag(c(.8, .3)), c(2, 2, 1))
  fit <- dCCA(var.xy = xy, equal_var.x = TRUE, equal_var.y = TRUE,
              c1 = 1.2, c2 = 1.2, ncc = 2, max_iter = 30, L2_sol_init = TRUE)
  expect_length(fit$result, 2)
  for (component in fit$result) {
    expect_lte(sum(abs(component$a)), 1.200001)
    expect_lte(sum(abs(component$b)), 1.200001)
    expect_equal(sum(component$a^2), 1)
    expect_true(all(is.finite(component$obj_val_record)))
  }
})

test_that("single-feature iterative fits and zero covariance are supported", {
  fit <- dCCA(var.xy = matrix(.5), max_iter = 5)
  expect_equal(fit$result$rho, .5)
  zero <- dCCA(var.xy = matrix(0, 2, 2), equal_var.x = TRUE,
               equal_var.y = TRUE, c1 = 1.2, c2 = 1.2, max_iter = 5)
  expect_equal(zero$result$obj_val, 0)
  expect_true(all(is.finite(zero$result$a)))
})
