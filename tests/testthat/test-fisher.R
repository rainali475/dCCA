test_that("Fisher transformation clips endpoints and uses raw correlations", {
  transformed <- dCCA:::safe_fisher_transform(c(-1, -.5, 0, .5, 1))
  expect_true(all(is.finite(transformed)))
  expect_equal(transformed[2:4], atanh(c(-.5, 0, .5)))
  expect_equal(transformed, -rev(transformed))

  fit <- dCCA(var.xy = array(c(.8, .4), c(1, 1, 2)), z = c(1, 1),
              equal_var.x = TRUE, equal_var.y = TRUE, fisher.transform = TRUE)
  expect_equal(fit$result$rho, c(.8, .4))
  expect_equal(fit$result$obj_val, sum(atanh(c(.8, .4))))
  expect_true(fit$parameters$fisher.transform)
  expect_equal(unname(fit$parameters$method), c("CG", "CG"))
})

test_that("Fisher analytic gradients match finite differences", {
  xy <- array(c(diag(c(.6, .2)), diag(c(.1, -.3))), c(2, 2, 2))
  vx <- array(rep(diag(c(2, 1)), 2), c(2, 2, 2))
  vy <- array(rep(diag(2), 2), c(2, 2, 2))
  a <- c(.7, -.2); b <- c(.4, .8); z <- c(-1, 1)
  numeric_gradient <- vapply(seq_along(a), function(i) {
    step <- c(0, 0); step[i] <- 1e-6
    (dCCA:::obj_fn(a + step, b, z, xy, vx, vy, TRUE) -
       dCCA:::obj_fn(a - step, b, z, xy, vx, vy, TRUE)) / 2e-6
  }, numeric(1))
  expect_equal(dCCA:::obj_gr(a, b, z, xy, vx, vy, TRUE), numeric_gradient,
               tolerance = 1e-6)

  # Above the clipping limit the objective is locally constant.
  xy[] <- 0
  xy[, , 1] <- diag(c(3, 0))
  expect_equal(dCCA:::obj_gr(c(1, 0), c(1, 0), c(1, 0), xy, vx, vy, TRUE),
               c(0, 0))
})

test_that("Fisher mode propagates through sparse starts and deflation", {
  for (solver in c("CG", "Nelder-Mead")) {
    fit <- dCCA(var.xy = diag(c(.7, .2)), equal_var.x = TRUE,
                equal_var.y = TRUE, c1 = 1.3, c2 = 1.3, ncc = 2,
                a_inits = list(c(1, .1), c(.1, 1)),
                b_inits = list(c(1, .1), c(.1, 1)), L2_sol_init = TRUE,
                fisher.transform = TRUE, solver = solver, max_iter = 10)
    expect_length(fit$result, 2)
    for (result in fit$result) {
      expect_equal(result$obj_val, sum(dCCA:::safe_fisher_transform(result$rho)))
      expect_true(all(is.finite(result$obj_val_record)))
      expect_lte(sum(abs(result$a)), 1.300001)
      expect_lte(sum(abs(result$b)), 1.300001)
    }
  }
})

test_that("Fisher flag is validated and Fisher mode is the default", {
  for (flag in list(NA, 1, c(TRUE, FALSE), NULL))
    expect_error(dCCA(var.xy = diag(2), fisher.transform = flag), "single logical")
  expect_identical(dCCA(var.xy = diag(c(.8, .3))),
                   dCCA(var.xy = diag(c(.8, .3)), fisher.transform = TRUE))
  fit <- dCCA(var.xy = matrix(.5), equal_var.x = TRUE, equal_var.y = TRUE)
  expect_true(fit$parameters$fisher.transform)
  expect_equal(fit$result$obj_val, atanh(.5))
  expect_equal(fit$result$rho, .5)
})
