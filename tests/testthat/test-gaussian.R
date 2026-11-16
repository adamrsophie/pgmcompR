test_that("Gaussian factor operations match closed forms", {
  S <- matrix(c(2, 0.6, 0.3, 0.6, 1, 0.2, 0.3, 0.2, 1.5), 3)
  mu <- c(1, -1, 0.5)
  f <- gaussian_factor(c("X", "Y", "Z"), mu, S)
  m <- gaussian_moments(gaussian_marginalize(f, "Z"))
  expect_equal(unname(m$covariance), S[1:2, 1:2])
  r <- gaussian_moments(gaussian_reduce(f, c(Z = 2)))
  expect_equal(unname(r$mean), mu[1:2] + S[1:2, 3] / S[3, 3] * (2 - mu[3]))
  expect_equal(unname(r$covariance), S[1:2, 1:2] - S[1:2, 3, drop = FALSE] %*% t(S[1:2, 3, drop = FALSE]) / S[3, 3])
  # Canonical path: marginalize without moments
  fc <- gaussian_factor(f$vars, precision = f$K, h = f$h)
  expect_equal(unname(gaussian_moments(gaussian_marginalize(fc, "Z"))$covariance), S[1:2, 1:2])
  expect_equal(unname(gaussian_moments(gaussian_reduce(fc, c(Z = 2)))$mean), unname(r$mean))
  expect_error(gaussian_factor(c("X", "Y"), c(0, 0), matrix(c(1, 1, 1, 1), 2)), class = "pgm_singular_matrix")
})

test_that("linear Gaussian networks: joint and conditioning", {
  g <- example_network("gaussian")
  j <- gaussian_moments(gaussian_joint(g))
  # X1 ~ N(1, 1); X2 = -0.5 + 0.8 X1 + e2; X3 = 0.2 + 1.5 X2 - 0.4 X1 + e3
  expect_equal(unname(j$mean), c(1, 0.3, 0.25))
  v2 <- 0.64 + 0.36
  c12 <- 0.8
  c13 <- 1.5 * 0.8 - 0.4
  v3 <- 1.5^2 * v2 + 0.16 - 2 * 1.5 * 0.4 * c12 + 0.25
  expect_equal(j$covariance["X2", "X2"], v2)
  expect_equal(j$covariance["X1", "X3"], c13)
  expect_equal(j$covariance["X3", "X3"], v3)
  q <- gaussian_moments(query(g, "X1", list(X3 = 2)))
  expect_equal(unname(q$mean), 1 + c13 / v3 * (2 - 0.25))
  # Canonical CPD factors multiply to the joint density
  prodf <- Reduce(gaussian_product, lapply(g$cpds, lg_cpd_factor))
  expect_equal(unname(gaussian_moments(prodf)$covariance[c("X1", "X2", "X3"), c("X1", "X2", "X3")]), unname(j$covariance))
  # Intervention cuts X1 -> X2
  dq <- gaussian_moments(gaussian_query(g, "X2", interventions = list(X1 = 3)))
  expect_equal(unname(dq$mean), -0.5 + 0.8 * 3)
  dq3 <- gaussian_moments(gaussian_query(g, "X1", interventions = list(X2 = 3)))
  expect_equal(unname(dq3$mean), 1)
})

test_that("conditional Gaussian inference is exact", {
  hb <- bayesian_network(list(c("Sex", "Weight")))
  hb <- add_cpd(hb, cpd("Sex", c(F = 0.4, M = 0.6)))
  hb <- add_cpd(hb, cg_cpd("Weight", "Sex", list(Sex = c("F", "M")), intercept = c(60, 80), sd = c(5, 10)))
  q <- query(hb, "Sex", list(Weight = 70))
  lf <- 0.4 * dnorm(70, 60, 5); lm <- 0.6 * dnorm(70, 80, 10)
  expect_equal(factor_values(q$discrete), c(lf, lm) / (lf + lm))
  w <- query(hb, "Weight")
  expect_equal(unname(w$mean), 0.4 * 60 + 0.6 * 80)
  expect_equal(unname(w$covariance[1, 1]), 0.4 * 25 + 0.6 * 100 + 0.4 * 0.6 * 20^2)
  ws <- query(hb, "Weight", list(Sex = "M"))
  expect_equal(unname(ws$mean), 80)
})

test_that("Gaussian BP gives exact means", {
  J <- matrix(c(3, -1, 0, -1, 3, -1, 0, -1, 3), 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  gm <- gaussian_markov_network(J, c(1, 2, 0))
  bp <- gaussian_bp(gm)
  expect_true(bp$converged && bp$exact)
  expect_equal(unname(bp$mean), as.numeric(solve(J, c(1, 2, 0))), tolerance = 1e-10)
  expect_equal(unname(bp$variance), unname(diag(solve(J))), tolerance = 1e-10)
})
