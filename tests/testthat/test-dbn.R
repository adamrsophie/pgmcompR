umbrella <- function() hmm(c("rain", "dry"), c("umbrella", "none"), c(0.5, 0.5),
                           matrix(c(0.7, 0.3, 0.3, 0.7), 2), matrix(c(0.9, 0.1, 0.2, 0.8), 2))

test_that("HMM filtering, smoothing and Viterbi match textbook values", {
  h <- umbrella()
  ev <- data.frame(Observation = c("umbrella", "umbrella"))
  f <- filter_dbn(h, ev)
  # Russell and Norvig: P(R1 | u1) = 0.818, P(R2 | u1, u2) = 0.883
  expect_equal(factor_values(f$marginals[[1]]$Hidden)[1], 0.818, tolerance = 1e-3)
  expect_equal(factor_values(f$marginals[[2]]$Hidden)[1], 0.883, tolerance = 1e-3)
  s <- smooth_dbn(h, ev)
  expect_equal(factor_values(s$marginals[[1]]$Hidden)[1], 0.883, tolerance = 1e-3)
  v <- viterbi_dbn(h, data.frame(Observation = c("umbrella", "umbrella", "none", "umbrella", "umbrella")))
  expect_identical(v$Hidden, c("rain", "rain", "dry", "rain", "rain"))
})

test_that("DBN inference agrees with the unrolled Bayesian network", {
  h <- umbrella()
  obs <- c("umbrella", "none", "umbrella")
  # Unroll three slices by hand
  st <- list(H1 = c("rain", "dry"), H2 = c("rain", "dry"), H3 = c("rain", "dry"),
             O1 = c("umbrella", "none"), O2 = c("umbrella", "none"), O3 = c("umbrella", "none"))
  bn <- bayesian_network(list(c("H1", "H2"), c("H2", "H3"), c("H1", "O1"), c("H2", "O2"), c("H3", "O3")), states = st)
  bn <- add_cpd(bn, cpd("H1", c(rain = 0.5, dry = 0.5)))
  for (k in 2:3) bn <- add_cpd(bn, cpd(paste0("H", k), c(0.7, 0.3, 0.3, 0.7), paste0("H", k - 1), st))
  for (k in 1:3) bn <- add_cpd(bn, cpd(paste0("O", k), c(0.9, 0.1, 0.2, 0.8), paste0("H", k), st))
  ev <- list(O1 = obs[1], O2 = obs[2], O3 = obs[3])
  s <- smooth_dbn(h, data.frame(Observation = obs))
  for (k in 1:3) expect_equal(factor_values(s$marginals[[k]]$Hidden), factor_values(query(bn, paste0("H", k), ev)), tolerance = 1e-12)
  expect_equal(filter_dbn(h, data.frame(Observation = obs))$loglik, evidence_probability(bn, ev, log = TRUE), tolerance = 1e-12)
  fc <- forecast_dbn(h, data.frame(Observation = obs), 1)
  expect_equal(sum(factor_values(fc[[1]]$Hidden)), 1)
  pf <- particle_filter(h, data.frame(Observation = obs), n_particles = 20000, seed = 1)
  expect_lt(max(abs(factor_values(pf$marginals[[3]]$Hidden) - factor_values(filter_dbn(h, data.frame(Observation = obs))$marginals[[3]]$Hidden))), 0.02)
})

test_that("simulate and fit a DBN; time varying transitions", {
  h <- umbrella()
  sim <- simulate_dbn(h, 3000, seed = 1)
  fh <- suppressWarnings(fit_dbn(h, sim))
  expect_equal(matrix(factor_values(fh$transition$Hidden), 2)[1, 1], 0.7, tolerance = 0.04)
  h2 <- add_dbn_cpd(h, "Hidden", matrix(c(0, 1, 0, 1), 2), "transition", time = 2)
  f <- filter_dbn(h2, data.frame(Observation = c("umbrella", "umbrella")))
  expect_equal(factor_values(f$marginals[[2]]$Hidden), c(0, 1))
  expect_error(filter_dbn(h2, data.frame(Observation = c("umbrella", NA))), NA)
})

test_that("Kalman filter matches the conditional Gaussian of the unrolled model", {
  y <- c(1, 2.5, NA, 3)
  kf <- kalman_filter(y, A = 0.9, C = 1, Q = 0.5, R = 0.25, m0 = 0, P0 = 2)
  ks <- kalman_smoother(kf)
  # Build the joint Gaussian of x1..x4, y1..y4 directly
  bn <- bayesian_network(list(c("x1", "x2"), c("x2", "x3"), c("x3", "x4"), c("x1", "y1"), c("x2", "y2"), c("x3", "y3"), c("x4", "y4")))
  bn <- add_cpd(bn, lg_cpd("x1", 0, sd = sqrt(2)))
  for (k in 2:4) bn <- add_cpd(bn, lg_cpd(paste0("x", k), 0, setNames(0.9, paste0("x", k - 1)), sqrt(0.5)))
  for (k in 1:4) bn <- add_cpd(bn, lg_cpd(paste0("y", k), 0, setNames(1, paste0("x", k)), 0.5))
  q <- gaussian_moments(query(bn, c("x1", "x2", "x3", "x4"), list(y1 = 1, y2 = 2.5, y4 = 3)))
  expect_equal(as.numeric(ks$m_smooth), unname(q$mean), tolerance = 1e-10)
  expect_equal(vapply(ks$P_smooth, function(p) p[1, 1], 1), unname(diag(q$covariance)), tolerance = 1e-10)
  qf <- gaussian_moments(query(bn, "x2", list(y1 = 1, y2 = 2.5)))
  expect_equal(kf$m[2, 1], unname(qf$mean), tolerance = 1e-10)
})
