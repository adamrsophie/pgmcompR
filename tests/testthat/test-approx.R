test_that("sampling methods approximate exact posteriors", {
  sp <- example_network("sprinkler")
  ex <- factor_values(query(sp, "Rain", list(WetGrass = "T")))
  for (m in c("rejection", "likelihood_weighting", "importance")) {
    q <- query(sp, "Rain", list(WetGrass = "T"), method = m, n = 20000, seed = 1)
    expect_s3_class(q, "pgm_approximate")
    d <- attr(q, "diagnostics")
    expect_true(all(abs(factor_values(q) - ex) < 5 * d$mcse + 1e-3))
  }
  g <- suppressWarnings(query(sp, "Rain", list(WetGrass = "T"), method = "gibbs", n = 4000, seed = 2))
  expect_true(max(abs(factor_values(g) - ex)) < 0.03)
  d <- attr(g, "diagnostics")
  expect_true(all(c("rhat", "ess", "mcse") %in% names(d)))
})

test_that("importance sampling with a different proposal is consistent", {
  bn <- example_network("student")
  prop <- bn
  prop <- add_cpd(prop, cpd("Intelligence", c(low = 0.4, high = 0.6)))
  ex <- factor_values(query(bn, "Intelligence", list(Letter = "weak")))
  q <- query(bn, "Intelligence", list(Letter = "weak"), method = "importance", proposal = prop, n = 20000, seed = 3)
  expect_lt(max(abs(factor_values(q) - ex)), 0.02)
})

test_that("seeds make sampling reproducible and do not leak", {
  sp <- example_network("sprinkler")
  set.seed(99); a <- runif(1)
  q1 <- query(sp, "Rain", list(WetGrass = "T"), method = "likelihood_weighting", n = 500, seed = 5)
  set.seed(99); invisible(query(sp, "Rain", list(WetGrass = "T"), method = "likelihood_weighting", n = 500, seed = 5)); b <- runif(1)
  expect_identical(a, b)
  q2 <- query(sp, "Rain", list(WetGrass = "T"), method = "likelihood_weighting", n = 500, seed = 5)
  expect_identical(factor_values(q1), factor_values(q2))
})

test_that("loopy BP is exact on trees and approximate on loops", {
  ca <- example_network("cancer")
  bp <- belief_propagation(ca, list(Xray = "positive"))
  expect_true(bp$exact)
  for (v in c("Smoker", "Pollution", "Dyspnoea")) {
    expect_equal(factor_values(bp$beliefs[[v]]), factor_values(query(ca, v, list(Xray = "positive"))), tolerance = 1e-8)
  }
  sp <- example_network("sprinkler")
  expect_error(belief_propagation(sp, list(WetGrass = "T")), class = "pgm_loopy_graph")
  lb <- loopy_bp(sp, list(WetGrass = "T"))
  expect_false(lb$exact)
  expect_lt(max(abs(factor_values(lb$beliefs$Rain) - factor_values(query(sp, "Rain", list(WetGrass = "T"))))), 0.1)
  mx <- belief_propagation(ca, list(Xray = "positive"), type = "max")
  expect_identical(unname(mx$assignment[names(mpe_query(ca, list(Xray = "positive"))$assignment)]),
                   unname(mpe_query(ca, list(Xray = "positive"))$assignment))
})

test_that("non convergence of loopy BP warns", {
  st <- list(A = c("0", "1"), B = c("0", "1"), C = c("0", "1"))
  pot <- function(a, b) discrete_factor(c(a, b), st[c(a, b)], c(50, 1, 1, 50))
  mn <- markov_network(list(pot("A", "B"), pot("B", "C"), pot("C", "A"), discrete_factor("A", st["A"], c(2, 1))))
  expect_warning(belief_propagation(mn, max_iter = 2, allow_loops = TRUE), class = "pgm_convergence_warning")
})

test_that("MCMC diagnostics behave sensibly", {
  set.seed(1)
  good <- matrix(rnorm(4000), 1000, 4)
  d <- pgmcompR:::mcmc_diagnostics(good)
  expect_lt(d$rhat, 1.01)
  expect_gt(d$ess, 2000)
  bad <- cbind(rnorm(1000), rnorm(1000, 5))
  expect_gt(pgmcompR:::mcmc_diagnostics(bad)$rhat, 1.5)
})
