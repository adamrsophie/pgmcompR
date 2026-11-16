test_that("MLE and Bayesian estimates match hand computations", {
  d <- data.frame(A = c("x", "x", "y", "y", "y"), B = c("u", "v", "v", "v", "u"))
  bn <- bayesian_network(list(c("A", "B")))
  f <- fit_parameters(bn, d)
  expect_equal(factor_values(get_cpd(f, "A")), c(2, 3) / 5)
  expect_equal(factor_values(get_cpd(f, "B")), c(1 / 2, 1 / 2, 1 / 3, 2 / 3))
  fb <- fit_parameters(bn, d, method = "bayes", prior = "bdeu", equivalent_sample_size = 4)
  # BDeu: alpha = 4 / (2 * 2) = 1 per cell for B, 4 / 2 = 2 for A
  expect_equal(factor_values(get_cpd(fb, "A")), c(2 + 2, 3 + 2) / 9)
  expect_equal(factor_values(get_cpd(fb, "B")), c(2 / 4, 2 / 4, 2 / 5, 3 / 5))
  fk <- fit_parameters(bn, d, method = "bayes", pseudo_counts = list(B = matrix(c(1, 0, 0, 1), 2)))
  expect_equal(factor_values(get_cpd(fk, "B")), c(2 / 3, 1 / 3, 1 / 4, 3 / 4))
})

test_that("unobserved parent configurations follow an explicit policy", {
  d <- data.frame(A = factor(c("x", "x"), levels = c("x", "y")), B = c("u", "v"))
  bn <- bayesian_network(list(c("A", "B")))
  expect_error(fit_parameters(bn, d, unobserved = "error"), class = "pgm_unobserved_configuration")
  expect_warning(f <- fit_parameters(bn, d), class = "pgm_unobserved_configuration")
  expect_equal(factor_values(get_cpd(f, "B"))[3:4], c(0.5, 0.5))
  fn <- suppressWarnings(fit_parameters(bn, d, unobserved = "na"))
  expect_true(anyNA(factor_values(get_cpd(fn, "B"))))
  expect_error(check_model(fn))
})

test_that("structural zeros are respected", {
  d <- data.frame(A = c("x", "y", "y"), B = c("u", "u", "v"))
  bn <- bayesian_network(list(c("A", "B")))
  z <- matrix(c(FALSE, TRUE, FALSE, FALSE), 2)
  f <- fit_parameters(bn, d, method = "bayes", prior = "k2", structural_zeros = list(B = z))
  m <- matrix(factor_values(get_cpd(f, "B")), 2)
  expect_identical(m[2, 1], 0)
  expect_equal(colSums(m), c(1, 1))
})

test_that("missing data require EM, which recovers parameters", {
  bn <- example_network("asia")
  d <- simulate_bn(bn, 4000, seed = 7)
  expect_error(fit_parameters(bayesian_network(pgmcompR::edges(bn)), replace(d, cbind(1, 1), NA)), class = "pgm_missing_data")
  set.seed(8)
  dm <- d
  dm$smoke[sample(4000, 1000)] <- NA
  em <- fit_parameters(bayesian_network(pgmcompR::edges(bn)), dm, method = "em", max_iter = 50, seed = 1)
  ll <- attr(em, "loglik")
  expect_true(all(diff(ll) > -1e-6))
  expect_lt(max(abs(factor_values(get_cpd(em, "smoke")) - 0.5)), 0.03)
  expect_lt(abs(matrix(factor_values(get_cpd(em, "lung")), 2)[2, 2] - 0.1), 0.03)
})

test_that("linear Gaussian MLE recovers coefficients", {
  g <- example_network("gaussian")
  d <- simulate_bn(g, 5000, seed = 1)
  f <- fit_parameters(bayesian_network(pgmcompR::edges(g)), d)
  x3 <- get_cpd(f, "X3")
  expect_equal(unname(x3$coefficients[c("X1", "X2")]), c(-0.4, 1.5), tolerance = 0.05)
  expect_equal(x3$sd, 0.5, tolerance = 0.05)
  ls <- lm(X3 ~ X1 + X2, data = d)
  expect_equal(unname(x3$coefficients[c("X1", "X2")]), unname(coef(ls)[c("X1", "X2")]))
})

test_that("CI tests behave as expected", {
  d <- simulate_bn(example_network("asia"), 3000, seed = 3)
  expect_lt(ci_test(d, "smoke", "lung")$p_value, 1e-6)
  expect_gt(ci_test(d, "asia", "smoke")$p_value, 0.001)
  expect_gt(ci_test(d, "bronc", "lung", "smoke")$p_value, 0.001)
  set.seed(1)
  x <- rnorm(500); y <- x + rnorm(500); z <- y + rnorm(500)
  dd <- data.frame(x, y, z)
  expect_lt(ci_test(dd, "x", "z", test = "cor")$p_value, 1e-6)
  expect_gt(ci_test(dd, "x", "z", "y", test = "cor")$p_value, 0.001)
})

test_that("structure learners recover a simple network", {
  sp <- example_network("sprinkler")
  d <- simulate_bn(sp, 5000, seed = 1)
  # Greedy hill climbing can stop in a local optimum (bnlearn's hc does too
  # on this sample); tabu search escapes it and MMHC avoids it.
  hc <- hill_climb(d)
  expect_equal(attr(hc, "score"), network_score(hc, d, "bic"))
  for (g in list(tabu_search(d), mmhc(d))) {
    expect_true(markov_equivalent(g, sp))
    expect_gte(attr(g, "score"), attr(hc, "score"))
  }
  p <- pc_algorithm(d)
  expect_true(identical(sort(paste(pmin(p$edges$from, p$edges$to), pmax(p$edges$from, p$edges$to))),
                        sort(paste(pmin(sp$graph$edges$from, sp$graph$edges$to), pmax(sp$graph$edges$from, sp$graph$edges$to)))))
  bl <- hill_climb(d, blacklist = list(c("Cloudy", "Rain"), c("Rain", "Cloudy")))
  expect_false(any(bl$edges$from %in% c("Cloudy", "Rain") & bl$edges$to %in% c("Cloudy", "Rain")))
  wl <- hill_climb(d, whitelist = list(c("Cloudy", "WetGrass")))
  expect_true(any(wl$edges$from == "Cloudy" & wl$edges$to == "WetGrass"))
  tl <- hill_climb(d, tiers = list("WetGrass", c("Cloudy", "Sprinkler", "Rain")))
  expect_false(any(tl$edges$to == "WetGrass"))
})

test_that("scores are decomposable and BDeu is score equivalent", {
  d <- simulate_bn(example_network("sprinkler"), 1000, seed = 4)
  g1 <- dag(list(c("Cloudy", "Rain")), nodes = c("Sprinkler", "WetGrass"))
  g2 <- dag(list(c("Rain", "Cloudy")), nodes = c("Sprinkler", "WetGrass"))
  expect_equal(network_score(g1, d, "bdeu"), network_score(g2, d, "bdeu"))
  expect_equal(network_score(g1, d, "bic"), network_score(g2, d, "bic"))
  expect_gt(network_score(example_network("sprinkler"), d, "bic"), network_score(g1, d, "bic"))
})

test_that("bootstrap edge strengths", {
  d <- simulate_bn(example_network("sprinkler"), 800, seed = 5)
  s <- bootstrap_structure(d, R = 6, seed = 2)
  expect_true(all(s$strength >= 0 & s$strength <= 1))
  expect_s3_class(averaged_network(s, 0.5), "pgm_dag")
})
