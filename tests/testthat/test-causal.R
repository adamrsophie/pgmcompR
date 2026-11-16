confounded <- function() {
  st <- list(Z = c("0", "1"), X = c("0", "1"), Y = c("0", "1"))
  bn <- bayesian_network(list(c("Z", "X"), c("Z", "Y"), c("X", "Y")), states = st)
  bn <- add_cpd(bn, cpd("Z", c(`0` = 0.6, `1` = 0.4)))
  bn <- add_cpd(bn, cpd("X", c(0.8, 0.2, 0.3, 0.7), "Z", st))
  bn <- add_cpd(bn, cpd("Y", c(0.9, 0.1, 0.6, 0.4, 0.5, 0.5, 0.2, 0.8), c("X", "Z"), st))
  bn
}

test_that("do queries use truncated factorization and differ from conditioning", {
  bn <- confounded()
  # P(Y = 1 | do(X = 1)) = sum_z P(Y = 1 | X = 1, z) P(z)
  truth <- 0.4 * 0.6 + 0.8 * 0.4
  expect_equal(factor_values(do_query(bn, "Y", list(X = "1")))[2], truth)
  for (m in c("jt", "enumeration")) expect_equal(factor_values(do_query(bn, "Y", list(X = "1"), method = m))[2], truth)
  obs <- factor_values(query(bn, "Y", list(X = "1")))[2]
  expect_false(isTRUE(all.equal(obs, truth)))
  # Back door adjustment from the observational distribution gives the same value
  adj <- sum(vapply(c("0", "1"), function(z) {
    factor_values(query(bn, "Y", list(X = "1", Z = z)))[2] * factor_values(query(bn, "Z"))[match(z, c("0", "1"))]
  }, 1))
  expect_equal(adj, truth)
  mut <- intervene(bn, list(X = "1"))
  expect_identical(parents(mut, "X"), character())
  s <- simulate_bn(bn, 2000, interventions = list(X = "1"), seed = 1)
  expect_equal(mean(s$Y == "1"), truth, tolerance = 0.04)
})

test_that("adjustment criteria", {
  g <- dag(list(c("Z", "X"), c("Z", "Y"), c("X", "Y")))
  expect_true(is_backdoor_set(g, "X", "Y", "Z"))
  expect_false(is_backdoor_set(g, "X", "Y", character()))
  expect_identical(adjustment_sets(g, "X", "Y"), list("Z"))
  # M-bias: adjusting for the collider M opens a path
  m <- dag(list(c("A", "X"), c("A", "M"), c("B", "M"), c("B", "Y"), c("X", "Y")))
  expect_true(is_backdoor_set(m, "X", "Y", character()))
  expect_false(is_backdoor_set(m, "X", "Y", "M"))
  expect_true(is_backdoor_set(m, "X", "Y", c("M", "A")))
  fd <- dag(list(c("U", "X"), c("U", "Y"), c("X", "M"), c("M", "Y")))
  expect_true(is_frontdoor_set(fd, "X", "Y", "M"))
  id <- identify_effect(fd, "X", "Y", latent = "U")
  expect_identical(id$strategy, "frontdoor")
  bow <- dag(list(c("U", "X"), c("U", "Y"), c("X", "Y")))
  expect_error(identify_effect(bow, "X", "Y", latent = "U"), class = "pgm_not_identifiable")
  # descendants of the treatment are never adjustment variables
  d <- dag(list(c("X", "M"), c("M", "Y")))
  expect_false(is_backdoor_set(d, "X", "Y", "M"))
})

test_that("front door formula reproduces the interventional distribution", {
  st <- list(U = c("0", "1"), X = c("0", "1"), M = c("0", "1"), Y = c("0", "1"))
  bn <- bayesian_network(list(c("U", "X"), c("U", "Y"), c("X", "M"), c("M", "Y")), states = st)
  bn <- add_cpd(bn, cpd("U", c(`0` = 0.5, `1` = 0.5)))
  bn <- add_cpd(bn, cpd("X", c(0.8, 0.2, 0.2, 0.8), "U", st))
  bn <- add_cpd(bn, cpd("M", c(0.9, 0.1, 0.25, 0.75), "X", st))
  bn <- add_cpd(bn, cpd("Y", c(0.9, 0.1, 0.4, 0.6, 0.6, 0.4, 0.15, 0.85), c("M", "U"), st))
  truth <- factor_values(do_query(bn, "Y", list(X = "1")))[2]
  pm <- function(m, x) factor_values(query(bn, "M", list(X = x)))[match(m, c("0", "1"))]
  py <- function(x, m) factor_values(query(bn, "Y", list(X = x, M = m)))[2]
  px <- factor_values(query(bn, "X"))
  fd <- sum(vapply(c("0", "1"), function(m) pm(m, "1") * sum(px * c(py("0", m), py("1", m))), 1))
  expect_equal(fd, truth)
})

test_that("SCM counterfactuals follow abduction, action, prediction", {
  m <- scm(exogenous = list(U = c(a = 0.6, b = 0.4), V = c(a = 0.7, b = 0.3)),
           parents = list(X = "U", Y = c("X", "V")), states = list(X = c("0", "1"), Y = c("0", "1")),
           equations = list(X = function(U) if (U == "b") "1" else "0",
                            Y = function(X, V) if (X == "1" || V == "b") "1" else "0"))
  # Given X = 1, Y = 1: U = b; V keeps its prior. Under do(X = 0), Y = 1 iff V = b.
  cf <- counterfactual(m, "Y", list(X = "1", Y = "1"), list(X = "0"))
  expect_equal(factor_values(cf), c(0.7, 0.3))
  expect_error(counterfactual(confounded(), "Y"), "structural causal model")
})
