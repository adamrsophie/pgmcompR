test_that("Bayesian network construction and validation", {
  bn <- bayesian_network(list(c("Smoking", "Cancer"), c("Cancer", "XRay")))
  bn <- add_cpd(bn, "Smoking", c(no = 0.7, yes = 0.3))
  expect_error(check_model(bn), class = "pgm_invalid_model")
  expect_length(check_model(bn, error = FALSE), 2L)
  bn <- add_cpd(bn, "Cancer", matrix(c(0.99, 0.01, 0.9, 0.1), 2, dimnames = list(c("no", "yes"), NULL)))
  bn <- add_cpd(bn, "XRay", matrix(c(0.95, 0.05, 0.2, 0.8), 2, dimnames = list(c("negative", "positive"), NULL)))
  expect_true(check_model(bn))
  expect_identical(variable_registry(bn)$cardinality, c(2L, 2L, 2L))
  expect_error(add_cpd(bn, cpd("Cancer", c(no = 0.5, yes = 0.5))), class = "pgm_invalid_cpd")
  # Parent states in a CPD must match the registry
  bad <- cpd("XRay", c(0.9, 0.1, 0.2, 0.8), "Cancer", list(XRay = c("negative", "positive"), Cancer = c("yes", "no")))
  expect_error(add_cpd(bn, bad), class = "pgm_state_mismatch")
  expect_s3_class(get_cpd(bn, "XRay"), "pgm_cpd")
  expect_error(check_model(remove_cpd(bn, "XRay")), class = "pgm_invalid_model")
  q <- query(bn, "Cancer", list(XRay = "positive"))
  p_yes <- 0.3 * 0.1 * 0.8 + 0.7 * 0.01 * 0.8
  p_no <- 0.3 * 0.9 * 0.05 + 0.7 * 0.99 * 0.05
  expect_equal(factor_values(q), c(p_no, p_yes) / (p_no + p_yes))
})

test_that("arrays with named dimnames define CPDs in any parent order", {
  st <- list(C = c("c0", "c1"), A = c("a0", "a1"), B = c("b0", "b1", "b2"))
  arr <- array(0, c(2, 2, 3), dimnames = st)
  arr[1, , ] <- 0.25; arr[2, , ] <- 0.75
  arr[, "a1", "b2"] <- c(0.9, 0.1)
  bn <- bayesian_network(list(c("B", "C"), c("A", "C")), states = st[c("A", "B")])
  bn <- add_cpd(bn, "A", c(a0 = 0.5, a1 = 0.5))
  bn <- add_cpd(bn, "B", c(b0 = 0.2, b1 = 0.3, b2 = 0.5))
  bn <- add_cpd(bn, "C", arr)
  q <- query(bn, "C", list(A = "a1", B = "b2"))
  expect_equal(factor_values(q), c(0.9, 0.1))
})

test_that("custom and parametric CPDs work in inference", {
  st <- list(Fever = c("no", "yes"), Flu = c("no", "yes"), Cold = c("no", "yes"))
  bn <- bayesian_network(list(c("Flu", "Fever"), c("Cold", "Fever")))
  bn <- add_cpd(bn, cpd("Flu", c(no = 0.9, yes = 0.1)))
  bn <- add_cpd(bn, cpd("Cold", c(no = 0.8, yes = 0.2)))
  bn <- add_cpd(bn, noisy_or_cpd("Fever", c("Flu", "Cold"), st, c(Flu = 0.8, Cold = 0.4), leak = 0.05))
  q <- query(bn, "Fever")
  pf <- 0.95 * (1 - 0.1 * 0.8) * (1 - 0.2 * 0.4)
  expect_equal(factor_values(q), c(pf, 1 - pf))
  d <- deterministic_cpd("Fever", c("Flu", "Cold"), st, function(p) if (any(p == "yes")) "yes" else "no")
  expect_equal(matrix(factor_values(d), 2)[2, ], c(0, 1, 1, 1))
})

test_that("simulation reproduces the joint distribution", {
  bn <- example_network("asia")
  d <- simulate_bn(bn, 40000, seed = 1)
  expect_identical(levels(d$lung), c("no", "yes"))
  expect_equal(mean(d$smoke == "yes"), 0.5, tolerance = 0.02)
  p <- oracle_query(bn, "either")
  expect_equal(mean(d$either == "yes"), p[2], tolerance = 0.1)
  di <- simulate_bn(bn, 2000, interventions = list(either = "yes"), seed = 2)
  expect_true(all(di$either == "yes"))
  expect_equal(mean(di$xray == "yes"), 0.98, tolerance = 0.02)
})

test_that("predict_bn returns posterior modes and probabilities", {
  bn <- example_network("asia")
  nd <- data.frame(smoke = c("yes", "no", NA), xray = c("yes", "no", "yes"))
  pr <- predict_bn(bn, nd, "lung", type = "prob")
  expect_equal(pr[1, ], setNames(factor_values(query(bn, "lung", list(smoke = "yes", xray = "yes"))), c("no", "yes")))
  expect_equal(pr[3, ], setNames(factor_values(query(bn, "lung", list(xray = "yes"))), c("no", "yes")))
  cl <- predict_bn(bn, nd, "lung")
  expect_s3_class(cl, "factor")
})
