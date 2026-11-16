roundtrip_equal <- function(a, b) {
  for (v in nodes(a)) {
    expect_equal(factor_values(get_cpd(b, v) |> factor_reorder(c(v, cpd_parents(get_cpd(a, v))))), factor_values(get_cpd(a, v)), tolerance = 1e-15)
    expect_identical(b$registry[[v]]$states, a$registry[[v]]$states)
  }
}

test_that("BIF, XMLBIF, UAI, JSON and RDS round trips preserve the model", {
  for (bn in list(example_network("asia"), example_network("student"), random_bn(5, seed = 3))) {
    f <- tempfile()
    write_bif(bn, f); roundtrip_equal(bn, read_bif(f))
    skip_if_not_installed("xml2")
    write_xmlbif(bn, f); roundtrip_equal(bn, read_xmlbif(f))
    write_uai(bn, f); roundtrip_equal(bn, read_uai(f, names = nodes(bn), states = model_states(bn)[nodes(bn)]))
    skip_if_not_installed("jsonlite")
    write_json_model(bn, f); roundtrip_equal(bn, read_json_model(f))
    save_model(bn, f); expect_identical(load_model(f), bn)
  }
})

test_that("JSON stores Gaussian, hybrid and Markov models", {
  skip_if_not_installed("jsonlite")
  f <- tempfile()
  g <- example_network("gaussian")
  write_json_model(g, f)
  expect_equal(gaussian_moments(gaussian_joint(read_json_model(f))), gaussian_moments(gaussian_joint(g)))
  st <- list(A = c("0", "1"), B = c("x", "y", "z"))
  mn <- markov_network(list(discrete_factor(c("A", "B"), st, 1:6)))
  write_json_model(mn, f)
  expect_equal(partition_function(read_json_model(f)), 21)
})

test_that("unsupported features are reported, not dropped", {
  f <- tempfile()
  expect_error(write_bif(example_network("gaussian"), f), class = "pgm_unsupported_feature")
  st <- list(A = c("has space", "b"))
  bn <- bayesian_network(nodes = "A", states = st)
  bn <- add_cpd(bn, cpd("A", c(0.5, 0.5), states = st))
  expect_error(write_bif(bn, f), class = "pgm_unsupported_feature")
})

test_that("BIF files from other tools parse", {
  f <- tempfile()
  writeLines(c("network test {", "}", "variable A {", "  type discrete [ 2 ] { t, f };", "}",
               "variable B {", "  type discrete [ 2 ] { t, f };", "}",
               "probability ( A ) {", "  table 0.2, 0.8;", "}",
               "probability ( B | A ) {", "  (t) 0.9, 0.1;", "  (f) 0.3, 0.7;", "}"), f)
  bn <- read_bif(f)
  expect_equal(factor_values(query(bn, "B")), c(0.2 * 0.9 + 0.8 * 0.3, 0.2 * 0.1 + 0.8 * 0.7))
})
