test_that("bnlearn conversions round trip", {
  skip_if_not_installed("bnlearn")
  bn <- example_network("asia")
  fit <- as_bnlearn(bn)
  expect_s3_class(fit, "bn.fit")
  back <- from_bnlearn(fit)
  expect_equal(query(back, "lung", list(xray = "yes")), query(bn, "lung", list(xray = "yes")))
  g <- example_network("gaussian")
  gb <- from_bnlearn(as_bnlearn(g))
  expect_equal(gaussian_moments(gaussian_joint(gb)), gaussian_moments(gaussian_joint(g)))
  s <- from_bnlearn(bnlearn::model2network("[A][B|A][C|B]"))
  expect_true(is_dag(s))
  d <- simulate_bn(example_network("sprinkler"), 2000, seed = 1)
  expect_equal(network_score(learn_structure(d, "hill_climb", engine = "bnlearn"), d, "bic"),
               network_score(learn_structure(d, "hill_climb"), d, "bic"), tolerance = 1e-8)
})

test_that("gRain conversions agree with gRain inference", {
  skip_if_not_installed("gRain")
  bn <- example_network("asia")
  gr <- as_grain(bn)
  ref <- gRain::querygrain(gRain::setEvidence(gr, evidence = list(xray = "yes", smoke = "yes")), "lung")$lung
  expect_equal(factor_values(query(bn, "lung", list(xray = "yes", smoke = "yes"))), as.numeric(ref[c("no", "yes")]), tolerance = 1e-10)
  back <- from_grain(gr)
  expect_equal(query(back, "dysp", list(xray = "yes")), query(bn, "dysp", list(xray = "yes")))
})

test_that("igraph conversions", {
  skip_if_not_installed("igraph")
  bn <- example_network("asia")
  ig <- as_igraph(bn)
  expect_true(igraph::is_dag(ig))
  expect_true(isTRUE(is_isomorphic(as_pgm(ig), bn)))
  u <- as_igraph(moralize(bn))
  expect_false(igraph::is_directed(u))
  expect_s3_class(as_pgm(u), "pgm_ugraph")
})

test_that("plots draw without error", {
  skip_if_not_installed("igraph")
  pdf(NULL)
  on.exit(dev.off())
  expect_s3_class(plot(example_network("asia"), evidence = "xray", interventions = "smoke"), "igraph")
  expect_s3_class(plot(as_factor_graph(example_network("sprinkler"))), "igraph")
  expect_s3_class(plot(junction_tree(example_network("asia"))), "igraph")
  expect_silent(plot(get_cpd(example_network("asia"), "dysp")))
})
