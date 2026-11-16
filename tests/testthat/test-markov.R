st <- list(A = c("0", "1"), B = c("0", "1"), C = c("0", "1"))
pot <- function(a, b, v = c(3, 1, 1, 3)) discrete_factor(c(a, b), st[c(a, b)], v)

test_that("partition function and queries in a Markov network", {
  mn <- markov_network(list(pot("A", "B"), pot("B", "C"), pot("C", "A")))
  grid <- expand.grid(st, stringsAsFactors = FALSE)
  w <- apply(grid, 1, function(r) prod(ifelse(r["A"] == r["B"], 3, 1), ifelse(r["B"] == r["C"], 3, 1), ifelse(r["C"] == r["A"], 3, 1)))
  expect_equal(partition_function(mn), sum(w))
  sel <- grid$B == "1"
  ref <- c(sum(w[sel & grid$A == "0"]), sum(w[sel & grid$A == "1"])) / sum(w[sel])
  for (m in c("ve", "jt", "enumeration")) expect_equal(factor_values(query(mn, "A", list(B = "1"), method = m)), ref)
  expect_equal(evidence_probability(mn, list(B = "1")), sum(w[sel]) / sum(w))
  expect_true(check_model(mn))
  expect_s3_class(add_factor(mn, discrete_factor("A", st["A"], c(1, 2))), "markov_network")
  g <- query(mn, "A", list(B = "1"), method = "gibbs", n = 3000, seed = 1)
  expect_lt(max(abs(factor_values(g) - ref)), 0.03)
})

test_that("factor graph conversions preserve the distribution", {
  bn <- example_network("asia")
  fg <- as_factor_graph(bn)
  mn <- as_markov_network(bn)
  expect_equal(partition_function(fg), 1)
  for (m in list(fg, mn)) {
    expect_equal(factor_values(query(m, "lung", list(xray = "yes"))), oracle_query(bn, "lung", list(xray = "yes")))
  }
  expect_s3_class(junction_tree(mn), "junction_tree")
})

test_that("BP is exact on a tree Markov network and MAP agrees", {
  mn <- markov_network(list(pot("A", "B", c(5, 1, 2, 4)), pot("B", "C", c(1, 3, 2, 2))))
  bp <- belief_propagation(mn)
  expect_true(bp$exact)
  for (v in c("A", "B", "C")) expect_equal(factor_values(bp$beliefs[[v]]), factor_values(query(mn, v)), tolerance = 1e-10)
  mx <- belief_propagation(mn, type = "max")
  expect_identical(unname(mx$assignment), unname(mpe_query(mn)$assignment[c("A", "B", "C")]))
})
