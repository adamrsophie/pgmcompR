test_that("exact methods match the brute force oracle on Asia", {
  bn <- example_network("asia")
  cases <- list(
    list(v = "lung", e = list(xray = "yes", smoke = "yes")),
    list(v = c("tub", "bronc"), e = list(dysp = "yes")),
    list(v = "either", e = list()),
    list(v = c("asia", "smoke"), e = list(xray = "yes", dysp = "no"))
  )
  for (cs in cases) {
    ref <- oracle_query(bn, cs$v, cs$e)
    for (m in c("ve", "enumeration", "jt")) {
      q <- query(bn, cs$v, cs$e, method = m)
      expect_equal(factor_values(factor_reorder(q, cs$v)), ref, tolerance = 1e-12, label = m)
    }
    expect_equal(factor_values(query(bn, cs$v, cs$e, log = TRUE)), log(ref), tolerance = 1e-10)
  }
})

test_that("exact methods match the oracle on random networks with zeros", {
  for (s in 1:6) {
    bn <- random_bn(7, max_card = 3, seed = s, zeros = TRUE)
    vars <- nodes(bn)
    q <- vars[c(2, 6)]
    ev <- list()
    ev[[vars[7]]] <- "s1"
    ref <- tryCatch(oracle_query(bn, q, ev), error = function(e) NULL)
    if (is.null(ref) || anyNA(ref)) next
    for (ord in c("min_fill", "min_degree", "weighted_min_fill", "min_weight")) {
      expect_equal(factor_values(query(bn, q, ev, elimination_order = ord)), ref, tolerance = 1e-12)
    }
    expect_equal(factor_values(query(bn, q, ev, method = "jt")), ref, tolerance = 1e-12)
    expect_equal(factor_values(query(bn, q, ev, method = "enumeration")), ref, tolerance = 1e-12)
    user_ord <- rev(setdiff(vars, c(q, names(ev))))
    expect_equal(factor_values(query(bn, q, ev, elimination_order = user_ord)), ref, tolerance = 1e-12)
  }
})

test_that("query distinguishes observation, intervention and soft evidence", {
  sp <- example_network("sprinkler")
  obs <- query(sp, "Rain", list(Sprinkler = "T"))
  do <- do_query(sp, "Rain", list(Sprinkler = "T"))
  prior <- query(sp, "Rain")
  expect_equal(factor_values(do), factor_values(prior))
  expect_false(isTRUE(all.equal(factor_values(obs), factor_values(prior))))
  # Analytic: P(Rain = T) = 0.5 * 0.2 + 0.5 * 0.8
  expect_equal(factor_values(prior)[2], 0.5)
  # Soft evidence with likelihood (0.3, 0.7) on WetGrass equals adding a
  # virtual child observed with those likelihoods.
  soft <- query(sp, "Rain", virtual_evidence = list(WetGrass = c(F = 0.3, T = 0.7)))
  j <- oracle_joint(sp)
  w <- ifelse(j$WetGrass == "T", 0.7, 0.3) * j$p
  ref <- c(sum(w[j$Rain == "F"]), sum(w[j$Rain == "T"])) / sum(w)
  expect_equal(factor_values(soft), ref, tolerance = 1e-12)
  expect_equal(factor_values(query(sp, "Rain", virtual_evidence = list(WetGrass = c(F = 0.3, T = 0.7)), method = "jt")), ref)
  # Known textbook value P(Rain = T | WetGrass = T) = 0.7079
  expect_equal(factor_values(query(sp, "Rain", list(WetGrass = "T")))[2], 0.70793, tolerance = 1e-4)
})

test_that("impossible and invalid evidence is reported", {
  bn <- example_network("asia")
  expect_error(query(bn, "lung", list(tub = "yes", either = "no")), class = "pgm_impossible_evidence")
  expect_error(query(bn, "lung", list(tub = "yes", either = "no"), method = "jt"), class = "pgm_impossible_evidence")
  expect_error(query(bn, "lung", list(tub = "maybe")), class = "pgm_unknown_state")
  expect_error(query(bn, "lung", list(nope = "yes")), class = "pgm_unknown_variable")
  expect_error(query(bn, "lung", list(lung = "yes")), "observed")
  expect_identical(evidence_probability(bn, list(tub = "yes", either = "no")), 0)
})

test_that("evidence probability and marginals", {
  bn <- example_network("asia")
  j <- oracle_joint(bn)
  ref <- sum(j$p[j$xray == "yes" & j$dysp == "yes"])
  expect_equal(evidence_probability(bn, list(xray = "yes", dysp = "yes")), ref, tolerance = 1e-12)
  expect_equal(joint_probability(bn, as.list(j[1, nodes(bn)])), j$p[1], tolerance = 1e-12)
  expect_equal(factor_values(marginal_probability(bn, "either")), oracle_query(bn, "either"))
  ml <- query(bn, c("lung", "tub"), list(xray = "yes"), joint = FALSE)
  expect_named(ml, c("lung", "tub"))
  expect_equal(factor_values(ml$tub), oracle_query(bn, "tub", list(xray = "yes")))
})

test_that("MAP and MPE match brute force maximization", {
  for (s in 1:4) {
    bn <- random_bn(6, seed = 20 + s)
    vars <- nodes(bn)
    ev <- setNames(list("s1"), vars[6])
    j <- oracle_joint(bn)
    jj <- j[j[[vars[6]]] == "s1", ]
    pe <- sum(jj$p)
    # MPE over all other variables
    mpe <- mpe_query(bn, ev)
    k <- which.max(jj$p)
    expect_equal(mpe$probability, jj$p[k] / pe, tolerance = 1e-10)
    expect_equal(unname(mpe$assignment[vars[1:5]]), unname(unlist(jj[k, vars[1:5]])))
    # MAP over two variables: maximize the marginal
    mv <- vars[c(1, 3)]
    post <- oracle_query(bn, mv, ev)
    mp <- map_query(bn, mv, ev)
    expect_equal(mp$probability, max(post), tolerance = 1e-10)
  }
  # MAP differs from the restriction of MPE in general
  bn <- example_network("asia")
  expect_s3_class(map_query(bn, "lung", list(xray = "yes")), "pgm_map")
})

test_that("junction tree satisfies running intersection and calibration", {
  bn <- example_network("asia")
  jt <- junction_tree(bn)
  # running intersection: for each variable, cliques containing it form a connected subtree
  for (v in nodes(bn)) {
    has <- which(vapply(jt$cliques, function(cl) v %in% cl, logical(1)))
    sub <- jt$tree[jt$tree[, 1] %in% has & jt$tree[, 2] %in% has, , drop = FALSE]
    expect_identical(nrow(sub), length(has) - 1L)
  }
  cal <- calibrate(jt, list(xray = "yes"))
  for (r in seq_len(nrow(jt$tree))) {
    a <- jt$tree[r, 1]; b <- jt$tree[r, 2]
    s <- intersect(jt$cliques[[a]], jt$cliques[[b]])
    ma <- factor_marginalize(cal$beliefs[[a]], setdiff(jt$cliques[[a]], s))
    mb <- factor_marginalize(cal$beliefs[[b]], setdiff(jt$cliques[[b]], s))
    expect_equal(ma, mb, tolerance = 1e-12)
  }
  expect_equal(factor_values(query(jt, "lung", list(xray = "yes"), method = "jt")),
               oracle_query(bn, "lung", list(xray = "yes")), tolerance = 1e-12)
})

test_that("disconnected networks and large cardinalities", {
  st <- list(A = paste0("a", 1:12), B = c("b1", "b2"), C = c("c1", "c2"))
  bn <- bayesian_network(list(c("A", "B")), nodes = "C", states = st)
  set.seed(4)
  pa <- runif(12); bn <- add_cpd(bn, cpd("A", pa / sum(pa), states = st["A"]))
  m <- matrix(runif(24), 2); m <- sweep(m, 2, colSums(m), "/")
  bn <- add_cpd(bn, cpd("B", as.numeric(m), "A", st))
  bn <- add_cpd(bn, cpd("C", c(c1 = 0.3, c2 = 0.7)))
  for (meth in c("ve", "jt", "enumeration")) {
    expect_equal(factor_values(query(bn, "A", list(B = "b2"), method = meth)), oracle_query(bn, "A", list(B = "b2")), tolerance = 1e-12)
    expect_equal(factor_values(query(bn, "C", list(B = "b2"), method = meth)), c(0.3, 0.7))
  }
})

test_that("numerical underflow is handled in log space", {
  # A long chain with tiny probabilities
  n <- 300
  vars <- paste0("X", 1:n)
  st <- setNames(rep(list(c("0", "1")), n), vars)
  bn <- bayesian_network(lapply(2:n, function(i) c(vars[i - 1], vars[i])), states = st)
  bn <- add_cpd(bn, cpd("X1", c(`0` = 0.5, `1` = 0.5)))
  for (i in 2:n) bn <- add_cpd(bn, cpd(vars[i], c(0.999, 0.001, 0.001, 0.999), vars[i - 1], st))
  ev <- as.list(setNames(rep(c("0", "1"), length.out = n - 1), vars[-1]))
  lp <- evidence_probability(bn, ev, log = TRUE)
  expect_true(is.finite(lp) && lp < -1000)
  q <- query(bn, "X1", ev, log = TRUE)
  expect_true(all(is.finite(factor_values(q))))
})
