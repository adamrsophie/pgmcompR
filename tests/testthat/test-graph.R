test_that("DAG construction rejects cycles and duplicates", {
  expect_error(dag(list(c("A", "B"), c("B", "A"))), class = "pgm_error")
  expect_error(dag(list(c("A", "B"), c("B", "C"), c("C", "A"))), class = "pgm_cycle_error")
  expect_error(dag(list(c("A", "A"))), "Self loops")
  expect_error(dag(list(c("A", "B"), c("A", "B"))), "one edge")
  g <- dag(list(c("A", "B")), nodes = c("C", "A"))
  expect_identical(nodes(g), c("C", "A", "B"))
})

test_that("neighborhoods and orderings", {
  g <- dag(list(c("A", "C"), c("B", "C"), c("C", "D"), c("E", "D")))
  expect_identical(parents(g, "C"), c("A", "B"))
  expect_identical(children(g, "C"), "D")
  expect_identical(ancestors(g, "D"), c("A", "B", "C", "E"))
  expect_identical(descendants(g, "A"), c("C", "D"))
  expect_identical(markov_blanket(g, "C"), c("A", "B", "E", "D"))
  ord <- topological_sort(g)
  for (k in seq_len(nrow(g$edges))) expect_lt(match(g$edges$from[k], ord), match(g$edges$to[k], ord))
  expect_true(is_dag(g))
  expect_error(add_edges(g, "D", "A"), class = "pgm_cycle_error")
  g2 <- remove_nodes(g, "C")
  expect_identical(nrow(edges(g2)), 1L)
})

test_that("d-separation matches textbook cases", {
  chain <- dag(list(c("A", "B"), c("B", "C")))
  fork <- dag(list(c("B", "A"), c("B", "C")))
  coll <- dag(list(c("A", "B"), c("C", "B"), c("B", "D")))
  expect_false(d_separated(chain, "A", "C"))
  expect_true(d_separated(chain, "A", "C", "B"))
  expect_true(d_separated(fork, "A", "C", "B"))
  expect_true(d_separated(coll, "A", "C"))
  expect_false(d_separated(coll, "A", "C", "B"))
  expect_false(d_separated(coll, "A", "C", "D"))
  asia <- example_network("asia")
  expect_true(d_separated(asia, "asia", "smoke"))
  expect_false(d_separated(asia, "asia", "smoke", "dysp"))
  expect_true(d_separated(asia, "xray", "dysp", "either"))
  expect_true(d_separated(asia, "tub", "bronc", "smoke"))
})

test_that("d-separation agrees with numerical independence on random networks", {
  for (s in 1:4) {
    bn <- random_bn(5, seed = s)
    vars <- nodes(bn)
    j <- oracle_joint(bn)
    for (x in vars) for (y in vars) {
      if (x >= y) next
      z <- setdiff(vars, c(x, y))[1]
      if (d_separated(bn, x, y, z)) {
        # P(x, y | z) = P(x | z) P(y | z) for every z
        for (zs in bn$registry[[z]]$states) {
          pxy <- oracle_query(bn, c(x, y), setNames(list(zs), z))
          px <- oracle_query(bn, x, setNames(list(zs), z))
          py <- oracle_query(bn, y, setNames(list(zs), z))
          expect_equal(pxy, as.numeric(outer(px, py)), tolerance = 1e-10)
        }
      }
    }
  }
})

test_that("moralization, triangulation and chordality", {
  g <- dag(list(c("A", "C"), c("B", "C")))
  m <- moralize(g)
  expect_s3_class(m, "pgm_ugraph")
  expect_identical(nrow(edges(m)), 3L)
  c4 <- undirected_graph(list(c("A", "B"), c("B", "C"), c("C", "D"), c("D", "A")))
  expect_false(is_chordal(c4))
  tri <- triangulate(c4)
  expect_true(is_chordal(tri$graph))
  expect_identical(length(tri$cliques), 2L)
  expect_true(separated(c4, "A", "C", c("B", "D")))
  expect_false(separated(c4, "A", "C", "B"))
})

test_that("v-structures, CPDAG and Markov equivalence", {
  g1 <- dag(list(c("A", "B"), c("B", "C")))
  g2 <- dag(list(c("C", "B"), c("B", "A")))
  g3 <- dag(list(c("A", "B"), c("C", "B")))
  expect_true(markov_equivalent(g1, g2))
  expect_false(markov_equivalent(g1, g3))
  expect_identical(nrow(v_structures(g3)), 1L)
  cp <- cpdag(g3)
  expect_true(all(edges(cp)$type == "directed"))
  cp1 <- cpdag(g1)
  expect_true(all(edges(cp1)$type == "undirected"))
  # Meek R1 compels B -> D after the v-structure A -> B <- C
  g4 <- dag(list(c("A", "B"), c("C", "B"), c("B", "D")))
  e <- edges(cpdag(g4))
  expect_true(any(e$from == "B" & e$to == "D" & e$type == "directed"))
})

test_that("comparison, isomorphism and adjacency", {
  a <- dag(list(c("A", "B"), c("B", "C")))
  b <- dag(list(c("B", "A"), c("B", "C"), c("A", "C")))
  cmp <- compare_graphs(a, b)
  expect_identical(cmp$shd, 2L)
  expect_true(is_isomorphic(a, dag(list(c("X", "Y"), c("Y", "Z")))))
  expect_false(is_isomorphic(a, dag(list(c("X", "Y"), c("Z", "Y")))))
  m <- adjacency_matrix(a)
  expect_identical(sum(m), 2L)
  expect_true(is_isomorphic(graph_from_adjacency(m), a))
})
