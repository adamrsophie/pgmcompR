st <- list(A = c("a0", "a1"), B = c("b0", "b1", "b2"), C = c("c0", "c1"))

test_that("factor construction checks lengths and values", {
  expect_error(discrete_factor("A", st["A"], c(0.1, 0.2, 0.3)), "length")
  expect_error(discrete_factor("A", st["A"], c(-1, 2)), "non-negative")
  f <- discrete_factor(c("A", "B"), st[c("A", "B")], 1:6)
  expect_identical(as.array(f)["a1", "b2"], 6)
  arr <- array(1:6, c(2, 3), dimnames = list(A = st$A, B = st$B))
  expect_equal(discrete_factor(values = arr), f)
})

test_that("product matches a direct elementwise computation", {
  f <- discrete_factor(c("A", "B"), st[c("A", "B")], 1:6)
  g <- discrete_factor(c("C", "B"), st[c("C", "B")], c(1, 10, 2, 20, 3, 30))
  p <- factor_product(f, g)
  expect_identical(factor_scope(p), c("A", "B", "C"))
  arr <- as.array(p)
  for (a in st$A) for (b in st$B) for (c in st$C) {
    expect_equal(arr[a, b, c], as.array(f)[a, b] * as.array(g)[c, b])
  }
  expect_equal(f * g, p)
})

test_that("mismatched states are rejected, never recycled", {
  f <- discrete_factor("B", st["B"], 1:3)
  g <- discrete_factor("B", list(B = c("b0", "b1")), 1:2)
  expect_error(factor_product(f, g), class = "pgm_state_mismatch")
  h <- discrete_factor("B", list(B = c("b2", "b1", "b0")), 3:1)
  expect_error(factor_product(f, h), class = "pgm_state_mismatch")
  expect_equal(factor_reorder_states(h, "B", st$B), f)
})

test_that("marginalize, maximize, reduce, reorder", {
  f <- discrete_factor(c("A", "B"), st[c("A", "B")], 1:6)
  expect_equal(factor_values(factor_marginalize(f, "B")), c(1 + 3 + 5, 2 + 4 + 6))
  expect_equal(factor_values(factor_marginalize(f, "A")), c(3, 7, 11))
  mx <- factor_maximize(f, "B", argmax = TRUE)
  expect_equal(factor_values(mx), c(5, 6))
  expect_identical(attr(mx, "argmax")$B, c("b2", "b2"))
  r <- factor_reduce(f, list(B = "b1"))
  expect_equal(factor_values(r), c(3, 4))
  expect_error(factor_reduce(f, list(B = "zz")), class = "pgm_unknown_state")
  ro <- factor_reorder(f, c("B", "A"))
  expect_equal(as.array(ro), t(as.array(f)))
  expect_equal(factor_values(factor_marginalize(f, c("A", "B"))), 21)
})

test_that("normalize, divide and zero handling", {
  f <- discrete_factor("A", st["A"], c(0, 0))
  expect_error(factor_normalize(f), class = "pgm_impossible_evidence")
  num <- discrete_factor(c("A", "B"), st[c("A", "B")], c(0, 1, 2, 3, 4, 5))
  den <- discrete_factor("A", st["A"], c(0, 2))
  expect_error(factor_divide(num, den), class = "pgm_division_by_zero")
  num2 <- discrete_factor(c("A", "B"), st[c("A", "B")], c(0, 1, 0, 3, 0, 5))
  q <- factor_divide(num2, den)
  expect_equal(factor_values(q), c(0, 0.5, 0, 1.5, 0, 2.5))
})

test_that("log space operations agree with linear space", {
  f <- discrete_factor(c("A", "B"), st[c("A", "B")], c(1e-200, 2e-200, 3, 4, 5, 6))
  g <- discrete_factor(c("B", "C"), st[c("B", "C")], c(1, 2, 3, 4, 5, 6))
  lin <- factor_marginalize(factor_product(f, g), "B")
  lg <- factor_marginalize(factor_product(factor_to_log(f), factor_to_log(g)), "B")
  expect_true(lg$log)
  expect_equal(exp(factor_values(lg)), factor_values(lin), tolerance = 1e-12)
  tiny <- discrete_factor("A", st["A"], c(-2000, -2001), log = TRUE)
  n <- factor_normalize(tiny)
  expect_equal(exp(factor_values(n)), c(1, exp(-1)) / (1 + exp(-1)))
})

test_that("information measures", {
  f <- discrete_factor(c("A", "C"), st[c("A", "C")], c(0.4, 0.1, 0.1, 0.4))
  expect_equal(factor_entropy(discrete_factor("A", st["A"], c(1, 1)), base = 2), 1)
  hx <- -2 * 0.5 * log(0.5)
  hxy <- -sum(c(0.4, 0.1, 0.1, 0.4) * log(c(0.4, 0.1, 0.1, 0.4)))
  expect_equal(mutual_information(f, "A", "C"), 2 * hx - hxy)
  p <- discrete_factor("A", st["A"], c(0.5, 0.5))
  q <- discrete_factor("A", st["A"], c(0.9, 0.1))
  expect_equal(kl_divergence(p, q), 0.5 * log(0.5 / 0.9) + 0.5 * log(0.5 / 0.1))
  expect_identical(kl_divergence(p, discrete_factor("A", st["A"], c(1, 0))), Inf)
})

test_that("CPDs validate normalization", {
  expect_error(cpd("A", c(a0 = 0.5, a1 = 0.6)), class = "pgm_invalid_cpd")
  expect_error(cpd("A", matrix(c(0, 0, 0.5, 0.5), 2), parents = "C", states = c(st["C"], list(A = st$A))), class = "pgm_invalid_cpd")
  x <- cpd("A", c(a0 = 2, a1 = 6), normalize = TRUE)
  expect_equal(factor_values(x), c(0.25, 0.75))
  no <- noisy_or_cpd("A", c("B2", "C"), list(A = c("off", "on"), B2 = c("n", "y"), C = c("n", "y")),
                     probs = c(B2 = 0.8, C = 0.5), leak = 0.1)
  m <- matrix(factor_values(no), 2)
  expect_equal(m[1, 4], 0.9 * 0.2 * 0.5)
  expect_equal(m[1, 1], 0.9)
})
