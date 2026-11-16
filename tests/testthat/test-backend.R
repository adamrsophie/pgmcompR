test_that("C++ kernels agree with the R reference implementation", {
  set.seed(2)
  st <- lapply(c(A = 3, B = 2, C = 4, D = 2), function(k) paste0("s", seq_len(k)))
  f <- discrete_factor(c("A", "B", "C"), st[c("A", "B", "C")], runif(24))
  g <- discrete_factor(c("D", "C", "A"), st[c("D", "C", "A")], runif(24))
  for (op in list(
    function() factor_product(f, g),
    function() factor_marginalize(factor_product(f, g), c("C", "A")),
    function() factor_marginalize(factor_to_log(factor_product(f, g)), "B"),
    function() factor_divide(factor_product(f, g), factor_marginalize(g, "D"))
  )) {
    r <- with_backend("R", op())
    cc <- with_backend("cpp", op())
    expect_equal(factor_values(cc), factor_values(r), tolerance = 1e-12)
  }
  bn <- random_bn(7, seed = 11)
  expect_equal(with_backend("cpp", query(bn, c("V7", "V3"), list(V1 = "s1"))),
               with_backend("R", query(bn, c("V7", "V3"), list(V1 = "s1"))), tolerance = 1e-12)
  expect_identical(pgm_backend(), "R")
})
