# Benchmarks of pgmcompR's R and C++ backends, with bnlearn and gRain for
# reference. Run from the package root with
#   Rscript inst/benchmarks/run-benchmarks.R
# The results are saved to inst/benchmarks/results.rds and shown in the
# "Performance and backends" vignette, which does not rerun them.

library(pgmcompR)
stopifnot(requireNamespace("bench", quietly = TRUE))

random_chain_net <- function(n, card, seed) {
  set.seed(seed)
  vars <- sprintf("V%03d", seq_len(n))
  st <- setNames(rep(list(paste0("s", seq_len(card))), n), vars)
  edges <- list()
  for (i in seq_len(n)[-1]) {
    for (p in unique(c(i - 1, sample(max(1, i - 4):(i - 1), 1)))) edges[[length(edges) + 1]] <- c(vars[p], vars[i])
  }
  bn <- bayesian_network(edges, nodes = vars, states = st)
  for (v in vars) {
    pa <- pgmcompR::parents(bn, v)
    m <- matrix(rgamma(card * card^length(pa), 1), card)
    bn <- add_cpd(bn, cpd(v, as.numeric(sweep(m, 2, colSums(m), "/")), pa, st[c(v, pa)]))
  }
  bn
}

results <- list()

# 1. Factor product and marginalization kernels
for (k in c(4, 8, 12)) {
  st <- lapply(setNames(nm = c("A", "B", "C", "D")), function(v) paste0("s", seq_len(k)))
  f <- discrete_factor(c("A", "B", "C"), st[c("A", "B", "C")], runif(k^3))
  g <- discrete_factor(c("C", "D", "A"), st[c("C", "D", "A")], runif(k^3))
  b <- bench::mark(
    R = with_backend("R", factor_marginalize(factor_product(f, g), "C")),
    cpp = with_backend("cpp", factor_marginalize(factor_product(f, g), "C")),
    check = function(a, b) isTRUE(all.equal(factor_values(a), factor_values(b), tolerance = 1e-12)),
    min_iterations = 20
  )
  results[[length(results) + 1]] <- data.frame(task = "product + marginalize", size = sprintf("%d^4 cells", k),
                                               engine = as.character(b$expression), median_ms = as.numeric(b$median) * 1000)
}

# 2. Exact inference on random networks
for (n in c(20, 50, 100)) {
  bn <- random_chain_net(n, 3, seed = n)
  vars <- pgmcompR::nodes(bn)
  ev <- setNames(list("s1", "s2"), vars[c(n, n %/% 2)])
  marks <- list(
    ve_R = function() with_backend("R", query(bn, vars[1], ev)),
    ve_cpp = function() with_backend("cpp", query(bn, vars[1], ev)),
    jt_R = function() with_backend("R", query(bn, vars[1], ev, method = "jt"))
  )
  if (requireNamespace("gRain", quietly = TRUE)) {
    gr <- as_grain(bn)
    marks$gRain <- function() gRain::querygrain(gRain::setEvidence(gr, evidence = ev), vars[1])
  }
  if (requireNamespace("bnlearn", quietly = TRUE)) {
    fit <- as_bnlearn(bn)
    marks$bnlearn_lw <- function() bnlearn::cpdist(fit, nodes = vars[1], evidence = ev, method = "lw", n = 5000)
  }
  b <- bench::mark(exprs = lapply(marks, function(f) bquote(.(f)())), check = FALSE, min_iterations = 5)
  results[[length(results) + 1]] <- data.frame(task = "query, 2 evidence nodes", size = sprintf("%d nodes", n),
                                               engine = names(marks), median_ms = as.numeric(b$median) * 1000)
}

# 3. Structure learning
d <- simulate_bn(example_network("asia"), 5000, seed = 1)
marks <- list(native_hc = function() hill_climb(d))
if (requireNamespace("bnlearn", quietly = TRUE)) marks$bnlearn_hc <- function() bnlearn::hc(d)
b <- bench::mark(exprs = lapply(marks, function(f) bquote(.(f)())), check = FALSE, min_iterations = 3)
results[[length(results) + 1]] <- data.frame(task = "hill climbing, Asia", size = "5000 rows",
                                             engine = names(marks), median_ms = as.numeric(b$median) * 1000)

out <- do.call(rbind, results)
attr(out, "session") <- list(R = R.version.string, platform = R.version$platform, date = Sys.Date(),
                             pgmcompR = as.character(utils::packageVersion("pgmcompR")))
saveRDS(out, file.path("inst", "benchmarks", "results.rds"))
print(out)
