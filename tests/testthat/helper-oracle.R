# Independent brute force oracle: enumerates every complete assignment and
# multiplies CPT entries looked up directly from the stored column major
# values. It does not use the factor engine.
oracle_joint <- function(bn) {
  vars <- bn$graph$nodes
  st <- lapply(bn$registry[vars], `[[`, "states")
  grid <- expand.grid(st, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  p <- rep(1, nrow(grid))
  for (v in vars) {
    x <- bn$cpds[[v]]
    fam <- c(x$child, x$parents)
    idx <- rep(1, nrow(grid))
    stride <- 1
    for (u in fam) {
      idx <- idx + (match(grid[[u]], x$states[[u]]) - 1) * stride
      stride <- stride * length(x$states[[u]])
    }
    p <- p * x$values[idx]
  }
  grid$p <- p
  grid
}

oracle_query <- function(bn, variables, evidence = list()) {
  j <- oracle_joint(bn)
  keep <- rep(TRUE, nrow(j))
  for (v in names(evidence)) keep <- keep & j[[v]] == evidence[[v]]
  j <- j[keep, , drop = FALSE]
  st <- lapply(bn$registry[variables], `[[`, "states")
  grid <- expand.grid(st, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  vals <- vapply(seq_len(nrow(grid)), function(i) {
    sel <- rep(TRUE, nrow(j))
    for (v in variables) sel <- sel & j[[v]] == grid[[v]][i]
    sum(j$p[sel])
  }, numeric(1))
  vals / sum(vals)
}

random_bn <- function(n = 6, max_card = 3, max_parents = 2, seed = 1, zeros = FALSE) {
  set.seed(seed)
  vars <- paste0("V", seq_len(n))
  edges <- list()
  for (i in seq_len(n)[-1]) {
    k <- sample(0:min(max_parents, i - 1), 1)
    if (k) for (p in sample(vars[seq_len(i - 1)], k)) edges[[length(edges) + 1]] <- c(p, vars[i])
  }
  card <- setNames(sample(2:max_card, n, replace = TRUE), vars)
  st <- lapply(card, function(k) paste0("s", seq_len(k)))
  bn <- bayesian_network(edges, nodes = vars, states = st)
  for (v in vars) {
    pa <- parents(bn, v)
    q <- prod(card[pa])
    m <- matrix(rgamma(card[[v]] * q, 1), card[[v]])
    if (zeros) m[sample(length(m), floor(length(m) / 4))] <- 0
    m[, colSums(m) == 0] <- 1
    m <- sweep(m, 2, colSums(m), "/")
    bn <- add_cpd(bn, cpd(v, as.numeric(m), parents = pa, states = st[c(v, pa)]))
  }
  bn
}
