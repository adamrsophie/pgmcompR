# Structure learning ------------------------------------------------------------------

learning_data <- function(data) {
  if (!is.data.frame(data)) pgm_abort("`data` must be a data frame.")
  if (anyNA(data)) pgm_abort("Structure learning requires complete data.", class = "pgm_missing_data")
  num <- vapply(data, is.numeric, logical(1))
  if (any(num) && !all(num)) pgm_abort("Structure learning supports all discrete or all continuous data.")
  if (!any(num)) data[] <- lapply(data, function(x) if (is.factor(x)) droplevels(x) else factor(x, levels = sort(unique(x))))
  data
}

# Conditional independence tests -------------------------------------------------------

#' Conditional independence tests
#'
#' Tests whether `x` and `y` are independent given `z` in `data`.
#'
#' * `"g2"`: likelihood ratio (G squared) test for discrete data;
#' * `"chi2"`: Pearson chi squared test for discrete data;
#' * `"cor"`: Fisher z test of zero partial correlation for Gaussian data.
#'
#' The discrete tests sum over strata of `z`; degrees of freedom are
#' (rx - 1)(ry - 1) times the number of strata with at least one observation.
#'
#' @param data A data frame.
#' @param x,y Variable names.
#' @param z Conditioning variables.
#' @param test Test name.
#' @return A list of class `pgm_ci_test` with `statistic`, `df` and `p_value`.
#' @examples
#' d <- simulate_bn(example_network("asia"), 2000, seed = 1)
#' ci_test(d, "smoke", "dysp")
#' ci_test(d, "tub", "smoke")
#' @export
ci_test <- function(data, x, y, z = character(), test = c("g2", "chi2", "cor")) {
  test <- match.arg(test)
  if (test == "cor") {
    v <- c(x, y, z)
    C <- stats::cor(as.matrix(data[v]))
    P <- tryCatch(solve(C), error = function(e) pgm_abort("The correlation matrix is singular.", class = "pgm_singular_matrix"))
    r <- -P[1, 2] / sqrt(P[1, 1] * P[2, 2])
    r <- max(min(r, 1 - 1e-15), -1 + 1e-15)
    n <- nrow(data)
    df <- n - length(z) - 3
    if (df <= 0) pgm_abort("Too few observations for the Fisher z test.")
    stat <- sqrt(df) * 0.5 * log((1 + r) / (1 - r))
    return(structure(list(statistic = stat, df = df, p_value = 2 * pnorm(-abs(stat)), test = test), class = "pgm_ci_test"))
  }
  fx <- as.factor(data[[x]]); fy <- as.factor(data[[y]])
  strata <- if (length(z)) interaction(data[z], drop = TRUE) else factor(rep(1, nrow(data)))
  tab <- table(fx, fy, strata)
  stat <- 0
  ns <- 0
  for (k in seq_len(dim(tab)[3])) {
    t2 <- tab[, , k]
    nk <- sum(t2)
    if (nk == 0) next
    ns <- ns + 1
    e <- outer(rowSums(t2), colSums(t2)) / nk
    pos <- t2 > 0
    stat <- stat + if (test == "g2") 2 * sum(t2[pos] * log(t2[pos] / e[pos])) else sum(((t2 - e)^2 / e)[e > 0])
  }
  df <- max(1, (nlevels(fx) - 1) * (nlevels(fy) - 1) * ns)
  structure(list(statistic = stat, df = df, p_value = pchisq(stat, df, lower.tail = FALSE), test = test), class = "pgm_ci_test")
}

#' @export
print.pgm_ci_test <- function(x, ...) {
  cat(sprintf("<CI test %s> statistic = %.4g, df = %g, p = %.4g\n", x$test, x$statistic, x$df, x$p_value))
  invisible(x)
}

# Background knowledge -------------------------------------------------------------------

pair_set <- function(m) {
  if (is.null(m)) return(character())
  m <- parse_edges(m)
  paste(m$from, m$to, sep = "\r")
}

tier_forbidden <- function(tiers, from, to) {
  if (is.null(tiers)) return(FALSE)
  tf <- tier_of(tiers, from); tt <- tier_of(tiers, to)
  !is.na(tf) && !is.na(tt) && tf > tt
}

tier_of <- function(tiers, v) {
  for (i in seq_along(tiers)) if (v %in% tiers[[i]]) return(i)
  NA_integer_
}

allowed_arc <- function(from, to, blacklist, tiers) {
  !(paste(from, to, sep = "\r") %in% blacklist) && !tier_forbidden(tiers, from, to)
}

# PC algorithm ------------------------------------------------------------------------------

#' The PC algorithm
#'
#' Learns a CPDAG by the order independent ("stable") PC algorithm (Spirtes,
#' Glymour and Scheines, 2000; Colombo and Maathuis, 2014): edges are removed
#' when a conditional independence test accepts independence at level
#' `alpha` given a subset of current neighbors, v-structures are oriented
#' from the separating sets, and Meek's rules complete the orientation.
#' Background knowledge: arcs in `blacklist` in both directions remove the
#' edge; one direction forbidden orients the edge the other way; `tiers`
#' (a list of character vectors, earlier tiers first) forbid arcs from later
#' to earlier tiers.
#'
#' Observational data identify a DAG only up to Markov equivalence, so the
#' result is a partially directed graph; undirected edges are not causal
#' claims, and directed ones rest on the faithfulness and causal sufficiency
#' assumptions.
#'
#' @param data A data frame (all discrete or all continuous).
#' @param test CI test (`"g2"`, `"chi2"` or `"cor"`); chosen from the data
#'   type by default.
#' @param alpha Significance level.
#' @param max_cond Largest conditioning set size.
#' @param blacklist,whitelist Arcs (as accepted by [dag()]) to forbid or
#'   require.
#' @param tiers Optional temporal tiers.
#' @return A `pgm_pdag` (or `pgm_dag` if fully oriented) with attribute
#'   `sepsets`.
#' @examples
#' d <- simulate_bn(example_network("sprinkler"), 3000, seed = 1)
#' pc_algorithm(d)
#' @export
pc_algorithm <- function(data, test = NULL, alpha = 0.05, max_cond = Inf, blacklist = NULL,
                         whitelist = NULL, tiers = NULL) {
  data <- learning_data(data)
  test <- test %||% if (is.numeric(data[[1]])) "cor" else "g2"
  vars <- names(data)
  bl <- pair_set(blacklist); wl <- pair_set(whitelist)
  adj <- matrix(TRUE, length(vars), length(vars), dimnames = list(vars, vars))
  diag(adj) <- FALSE
  for (i in vars) for (j in vars) if (i != j) {
    if (!allowed_arc(i, j, bl, tiers) && !allowed_arc(j, i, bl, tiers)) adj[i, j] <- adj[j, i] <- FALSE
  }
  sepset <- list()
  l <- 0
  repeat {
    if (l > max_cond) break
    frozen <- adj
    any_test <- FALSE
    for (i in vars) for (j in vars) {
      if (i >= j || !adj[i, j]) next
      if (paste(i, j, sep = "\r") %in% wl || paste(j, i, sep = "\r") %in% wl) next
      for (side in list(c(i, j), c(j, i))) {
        a <- side[1]; b <- side[2]
        nb <- setdiff(vars[frozen[a, ]], b)
        if (length(nb) < l) next
        any_test <- TRUE
        sets <- if (l == 0) list(character()) else utils::combn(nb, l, simplify = FALSE)
        done <- FALSE
        for (s in sets) {
          p <- ci_test(data, a, b, s, test)$p_value
          if (p > alpha) {
            adj[i, j] <- adj[j, i] <- FALSE
            sepset[[paste(i, j, sep = "\r")]] <- s
            sepset[[paste(j, i, sep = "\r")]] <- s
            done <- TRUE
            break
          }
        }
        if (done) break
      }
    }
    if (!any_test) break
    l <- l + 1
  }
  # Build PDAG with v-structures
  dir <- matrix(FALSE, length(vars), length(vars), dimnames = list(vars, vars))
  und <- adj
  orient <- function(a, b) {
    if (und[a, b] && allowed_arc(a, b, bl, tiers)) {
      und[a, b] <<- und[b, a] <<- FALSE
      dir[a, b] <<- TRUE
    }
  }
  for (w in wl) {
    p <- strsplit(w, "\r", fixed = TRUE)[[1]]
    if (adj[p[1], p[2]]) orient(p[1], p[2])
  }
  for (i in vars) for (j in vars) if (adj[i, j] && und[i, j]) {
    if (!allowed_arc(i, j, bl, tiers) && allowed_arc(j, i, bl, tiers)) orient(j, i)
  }
  for (k in vars) {
    nb <- vars[adj[k, ]]
    if (length(nb) < 2) next
    cmb <- utils::combn(nb, 2)
    for (c in seq_len(ncol(cmb))) {
      i <- cmb[1, c]; j <- cmb[2, c]
      if (adj[i, j]) next
      s <- sepset[[paste(i, j, sep = "\r")]]
      if (!(k %in% s)) {
        if ((und[i, k] || dir[i, k]) && (und[j, k] || dir[j, k]) && !dir[k, i] && !dir[k, j]) {
          if (allowed_arc(i, k, bl, tiers) && allowed_arc(j, k, bl, tiers)) {
            orient(i, k); orient(j, k)
          }
        }
      }
    }
  }
  from <- to <- type <- character()
  for (i in seq_along(vars)) for (j in seq_along(vars)) {
    if (dir[i, j]) { from <- c(from, vars[i]); to <- c(to, vars[j]); type <- c(type, "directed") }
    else if (und[i, j] && i < j) { from <- c(from, vars[i]); to <- c(to, vars[j]); type <- c(type, "undirected") }
  }
  p <- new_pgm_graph(vars, edge_frame(from, to, type), "pgm_pdag")
  if (has_directed_cycle(p)) {
    pgm_warn("Conflicting v-structures produced a directed cycle; the result is not a valid PDAG.", class = "pgm_learning_warning")
  } else p <- meek_rules(p)
  attr(p, "sepsets") <- sepset
  p
}

# Decomposable scores -------------------------------------------------------------------------

local_counts <- function(node, parents, data) {
  x <- data[[node]]
  r <- nlevels(x)
  if (!length(parents)) return(matrix(tabulate(as.integer(x), r), r))
  cfg <- interaction(data[parents], drop = FALSE, lex.order = FALSE)
  t <- table(x, cfg)
  matrix(as.numeric(t), r)
}

local_score_loglik <- function(node, parents, data, ...) {
  if (is.numeric(data[[node]])) return(gaussian_local(node, parents, data, 0))
  n <- local_counts(node, parents, data)
  cs <- colSums(n)
  pos <- n > 0
  sum(n[pos] * log((n / rep(cs, each = nrow(n)))[pos]))
}

n_params <- function(node, parents, data) {
  if (is.numeric(data[[node]])) return(length(parents) + 2)
  r <- nlevels(data[[node]])
  q <- prod(vapply(parents, function(p) nlevels(data[[p]]), integer(1)))
  (r - 1) * q
}

local_score_bic <- function(node, parents, data, ...) {
  local_score_loglik(node, parents, data) - 0.5 * log(nrow(data)) * n_params(node, parents, data)
}

local_score_aic <- function(node, parents, data, ...) {
  local_score_loglik(node, parents, data) - n_params(node, parents, data)
}

local_score_bic_g <- function(node, parents, data, ...) local_score_bic(node, parents, data)

gaussian_local <- function(node, parents, data, penalty) {
  y <- data[[node]]
  X <- cbind(1, as.matrix(data[parents]))
  fit <- lm.fit(X, y)
  n <- length(y)
  s2 <- sum(fit$residuals^2) / n
  if (!(s2 > 0)) return(Inf)
  -0.5 * n * (log(2 * pi * s2) + 1)
}

bd_score <- function(n, alpha) {
  a_j <- colSums(alpha)
  n_j <- colSums(n)
  sum(lgamma(a_j) - lgamma(a_j + n_j)) + sum(lgamma(alpha + n) - lgamma(alpha))
}

local_score_bdeu <- function(node, parents, data, iss = 1, ...) {
  if (is.numeric(data[[node]])) pgm_abort("BDeu needs discrete data.")
  n <- local_counts(node, parents, data)
  bd_score(n, matrix(iss / length(n), nrow(n), ncol(n)))
}

local_score_k2 <- function(node, parents, data, ...) {
  if (is.numeric(data[[node]])) pgm_abort("K2 needs discrete data.")
  n <- local_counts(node, parents, data)
  bd_score(n, matrix(1, nrow(n), ncol(n)))
}

#' Network scores
#'
#' Computes a decomposable score of a DAG for `data`: the sum over nodes of a
#' local score of the node given its parents. Built in scores: `"loglik"`,
#' `"aic"`, `"bic"` (discrete or Gaussian data), `"bdeu"` (with imaginary
#' sample size `iss`) and `"k2"` (discrete data). Larger is better. Custom
#' scores can be registered with [register_score()].
#'
#' @param graph A DAG or Bayesian network.
#' @param data A complete data frame.
#' @param score Score name.
#' @param ... Passed to the local score (for example `iss`).
#' @return A numeric scalar.
#' @examples
#' bn <- example_network("asia")
#' d <- simulate_bn(bn, 1000, seed = 1)
#' network_score(bn, d, "bic")
#' @export
network_score <- function(graph, data, score = "bic", ...) {
  g <- as_pgm_graph(graph)
  data <- learning_data(data[g$nodes])
  fun <- registry_get("score", score)
  sum(vapply(g$nodes, function(v) fun(v, parents(g, v), data, ...), numeric(1)))
}

# Score based search ------------------------------------------------------------------------

reaches <- function(A, from, to) {
  # Is `to` reachable from `from` along arcs of adjacency matrix A?
  seen <- from
  fr <- from
  while (length(fr)) {
    nx <- setdiff(which(colSums(A[fr, , drop = FALSE]) > 0), seen)
    if (to %in% nx) return(TRUE)
    seen <- c(seen, nx)
    fr <- nx
  }
  FALSE
}

#' Score based structure learning
#'
#' `hill_climb()` performs greedy search over DAGs with arc addition,
#' deletion and reversal, moving to the best scoring neighbor until no move
#' improves the score. `tabu_search()` continues past local optima, allowing
#' the best non tabu move even when it lowers the score and keeping the best
#' DAG seen; recently visited DAGs are tabu. `mmhc()` restricts the search to
#' candidate neighbors found by the max min parents and children algorithm
#' (Tsamardinos, Brown and Aliferis, 2006). Local scores are cached.
#'
#' @param data A complete data frame (all discrete or all continuous).
#' @param score Score name (see [network_score()]).
#' @param start Optional starting DAG.
#' @param blacklist,whitelist Arcs to forbid or require.
#' @param tiers Optional temporal tiers (no arc from a later to an earlier
#'   tier).
#' @param max_parents Maximum in-degree.
#' @param max_iter Maximum number of moves.
#' @param tabu_length Number of recent DAGs kept tabu.
#' @param max_no_improve Moves without improvement before `tabu_search()`
#'   stops.
#' @param candidates Optional named list of allowed neighbors per node.
#' @param alpha,test,max_cond Settings for the MMPC phase of `mmhc()`.
#' @param ... Passed to the score.
#' @return A `pgm_dag` with attribute `score`.
#' @examples
#' d <- simulate_bn(example_network("asia"), 3000, seed = 1)
#' hill_climb(d)
#' @export
hill_climb <- function(data, score = "bic", start = NULL, blacklist = NULL, whitelist = NULL, tiers = NULL,
                       max_parents = Inf, max_iter = 1000L, candidates = NULL, ...) {
  score_search(data, score, start, blacklist, whitelist, tiers, max_parents, max_iter, 0L, 0L, candidates, ...)
}

#' @rdname hill_climb
#' @export
tabu_search <- function(data, score = "bic", start = NULL, blacklist = NULL, whitelist = NULL, tiers = NULL,
                        max_parents = Inf, max_iter = 1000L, tabu_length = 10L, max_no_improve = 10L, ...) {
  score_search(data, score, start, blacklist, whitelist, tiers, max_parents, max_iter, tabu_length, max_no_improve, NULL, ...)
}

score_search <- function(data, score, start, blacklist, whitelist, tiers, max_parents, max_iter,
                         tabu_length, max_no_improve, candidates, ...) {
  data <- learning_data(data)
  vars <- names(data)
  p <- length(vars)
  fun <- registry_get("score", score)
  cache <- new.env(hash = TRUE)
  local <- function(j, pa) {
    key <- paste(j, paste(sort(pa), collapse = ","))
    val <- cache[[key]]
    if (is.null(val)) {
      val <- fun(vars[j], vars[sort(pa)], data, ...)
      assign(key, val, envir = cache)
    }
    val
  }
  bl <- pair_set(blacklist); wl <- pair_set(whitelist)
  allowed <- matrix(TRUE, p, p); diag(allowed) <- FALSE
  for (i in seq_len(p)) for (j in seq_len(p)) if (i != j) {
    allowed[i, j] <- allowed_arc(vars[i], vars[j], bl, tiers)
    if (!is.null(candidates) && !(vars[i] %in% candidates[[vars[j]]] || vars[j] %in% candidates[[vars[i]]])) allowed[i, j] <- FALSE
  }
  required <- matrix(FALSE, p, p)
  for (w in wl) {
    e <- strsplit(w, "\r", fixed = TRUE)[[1]]
    required[match(e[1], vars), match(e[2], vars)] <- TRUE
  }
  A <- matrix(0L, p, p, dimnames = list(vars, vars))
  if (!is.null(start)) {
    sg <- as_pgm_graph(start)
    A[cbind(match(sg$edges$from, vars), match(sg$edges$to, vars))] <- 1L
  }
  A[required] <- 1L
  if (has_directed_cycle(graph_from_adjacency(A))) pgm_abort("The whitelist or start graph contains a directed cycle.")
  node_score <- vapply(seq_len(p), function(j) local(j, which(A[, j] == 1L)), numeric(1))
  best_A <- A; best_total <- sum(node_score)
  tabu <- character()
  no_improve <- 0L
  for (it in seq_len(max_iter)) {
    best_delta <- -Inf; best_move <- NULL
    keyA <- function(M) paste(which(M == 1L), collapse = ",")
    for (i in seq_len(p)) for (j in seq_len(p)) {
      if (i == j) next
      if (A[i, j] == 1L) {
        # deletion
        if (!required[i, j]) {
          pa <- setdiff(which(A[, j] == 1L), i)
          d <- local(j, pa) - node_score[j]
          if (d > best_delta) {
            M <- A; M[i, j] <- 0L
            if (!(keyA(M) %in% tabu)) { best_delta <- d; best_move <- list("del", i, j, d, NULL) }
          }
        }
        # reversal
        if (!required[i, j] && allowed[j, i] && sum(A[, i]) < max_parents) {
          M <- A; M[i, j] <- 0L
          if (!reaches(M, i, j)) {
            paj <- setdiff(which(A[, j] == 1L), i)
            pai <- c(which(A[, i] == 1L), j)
            d <- local(j, paj) - node_score[j] + local(i, pai) - node_score[i]
            if (d > best_delta) {
              M[j, i] <- 1L
              if (!(keyA(M) %in% tabu)) { best_delta <- d; best_move <- list("rev", i, j, d, NULL) }
            }
          }
        }
      } else if (A[j, i] == 0L && allowed[i, j] && sum(A[, j]) < max_parents) {
        if (!reaches(A, j, i)) {
          pa <- c(which(A[, j] == 1L), i)
          d <- local(j, pa) - node_score[j]
          if (d > best_delta) {
            M <- A; M[i, j] <- 1L
            if (!(keyA(M) %in% tabu)) { best_delta <- d; best_move <- list("add", i, j, d, NULL) }
          }
        }
      }
    }
    if (is.null(best_move)) break
    if (best_delta <= 1e-10 && tabu_length == 0L) break
    i <- best_move[[2]]; j <- best_move[[3]]
    if (best_move[[1]] == "add") A[i, j] <- 1L
    if (best_move[[1]] == "del") A[i, j] <- 0L
    if (best_move[[1]] == "rev") { A[i, j] <- 0L; A[j, i] <- 1L }
    node_score[j] <- local(j, which(A[, j] == 1L))
    node_score[i] <- local(i, which(A[, i] == 1L))
    total <- sum(node_score)
    if (tabu_length > 0L) {
      tabu <- utils::tail(c(tabu, paste(which(A == 1L), collapse = ",")), tabu_length)
      if (total > best_total + 1e-10) { best_total <- total; best_A <- A; no_improve <- 0L }
      else { no_improve <- no_improve + 1L; if (no_improve >= max_no_improve) break }
    } else {
      best_A <- A; best_total <- total
    }
  }
  g <- graph_from_adjacency(best_A)
  g <- dag(g$edges[, c("from", "to")], nodes = vars)
  attr(g, "score") <- best_total
  g
}

# Max min parents and children (MMPC) for one target.
mmpc <- function(data, target, test, alpha, max_cond) {
  vars <- setdiff(names(data), target)
  cpc <- character()
  assoc_min <- function(x, cpc) {
    sets <- list(character())
    k <- min(length(cpc), max_cond)
    if (k > 0) for (s in seq_len(k)) sets <- c(sets, utils::combn(cpc, s, simplify = FALSE))
    max(vapply(sets, function(s) ci_test(data, target, x, s, test)$p_value, numeric(1)))
  }
  repeat {
    rest <- setdiff(vars, cpc)
    if (!length(rest)) break
    pv <- vapply(rest, assoc_min, numeric(1), cpc = cpc)
    best <- which.min(pv)
    if (pv[best] > alpha) break
    cpc <- c(cpc, rest[best])
  }
  for (x in cpc) {
    others <- setdiff(cpc, x)
    if (assoc_min(x, others) > alpha) cpc <- setdiff(cpc, x)
  }
  cpc
}

#' @rdname hill_climb
#' @export
mmhc <- function(data, score = "bic", alpha = 0.05, test = NULL, max_cond = 3L, blacklist = NULL,
                 whitelist = NULL, tiers = NULL, max_parents = Inf, ...) {
  data <- learning_data(data)
  test <- test %||% if (is.numeric(data[[1]])) "cor" else "g2"
  vars <- names(data)
  cpc <- lapply(setNames(vars, vars), function(v) mmpc(data, v, test, alpha, max_cond))
  # Symmetry correction: keep x in CPC(t) only if t is in CPC(x).
  cand <- lapply(setNames(vars, vars), function(v) Filter(function(x) v %in% cpc[[x]], cpc[[v]]))
  for (w in pair_set(whitelist)) {
    e <- strsplit(w, "\r", fixed = TRUE)[[1]]
    cand[[e[2]]] <- union(cand[[e[2]]], e[1]); cand[[e[1]]] <- union(cand[[e[1]]], e[2])
  }
  g <- score_search(data, score, NULL, blacklist, whitelist, tiers, max_parents, 1000L, 0L, 0L, cand, ...)
  attr(g, "candidates") <- cand
  g
}

#' Structure learning front end
#'
#' Dispatches to the native learners or, with `engine = "bnlearn"`, to the
#' corresponding 'bnlearn' function (`hc`, `tabu`, `pc.stable`, `mmhc`),
#' returning the learned graph as a pgmcompR graph.
#'
#' @param data A complete data frame.
#' @param algorithm `"hill_climb"`, `"tabu"`, `"pc"` or `"mmhc"`.
#' @param engine `"native"` or `"bnlearn"`.
#' @param ... Passed to the learner.
#' @return A `pgm_dag` or `pgm_pdag`.
#' @export
learn_structure <- function(data, algorithm = c("hill_climb", "tabu", "pc", "mmhc"), engine = c("native", "bnlearn"), ...) {
  algorithm <- match.arg(algorithm)
  engine <- match.arg(engine)
  if (engine == "bnlearn") {
    if (!requireNamespace("bnlearn", quietly = TRUE)) pgm_abort("The 'bnlearn' package is required for engine = 'bnlearn'.")
    data <- learning_data(data)
    fun <- switch(algorithm, hill_climb = bnlearn::hc, tabu = bnlearn::tabu, pc = bnlearn::pc.stable, mmhc = bnlearn::mmhc)
    return(bnlearn_graph(fun(data, ...)))
  }
  switch(algorithm, hill_climb = hill_climb(data, ...), tabu = tabu_search(data, ...),
         pc = pc_algorithm(data, ...), mmhc = mmhc(data, ...))
}

#' Bootstrap edge stability
#'
#' Relearns the structure on `R` nonparametric bootstrap samples and reports,
#' for each pair of variables, the proportion of networks containing an edge
#' (`strength`) and, among those, the proportion oriented `from -> to`
#' (`direction`). `averaged_network()` keeps edges with strength at least
#' `threshold`, oriented by majority, skipping arcs that would create cycles.
#'
#' @param data A complete data frame.
#' @param R Number of bootstrap samples.
#' @param algorithm Learner passed to [learn_structure()].
#' @param seed Optional seed.
#' @param ... Passed to the learner.
#' @param strength Output of `bootstrap_structure()`.
#' @param threshold Minimum strength.
#' @return `bootstrap_structure()` returns a data frame with columns `from`,
#'   `to`, `strength`, `direction`; `averaged_network()` a `pgm_dag`.
#' @examples
#' d <- simulate_bn(example_network("sprinkler"), 500, seed = 1)
#' s <- bootstrap_structure(d, R = 10, seed = 1)
#' averaged_network(s)
#' @export
bootstrap_structure <- function(data, R = 100L, algorithm = "hill_climb", seed = NULL, ...) {
  vars <- names(data)
  cnt_dir <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  cnt_und <- cnt_dir
  with_seed(seed, {
    for (r in seq_len(R)) {
      b <- data[sample.int(nrow(data), replace = TRUE), , drop = FALSE]
      g <- learn_structure(b, algorithm = algorithm, ...)
      e <- g$edges
      for (k in seq_len(nrow(e))) {
        if (e$type[k] == "directed") cnt_dir[e$from[k], e$to[k]] <- cnt_dir[e$from[k], e$to[k]] + 1
        else {
          cnt_und[e$from[k], e$to[k]] <- cnt_und[e$from[k], e$to[k]] + 1
          cnt_und[e$to[k], e$from[k]] <- cnt_und[e$to[k], e$from[k]] + 1
        }
      }
    }
  })
  out <- list()
  for (i in seq_along(vars)) for (j in seq_along(vars)) {
    if (i >= j) next
    tot <- cnt_dir[i, j] + cnt_dir[j, i] + cnt_und[i, j]
    if (tot == 0) next
    dij <- (cnt_dir[i, j] + cnt_und[i, j] / 2) / tot
    out[[length(out) + 1L]] <- data.frame(from = c(vars[i], vars[j]), to = c(vars[j], vars[i]),
                                          strength = tot / R, direction = c(dij, 1 - dij), stringsAsFactors = FALSE)
  }
  res <- if (length(out)) do.call(rbind, out) else data.frame(from = character(), to = character(), strength = numeric(), direction = numeric())
  attr(res, "nodes") <- vars
  attr(res, "R") <- R
  res
}

#' @rdname bootstrap_structure
#' @export
averaged_network <- function(strength, threshold = 0.5) {
  vars <- attr(strength, "nodes") %||% unique(c(strength$from, strength$to))
  s <- strength[strength$strength >= threshold & strength$direction > 0.5, , drop = FALSE]
  s <- s[order(-s$strength, -s$direction), , drop = FALSE]
  g <- dag(nodes = vars)
  for (k in seq_len(nrow(s))) {
    g2 <- tryCatch(add_edges(g, s$from[k], s$to[k]), pgm_cycle_error = function(e) NULL)
    if (!is.null(g2)) g <- g2
  }
  g
}
