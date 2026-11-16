# Junction trees -------------------------------------------------------------------

#' Junction trees
#'
#' Compiles a discrete model into a junction tree (clique tree): the
#' interaction graph (the moral graph for a Bayesian network) is
#' triangulated, its maximal cliques are joined by a maximum weight spanning
#' tree on separator sizes, which satisfies the running intersection
#' property, and each factor is assigned to one clique containing its scope.
#' `calibrate()` enters evidence and runs the two pass Shafer and Shenoy
#' sum product schedule so that every clique belief is proportional to the
#' joint posterior of the clique's variables.
#'
#' @param model A discrete Bayesian network, Markov network, factor graph or
#'   list of factors.
#' @param heuristic Triangulation heuristic (see [triangulate()]).
#' @param order Optional explicit elimination order.
#' @param jt A `junction_tree`.
#' @param evidence Named list of observed states.
#' @param virtual_evidence Named list of likelihood vectors.
#' @return `junction_tree()` and `calibrate()` return a `junction_tree`;
#'   `clique_beliefs()` returns the list of normalized clique beliefs.
#' @examples
#' jt <- junction_tree(example_network("asia"))
#' jt
#' jt <- calibrate(jt, evidence = list(xray = "yes"))
#' clique_beliefs(jt)[[1]]
#' @export
junction_tree <- function(model, heuristic = "min_fill", order = NULL) {
  factors <- unname(model_factors(model))
  reg <- factor_registry(factors)
  states <- if (!is.null(model$registry)) model_states(model) else lapply(reg, `[[`, "states")
  vars <- names(states)
  adj <- interaction_adjacency(factors, vars)
  pairs_from <- unlist(lapply(names(adj), function(v) rep(v, length(adj[[v]]))), use.names = FALSE)
  pairs_to <- unlist(adj, use.names = FALSE)
  ug <- undirected_from_pairs(vars, pairs_from, pairs_to)
  card <- vapply(states, length, integer(1))
  tri <- triangulate(ug, order = order, heuristic = heuristic, cardinality = card)
  cliques <- tri$cliques
  k <- length(cliques)
  # Maximum weight spanning forest (Kruskal), weight = separator size.
  tree <- matrix(character(), 0, 2)
  if (k > 1L) {
    cand <- t(combn(k, 2L))
    w <- apply(cand, 1L, function(p) length(intersect(cliques[[p[1]]], cliques[[p[2]]])))
    cand <- cand[w > 0, , drop = FALSE]
    w <- w[w > 0]
    ordw <- order(-w, cand[, 1], cand[, 2])
    comp <- seq_len(k)
    edges_t <- list()
    for (i in ordw) {
      a <- cand[i, 1]; b <- cand[i, 2]
      if (comp[a] != comp[b]) {
        old <- comp[b]
        comp[comp == old] <- comp[a]
        edges_t[[length(edges_t) + 1L]] <- c(a, b)
      }
    }
    tree <- if (length(edges_t)) do.call(rbind, edges_t) else matrix(integer(), 0, 2)
  } else tree <- matrix(integer(), 0, 2)
  # Assign factors to cliques.
  assign <- vapply(factors, function(f) {
    hit <- which(vapply(cliques, function(cl) all(f$vars %in% cl), logical(1)))
    if (!length(hit)) pgm_abort("Internal error: a factor scope is not covered by any clique.")
    hit[which.min(vapply(cliques[hit], length, integer(1)))]
  }, numeric(1))
  potentials <- lapply(seq_len(k), function(i) {
    fs <- factors[assign == i]
    unit <- new_factor(cliques[[i]], states[cliques[[i]]], rep(1, prod(card[cliques[[i]]])))
    factor_reorder(factor_product(c(list(unit), fs)), cliques[[i]])
  })
  structure(list(cliques = cliques, tree = tree, potentials = potentials, states = states,
                 factors = factors, registry = lapply(states, function(s) list(type = "discrete", states = s)),
                 beliefs = NULL, evidence = list(), order = tri$order),
            class = c("junction_tree", "pgm_model"))
}

#' @export
nodes.junction_tree <- function(x) names(x$states)

#' @export
model_factors.junction_tree <- function(model) model$factors

#' @export
print.junction_tree <- function(x, ...) {
  cat(sprintf("<Junction tree: %d cliques, max clique size %d, %s>\n", length(x$cliques),
              max(vapply(x$cliques, length, integer(1))), if (is.null(x$beliefs)) "not calibrated" else "calibrated"))
  for (i in seq_along(x$cliques)) cat(sprintf("  C%d: {%s}\n", i, paste(x$cliques[[i]], collapse = ", ")))
  if (nrow(x$tree)) {
    sep <- apply(x$tree, 1L, function(p) paste(intersect(x$cliques[[p[1]]], x$cliques[[p[2]]]), collapse = ", "))
    cat(paste0("  C", x$tree[, 1], " -- C", x$tree[, 2], "  sep {", sep, "}", collapse = "\n"), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname junction_tree
#' @export
calibrate <- function(jt, evidence = list(), virtual_evidence = NULL) {
  if (!inherits(jt, "junction_tree")) pgm_abort("`jt` must be a junction_tree.")
  evidence <- normalize_evidence(evidence)
  for (v in names(evidence)) {
    if (is.null(jt$states[[v]])) pgm_abort(sprintf("Unknown variable '%s'.", v), class = "pgm_unknown_variable")
    if (!(evidence[[v]] %in% jt$states[[v]])) pgm_abort(sprintf("'%s' is not a state of '%s'.", evidence[[v]], v), class = "pgm_unknown_state")
  }
  k <- length(jt$cliques)
  pot <- jt$potentials
  home <- function(v) which(vapply(jt$cliques, function(cl) v %in% cl, logical(1)))[1L]
  for (v in names(evidence)) {
    i <- home(v)
    ind <- discrete_factor(v, jt$states[v], as.numeric(jt$states[[v]] == evidence[[v]]))
    pot[[i]] <- factor_product(pot[[i]], ind)
  }
  for (f in virtual_evidence_factors(virtual_evidence, jt$states)) {
    i <- home(f$vars)
    pot[[i]] <- factor_product(pot[[i]], f)
  }
  nb <- lapply(seq_len(k), function(i) c(jt$tree[jt$tree[, 1] == i, 2], jt$tree[jt$tree[, 2] == i, 1]))
  msg <- list()
  key <- function(i, j) paste(i, j)
  send <- function(i, j) {
    incoming <- lapply(setdiff(nb[[i]], j), function(m) msg[[key(m, i)]])
    prodf <- factor_product(c(list(pot[[i]]), incoming))
    sepv <- intersect(jt$cliques[[i]], jt$cliques[[j]])
    out <- factor_marginalize(prodf, setdiff(prodf$vars, sepv))
    s <- sum(out$values)
    if (s > 0) out$values <- out$values / s
    out
  }
  visited <- rep(FALSE, k)
  for (root in seq_len(k)) {
    if (visited[root]) next
    # Order nodes of this component by BFS from root.
    ord <- root
    par <- setNames(NA_integer_, root)
    q <- root
    visited[root] <- TRUE
    while (length(q)) {
      i <- q[1]; q <- q[-1]
      for (j in nb[[i]]) if (!visited[j]) {
        visited[j] <- TRUE
        ord <- c(ord, j)
        par[as.character(j)] <- i
        q <- c(q, j)
      }
    }
    for (j in rev(ord[-1])) {
      p <- par[[as.character(j)]]
      msg[[key(j, p)]] <- send(j, p)
    }
    for (j in ord[-1]) {
      p <- par[[as.character(j)]]
      msg[[key(p, j)]] <- send(p, j)
    }
  }
  beliefs <- lapply(seq_len(k), function(i) {
    b <- factor_product(c(list(pot[[i]]), lapply(nb[[i]], function(m) msg[[key(m, i)]])))
    if (!(sum(b$values) > 0)) impossible_evidence(evidence)
    factor_reorder(factor_normalize(b), jt$cliques[[i]])
  })
  jt$beliefs <- beliefs
  jt$evidence <- evidence
  jt$virtual_evidence <- virtual_evidence
  jt
}

#' @rdname junction_tree
#' @export
clique_beliefs <- function(jt) {
  if (is.null(jt$beliefs)) pgm_abort("The junction tree is not calibrated; call `calibrate()` first.")
  jt$beliefs
}

infer_jt <- function(model, variables, evidence, joint = TRUE, virtual_evidence = NULL,
                     interventions = NULL, elimination_order = "min_fill", log = FALSE, ...) {
  if (length(normalize_evidence(interventions, "interventions"))) {
    model <- intervene(model, normalize_evidence(interventions, "interventions"))
  }
  jt <- if (inherits(model, "junction_tree")) model else junction_tree(model)
  prep <- prepare_discrete(jt, variables, evidence, virtual_evidence, prune = FALSE)
  same_ev <- !is.null(jt$beliefs) && identical(jt$evidence, prep$evidence) && is.null(virtual_evidence) && is.null(jt$virtual_evidence)
  if (!same_ev) jt <- calibrate(jt, prep$evidence, virtual_evidence)
  hit <- which(vapply(jt$cliques, function(cl) all(prep$variables %in% cl), logical(1)))
  if (length(hit)) {
    i <- hit[which.min(vapply(jt$cliques[hit], length, integer(1)))]
    b <- jt$beliefs[[i]]
    post <- factor_reorder(factor_marginalize(b, setdiff(b$vars, prep$variables)), prep$variables)
    post <- if (log) factor_to_log(post) else post
    return(finish_posterior(post, prep, log, "jt", joint))
  }
  # Variables span several cliques: the calibrated beliefs and separator
  # marginals define the joint; eliminate the other variables exactly.
  fs <- jt$beliefs
  seps <- lapply(seq_len(nrow(jt$tree)), function(r) {
    a <- jt$tree[r, 1]; b <- jt$tree[r, 2]
    s <- intersect(jt$cliques[[a]], jt$cliques[[b]])
    m <- factor_marginalize(jt$beliefs[[a]], setdiff(jt$beliefs[[a]]$vars, s))
    m
  })
  pool <- fs
  for (s in seps) {
    i <- which(vapply(pool, function(f) all(s$vars %in% f$vars), logical(1)))[1L]
    pool[[i]] <- factor_divide(pool[[i]], s)
  }
  scope <- unique(unlist(lapply(pool, `[[`, "vars"), use.names = FALSE))
  ord <- resolve_order(pool, setdiff(scope, prep$variables), elimination_order)
  res <- eliminate_vars(pool, ord, "sum")
  j <- factor_reorder(factor_product(res$factors), prep$variables)
  if (log) j <- factor_to_log(j)
  finish_posterior(j, prep, log, "jt", joint)
}
