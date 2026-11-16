# Interoperability ---------------------------------------------------------------------

#' Convert to and from igraph
#'
#' `as_igraph()` converts a graph or model to an 'igraph' object (directed
#' when every edge is directed; a mixed graph becomes a directed graph with
#' undirected edges represented in both directions and marked by the edge
#' attribute `type`). `as_pgm()` converts an 'igraph' graph back: directed
#' acyclic graphs become DAGs, undirected graphs undirected graphs.
#'
#' @param x A pgmcompR graph or model, or an 'igraph' object.
#' @return An 'igraph' object or a `pgm_graph`.
#' @export
as_igraph <- function(x) {
  if (!requireNamespace("igraph", quietly = TRUE)) pgm_abort("The 'igraph' package is required.")
  g <- as_pgm_graph(x)
  e <- g$edges
  if (all(e$type == "undirected")) {
    out <- igraph::graph_from_data_frame(e[, c("from", "to")], directed = FALSE, vertices = data.frame(name = g$nodes))
  } else {
    u <- e$type == "undirected"
    ee <- rbind(e, data.frame(from = e$to[u], to = e$from[u], type = e$type[u], stringsAsFactors = FALSE))
    out <- igraph::graph_from_data_frame(ee, directed = TRUE, vertices = data.frame(name = g$nodes))
  }
  out
}

#' @rdname as_igraph
#' @export
as_pgm <- function(x) {
  if (inherits(x, "igraph")) {
    if (!requireNamespace("igraph", quietly = TRUE)) pgm_abort("The 'igraph' package is required.")
    el <- igraph::as_edgelist(x, names = TRUE)
    nm <- igraph::V(x)$name %||% as.character(seq_len(igraph::vcount(x)))
    if (!igraph::is_directed(x)) return(undirected_graph(if (nrow(el)) el else NULL, nodes = nm))
    return(pgm_graph(if (nrow(el)) el else NULL, nodes = nm))
  }
  if (inherits(x, "bn") || inherits(x, "bn.fit")) return(from_bnlearn(x))
  if (inherits(x, "grain")) return(from_grain(x))
  pgm_abort("Cannot convert this object to a pgmcompR model.")
}

bnlearn_graph <- function(x) {
  arcs <- bnlearn::arcs(x)
  nm <- bnlearn::nodes(x)
  und <- bnlearn::undirected.arcs(x)
  dir <- bnlearn::directed.arcs(x)
  if (nrow(und)) {
    keep <- und[, 1] < und[, 2]
    return(pdag(if (nrow(dir)) dir else NULL, und[keep, , drop = FALSE], nodes = nm))
  }
  dag(if (nrow(arcs)) arcs else NULL, nodes = nm)
}

#' Convert to and from bnlearn
#'
#' `from_bnlearn()` converts a `bn` (structure) or `bn.fit` (fitted network
#' with `bn.fit.dnode`, `bn.fit.onode`, `bn.fit.gnode` or `bn.fit.cgnode`
#' nodes) into pgmcompR. `as_bnlearn()` converts a pgmcompR graph to `bn`
#' and a fitted discrete or linear Gaussian network to `bn.fit`. Ordinal
#' nodes become ordinary discrete nodes (the ordering flag is reported in a
#' message). Conditional Gaussian networks are converted from bnlearn but
#' not back. Both directions preserve state order and parameter values.
#'
#' @param x A 'bnlearn' object, or a pgmcompR graph or network.
#' @return A pgmcompR object or a 'bnlearn' object.
#' @export
from_bnlearn <- function(x) {
  if (!requireNamespace("bnlearn", quietly = TRUE)) pgm_abort("The 'bnlearn' package is required.")
  if (inherits(x, "bn")) return(bnlearn_graph(x))
  if (!inherits(x, "bn.fit")) pgm_abort("`x` must be a bnlearn `bn` or `bn.fit` object.")
  nm <- names(x)
  edges <- list()
  for (v in nm) for (p in x[[v]]$parents) edges[[length(edges) + 1L]] <- c(p, v)
  bn <- bayesian_network(edges, nodes = nm)
  ordinal <- character()
  for (v in bnlearn::node.ordering(x)) {
    node <- x[[v]]
    if (inherits(node, "bn.fit.onode")) ordinal <- c(ordinal, v)
    if (inherits(node, c("bn.fit.dnode", "bn.fit.onode"))) {
      tab <- node$prob
      pa <- node$parents
      dn <- dimnames(tab)
      st <- dn
      if (is.null(names(st))) names(st) <- c(v, pa)
      bn <- add_cpd(bn, cpd(v, as.numeric(tab), parents = pa, states = st))
    } else if (inherits(node, "bn.fit.gnode")) {
      co <- node$coefficients
      bn <- add_cpd(bn, lg_cpd(v, intercept = co[[1L]], coefficients = if (length(co) > 1L) co[-1L] else NULL, sd = node$sd))
    } else if (inherits(node, "bn.fit.cgnode")) {
      dp <- node$parents[node$dparents]
      cp <- node$parents[node$gparents]
      co <- node$coefficients
      dst <- lapply(setNames(dp, dp), function(p) bn$registry[[p]]$states %||% levels(node$dlevels[[p]]) %||% node$dlevels[[p]])
      coef <- t(co[-1L, , drop = FALSE])
      colnames(coef) <- cp
      bn <- add_cpd(bn, cg_cpd(v, dp, dst, intercept = co[1L, ], coefficients = coef, sd = node$sd, continuous_parents = cp))
    } else pgm_abort(sprintf("Unsupported bnlearn node class for '%s'.", v))
  }
  if (length(ordinal)) message("Ordinal nodes were converted to unordered discrete nodes: ", paste(ordinal, collapse = ", "), ".")
  bn
}

#' @rdname from_bnlearn
#' @export
as_bnlearn <- function(x) {
  if (!requireNamespace("bnlearn", quietly = TRUE)) pgm_abort("The 'bnlearn' package is required.")
  g <- as_pgm_graph(x)
  if (any(g$edges$type != "directed")) pgm_abort("Only DAGs can be converted to bnlearn structures here.")
  s <- bnlearn::empty.graph(g$nodes)
  if (nrow(g$edges)) bnlearn::arcs(s) <- as.matrix(g$edges[, c("from", "to")])
  if (!inherits(x, "bayesian_network") || !length(x$cpds)) return(s)
  check_model(x)
  if (is_discrete_bn(x)) {
    dist <- lapply(x$cpds, function(c) {
      arr <- array(c$values, dim = unname(c$card), dimnames = c$states)
      if (length(c$parents)) arr else as.table(setNames(c$values, c$states[[1L]]))
    })
    return(bnlearn::custom.fit(s, dist = dist))
  }
  if (is_gaussian_bn(x)) {
    dist <- lapply(x$cpds, function(c) list(coef = c(`(Intercept)` = c$intercept, c$coefficients), sd = c$sd))
    return(bnlearn::custom.fit(s, dist = dist))
  }
  pgm_abort("Conditional Gaussian networks cannot be exported to bnlearn here.", class = "pgm_unsupported_feature")
}

#' Convert to and from gRain
#'
#' `as_grain()` compiles a discrete pgmcompR Bayesian network into a 'gRain'
#' object; `from_grain()` converts a 'gRain' network (its CPTs) back.
#'
#' @param x A discrete `bayesian_network` or a `grain` object.
#' @return A `grain` object or a `bayesian_network`.
#' @export
as_grain <- function(x) {
  if (!requireNamespace("gRain", quietly = TRUE)) pgm_abort("The 'gRain' package is required.")
  without_grbase_attach(as_grain_impl(x))
}

as_grain_impl <- function(x) {
  require_discrete_bn(x, "gRain")
  tabs <- lapply(x$cpds, function(c) {
    vp <- c(c$child, c$parents)
    gRain::cptable(vp, values = c$values, levels = c$states[[1L]])
  })
  # cptable needs parent levels from other tables; compile together.
  gRain::grain(gRain::compileCPT(tabs))
}

#' @rdname as_grain
#' @export
from_grain <- function(x) {
  if (!requireNamespace("gRain", quietly = TRUE)) pgm_abort("The 'gRain' package is required.")
  without_grbase_attach(from_grain_impl(x))
}

from_grain_impl <- function(x) {
  cl <- x$cptlist %||% gRain::getgrain(x, "cpt")
  edges <- list(); cps <- list(); states <- list()
  for (v in names(cl)) {
    tab <- cl[[v]]
    dn <- dimnames(tab)
    vars <- names(dn)
    pa <- vars[-1L]
    for (k in seq_along(vars)) states[[vars[k]]] <- dn[[k]]
    cps[[v]] <- cpd(v, as.numeric(tab), parents = pa, states = dn)
    for (p in pa) edges[[length(edges) + 1L]] <- c(p, v)
  }
  bn <- bayesian_network(edges, nodes = names(cl), states = states)
  for (c in cps) bn <- add_cpd(bn, c)
  bn
}

# gRain attaches gRbase as a side effect, which masks several pgmcompR
# functions (nodes(), edges(), parents(), ...). Undo that attachment when
# our call caused it.
without_grbase_attach <- function(code) {
  was_attached <- "package:gRbase" %in% search()
  on.exit(if (!was_attached && "package:gRbase" %in% search()) {
    detach("package:gRbase", character.only = TRUE)
  })
  code
}
