# Belief propagation on factor graphs ----------------------------------------------

fg_is_tree <- function(factors, vars) {
  # The bipartite graph is a forest iff #edges = #nodes - #components.
  n_edges <- sum(vapply(factors, function(f) length(f$vars), integer(1)))
  nodes_all <- c(vars, paste0("\rf", seq_along(factors)))
  parent <- setNames(seq_along(nodes_all), nodes_all)
  find <- function(x) {
    while (parent[[x]] != match(x, nodes_all)) x <- nodes_all[parent[[x]]]
    x
  }
  comps <- length(nodes_all)
  for (i in seq_along(factors)) {
    fn <- paste0("\rf", i)
    for (v in factors[[i]]$vars) {
      a <- find(fn); b <- find(v)
      if (a != b) {
        parent[[a]] <- match(b, nodes_all)
        comps <- comps - 1L
      }
    }
  }
  n_edges == length(nodes_all) - comps
}

#' Belief propagation
#'
#' Runs sum product (marginals) or max product (max marginals) message
#' passing on the factor graph of a discrete model with a flooding schedule.
#' On a tree structured factor graph the messages converge after a number of
#' iterations bounded by the diameter and the beliefs are exact. On graphs
#' with cycles the result is loopy belief propagation, an approximation that
#' may fail to converge; damping (a convex combination of the old and new
#' factor to variable messages) often helps. Non convergence within
#' `max_iter` iterations raises a warning of class `pgm_convergence_warning`.
#'
#' `loopy_bp()` is `belief_propagation()` with `allow_loops = TRUE` and
#' default damping 0.5.
#'
#' @param model A discrete model or factor graph.
#' @param evidence Named list of observed states.
#' @param type `"sum"` for marginals or `"max"` for max marginals.
#' @param max_iter Maximum number of iterations.
#' @param tol Convergence tolerance on the largest absolute change of a
#'   normalized message.
#' @param damping Damping weight in `[0, 1)` on the previous message.
#' @param allow_loops If `FALSE` (default for `belief_propagation()`) a
#'   factor graph with cycles is an error, so exact and approximate results
#'   are never confused.
#' @param virtual_evidence Named list of likelihood vectors.
#' @return A list of class `pgm_bp` with `beliefs` (named list of normalized
#'   single variable factors), `converged`, `iterations`, `max_residual`,
#'   `exact` (TRUE on trees) and, for `type = "max"`, `assignment`.
#' @examples
#' bn <- example_network("asia")
#' bp <- belief_propagation(as_factor_graph(bn), evidence = list(xray = "yes"), allow_loops = TRUE)
#' bp$exact
#' st <- list(A = c("0", "1"), B = c("0", "1"), C = c("0", "1"))
#' pot <- function(a, b) discrete_factor(c(a, b), st[c(a, b)], c(3, 1, 1, 3))
#' mn <- markov_network(list(pot("A", "B"), pot("B", "C"), pot("C", "A")))
#' loopy_bp(mn)$beliefs$A
#' @export
belief_propagation <- function(model, evidence = list(), type = c("sum", "max"), max_iter = 100L,
                               tol = 1e-10, damping = 0, allow_loops = FALSE, virtual_evidence = NULL) {
  type <- match.arg(type)
  if (damping < 0 || damping >= 1) pgm_abort("`damping` must lie in [0, 1).")
  prep <- prepare_discrete(model, character(), evidence, virtual_evidence, prune = FALSE)
  fs <- reduce_all(prep$factors, prep$evidence)
  const <- vapply(fs, function(f) length(f$vars) == 0L, logical(1))
  if (any(const) && any(vapply(fs[const], function(f) f$values == 0, logical(1)))) impossible_evidence(prep$evidence)
  fs <- fs[!const]
  vars <- setdiff(nodes(model), names(prep$evidence))
  vars <- intersect(vars, unique(unlist(lapply(fs, `[[`, "vars"), use.names = FALSE)))
  tree <- fg_is_tree(fs, vars)
  if (!tree && !allow_loops) {
    pgm_abort("The factor graph has cycles, so belief propagation is not exact. Use `loopy_bp()` (approximate) or an exact method such as variable elimination or the junction tree.",
              class = "pgm_loopy_graph")
  }
  st <- prep$states
  unif <- function(v) new_factor(v, st[v], rep(1 / length(st[[v]]), length(st[[v]])))
  nf <- length(fs)
  fvars <- lapply(fs, `[[`, "vars")
  v2f <- list(); f2v <- list()
  for (i in seq_len(nf)) for (v in fvars[[i]]) {
    v2f[[paste(v, i)]] <- unif(v)
    f2v[[paste(i, v)]] <- unif(v)
  }
  var_fac <- lapply(setNames(vars, vars), function(v) which(vapply(fvars, function(s) v %in% s, logical(1))))
  converged <- FALSE
  resid <- Inf
  it <- 0L
  normalize_msg <- function(m) {
    s <- sum(m$values)
    if (!(s > 0)) impossible_evidence(prep$evidence)
    m$values <- m$values / s
    m
  }
  for (it in seq_len(max_iter)) {
    resid <- 0
    new_f2v <- f2v
    for (i in seq_len(nf)) {
      for (v in fvars[[i]]) {
        others <- setdiff(fvars[[i]], v)
        p <- factor_product(c(list(fs[[i]]), lapply(others, function(u) v2f[[paste(u, i)]])))
        m <- if (type == "sum") factor_marginalize(p, others) else factor_maximize(p, others)
        m <- normalize_msg(m)
        old <- f2v[[paste(i, v)]]
        if (damping > 0) m$values <- (1 - damping) * m$values + damping * old$values
        resid <- max(resid, max(abs(m$values - old$values)))
        new_f2v[[paste(i, v)]] <- m
      }
    }
    f2v <- new_f2v
    for (v in vars) {
      fi <- var_fac[[v]]
      for (i in fi) {
        inc <- lapply(setdiff(fi, i), function(j) f2v[[paste(j, v)]])
        m <- if (length(inc)) normalize_msg(factor_product(c(list(unif(v)), inc))) else unif(v)
        v2f[[paste(v, i)]] <- m
      }
    }
    if (resid < tol) {
      converged <- TRUE
      break
    }
  }
  beliefs <- lapply(setNames(vars, vars), function(v) {
    inc <- lapply(var_fac[[v]], function(j) f2v[[paste(j, v)]])
    factor_normalize(factor_product(c(list(unif(v)), inc)))
  })
  # Variables with no factor after reduction and no evidence are uniform.
  for (v in setdiff(setdiff(nodes(model), names(prep$evidence)), vars)) beliefs[[v]] <- unif(v)
  beliefs <- beliefs[intersect(nodes(model), names(beliefs))]
  if (!converged) {
    pgm_warn(sprintf("Belief propagation did not converge in %d iterations (max residual %.3g).", max_iter, resid),
             class = "pgm_convergence_warning")
  }
  out <- list(beliefs = beliefs, converged = converged, iterations = it, max_residual = resid,
              exact = tree && converged, type = type, evidence = prep$evidence)
  if (type == "max") {
    out$assignment <- vapply(beliefs, function(b) b$states[[1L]][which.max(b$values)], "")
  }
  class(out) <- "pgm_bp"
  out
}

#' @rdname belief_propagation
#' @param ... Passed to `belief_propagation()`.
#' @export
loopy_bp <- function(model, evidence = list(), damping = 0.5, max_iter = 200L, ...) {
  belief_propagation(model, evidence = evidence, damping = damping, max_iter = max_iter, allow_loops = TRUE, ...)
}

#' @export
print.pgm_bp <- function(x, ...) {
  cat(sprintf("<Belief propagation (%s product): %s, %d iterations, max residual %.2g, %s>\n",
              x$type, if (x$converged) "converged" else "NOT converged", x$iterations, x$max_residual,
              if (x$exact) "exact (tree)" else "approximate"))
  for (v in names(x$beliefs)) {
    b <- x$beliefs[[v]]
    cat(sprintf("  %s: %s\n", v, paste(sprintf("%s=%.4f", b$states[[1L]], b$values), collapse = "  ")))
  }
  invisible(x)
}

infer_bp <- function(model, variables, evidence, joint = TRUE, virtual_evidence = NULL,
                     interventions = NULL, elimination_order = NULL, log = FALSE, n = NULL, seed = NULL, ...) {
  interventions <- normalize_evidence(interventions, "interventions")
  if (length(interventions)) model <- intervene(model, interventions)
  if (length(variables) > 1L && joint) {
    pgm_abort("Belief propagation returns single variable marginals; use `joint = FALSE` or method = 've' for joint queries.")
  }
  res <- belief_propagation(model, evidence, virtual_evidence = virtual_evidence, ...)
  bp_result(res, variables, joint, log, "bp")
}

infer_loopy_bp <- function(model, variables, evidence, joint = TRUE, virtual_evidence = NULL,
                     interventions = NULL, elimination_order = NULL, log = FALSE, n = NULL, seed = NULL, ...) {
  interventions <- normalize_evidence(interventions, "interventions")
  if (length(interventions)) model <- intervene(model, interventions)
  if (length(variables) > 1L && joint) {
    pgm_abort("Loopy belief propagation returns single variable marginals; use `joint = FALSE`.")
  }
  res <- loopy_bp(model, evidence, virtual_evidence = virtual_evidence, ...)
  out <- bp_result(res, variables, joint, log, "loopy_bp")
  mark <- function(f) {
    if (!res$exact) class(f) <- c("pgm_approximate", class(f))
    attr(f, "diagnostics") <- list(converged = res$converged, iterations = res$iterations,
                                   max_residual = res$max_residual, exact = res$exact)
    f
  }
  if (inherits(out, "pgm_factor")) mark(out) else lapply(out, mark)
}

bp_result <- function(res, variables, joint, log, method) {
  ev <- intersect(variables, names(res$evidence))
  if (length(ev)) pgm_abort(paste0("Query variables cannot also be observed: ", paste(ev, collapse = ", "), "."))
  out <- lapply(setNames(variables, variables), function(v) {
    b <- res$beliefs[[v]]
    if (log) b <- factor_to_log(b)
    attr(b, "method") <- method
    b
  })
  if (length(variables) == 1L) out[[1L]] else out
}
