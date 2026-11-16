# Causal graphical models -------------------------------------------------------------

#' Interventions by graph mutilation
#'
#' `intervene()` returns the mutilated network for do(X = x): incoming arcs
#' of each intervened variable are removed and its CPD is replaced by a point
#' mass on the set value (truncated factorization). For a continuous
#' variable the CPD becomes a root linear Gaussian CPD; queries then
#' condition on the set value, which is equivalent in the mutilated graph.
#'
#' `do_query()` computes P(variables | do(interventions), evidence) by
#' inference in the mutilated network. It differs from conditioning,
#' P(variables | X = x), whenever X has a back door path to the query
#' variables. These computations assume the network is a causal Bayesian
#' network: its arcs are causal and there are no unmodeled confounders.
#'
#' @param model A `bayesian_network`.
#' @param interventions Named list of states (or numbers for continuous
#'   variables).
#' @param variables Query variables.
#' @param evidence Named list of observations made in the intervened world.
#' @param method Inference method for discrete networks.
#' @param ... Passed to [query()].
#' @return `intervene()` returns a `bayesian_network`; `do_query()` returns
#'   the same type as [query()].
#' @examples
#' bn <- example_network("sprinkler")
#' # Observing the sprinkler on lowers belief in rain
#' query(bn, "Rain", evidence = list(Sprinkler = "T"))
#' # Turning it on does not change rain
#' do_query(bn, "Rain", interventions = list(Sprinkler = "T"))
#' @export
intervene <- function(model, interventions) {
  if (!inherits(model, "bayesian_network")) {
    pgm_abort("Interventions are defined for (causal) Bayesian networks only.")
  }
  if (!length(interventions)) return(model)
  if (is.null(names(interventions)) || any(!nzchar(names(interventions)))) pgm_abort("`interventions` must be named.")
  check_nodes_in(model$graph, names(interventions), "interventions")
  g <- model$graph
  for (v in names(interventions)) {
    pa <- parents(g, v)
    if (length(pa)) g <- remove_edges(g, pa, rep(v, length(pa)))
  }
  model$graph <- g
  for (v in names(interventions)) {
    val <- interventions[[v]]
    r <- model$registry[[v]]
    if (identical(r$type, "continuous")) {
      if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) pgm_abort(sprintf("Intervention on '%s' needs a number.", v))
      model$cpds[[v]] <- lg_cpd(v, intercept = val, sd = 1)
    } else {
      val <- as.character(val)
      if (!val %in% r$states) pgm_abort(sprintf("'%s' is not a state of '%s'.", val, v), class = "pgm_unknown_state")
      model$cpds[[v]] <- new_cpd(v, character(), list(r$states), as.numeric(r$states == val))
    }
  }
  model$interventions <- c(model$interventions, interventions)
  model
}

#' @rdname intervene
#' @export
do_query <- function(model, variables, interventions, evidence = list(), method = "ve", ...) {
  if (!length(interventions)) pgm_abort("`interventions` must not be empty.")
  query(model, variables, evidence = evidence, interventions = interventions, method = method, ...)
}

# Adjustment criteria -----------------------------------------------------------------

remove_outgoing <- function(g, x) {
  e <- g$edges
  keep <- !(e$from %in% x & e$type == "directed")
  g$edges <- e[keep, , drop = FALSE]
  g
}

#' Back door and front door criteria
#'
#' `is_backdoor_set()` tests Pearl's back door criterion: no element of `z`
#' is a descendant of `x`, and `z` blocks every path between `x` and `y` that
#' starts with an arrow into `x` (d-separation in the graph with the arcs out
#' of `x` removed). `is_frontdoor_set()` tests the front door criterion:
#' `z` intercepts every directed path from `x` to `y`, there is no unblocked
#' back door path from `x` to `z`, and every back door path from `z` to `y`
#' is blocked by `x`. `adjustment_sets()` enumerates back door sets drawn
#' from observed non descendants of `x`; with `minimal = TRUE` only sets with
#' no valid proper subset are returned.
#'
#' @param g A DAG or Bayesian network (may include latent nodes).
#' @param x,y Treatment and outcome nodes.
#' @param z Candidate adjustment set.
#' @param latent Nodes that are unobserved and cannot be adjusted for.
#' @param minimal Return only minimal sets.
#' @param max_size Largest set size considered.
#' @return Logical scalar, or a list of character vectors.
#' @examples
#' g <- dag(list(c("Z", "X"), c("Z", "Y"), c("X", "Y")))
#' is_backdoor_set(g, "X", "Y", "Z")
#' adjustment_sets(g, "X", "Y")
#' @export
is_backdoor_set <- function(g, x, y, z = character()) {
  g <- as_pgm_graph(g)
  check_nodes_in(g, c(x, y, z))
  if (length(intersect(z, c(x, y)))) return(FALSE)
  if (any(z %in% descendants(g, x))) return(FALSE)
  d_separated(remove_outgoing(g, x), x, y, z)
}

#' @rdname is_backdoor_set
#' @export
is_frontdoor_set <- function(g, x, y, z) {
  g <- as_pgm_graph(g)
  check_nodes_in(g, c(x, y, z))
  if (!length(z) || length(intersect(z, c(x, y)))) return(FALSE)
  # (1) z intercepts all directed paths from x to y
  gz <- remove_outgoing(g, z)
  if (y %in% descendants(gz, x)) return(FALSE)
  # (2) no unblocked back door path from x to z
  if (!d_separated(remove_outgoing(g, x), x, z, character())) return(FALSE)
  # (3) all back door paths from z to y are blocked by x
  d_separated(remove_outgoing(g, z), z, y, x)
}

#' @rdname is_backdoor_set
#' @export
adjustment_sets <- function(g, x, y, latent = character(), minimal = TRUE, max_size = 8L) {
  g <- as_pgm_graph(g)
  check_nodes_in(g, c(x, y, latent))
  cand <- setdiff(g$nodes, c(x, y, latent, descendants(g, x)))
  out <- list()
  for (k in 0:min(max_size, length(cand))) {
    sets <- if (k == 0L) list(character()) else utils::combn(cand, k, simplify = FALSE)
    for (s in sets) {
      if (minimal && any(vapply(out, function(o) all(o %in% s), logical(1)))) next
      if (is_backdoor_set(g, x, y, s)) out[[length(out) + 1L]] <- s
    }
  }
  out
}

#' Identify a causal effect
#'
#' Looks for an identification strategy for P(y | do(x)) in a causal DAG with
#' possibly latent nodes, trying back door adjustment (smallest minimal set
#' first) and then the front door criterion with observed mediators. The
#' result gives the estimand as a formula. When neither applies, an error of
#' class `pgm_not_identifiable` reports that the effect is not identified by
#' these criteria (the general ID algorithm is not implemented, so other
#' identifiable effects are reported as unresolved rather than claimed
#' non identifiable).
#'
#' @param g A DAG or Bayesian network.
#' @param x,y Single treatment and outcome nodes.
#' @param latent Unobserved nodes.
#' @return A list of class `pgm_identification` with `strategy`, `set` and
#'   `estimand`.
#' @examples
#' g <- dag(list(c("U", "X"), c("U", "Y"), c("X", "M"), c("M", "Y")))
#' identify_effect(g, "X", "Y", latent = "U")
#' @export
identify_effect <- function(g, x, y, latent = character()) {
  g <- as_pgm_graph(g)
  check_string(x); check_string(y)
  sets <- adjustment_sets(g, x, y, latent)
  if (length(sets)) {
    s <- sets[[which.min(vapply(sets, length, integer(1)))]]
    est <- if (length(s)) sprintf("sum_{%s} P(%s | %s, %s) P(%s)", paste(s, collapse = ","), y, x, paste(s, collapse = ","), paste(s, collapse = ","))
    else sprintf("P(%s | %s)", y, x)
    return(structure(list(strategy = "backdoor", set = s, estimand = est, x = x, y = y), class = "pgm_identification"))
  }
  cand <- setdiff(g$nodes, c(x, y, latent))
  for (k in seq_len(min(4L, length(cand)))) {
    for (s in utils::combn(cand, k, simplify = FALSE)) {
      if (is_frontdoor_set(g, x, y, s)) {
        zs <- paste(s, collapse = ",")
        est <- sprintf("sum_{%s} P(%s | %s) sum_{%s'} P(%s | %s', %s) P(%s')", zs, zs, x, x, y, x, zs, x)
        return(structure(list(strategy = "frontdoor", set = s, estimand = est, x = x, y = y), class = "pgm_identification"))
      }
    }
  }
  pgm_abort(sprintf("P(%s | do(%s)) is not identified by the back door or front door criterion given the latent nodes {%s}.",
                    y, x, paste(latent, collapse = ", ")), class = "pgm_not_identifiable")
}

#' @export
print.pgm_identification <- function(x, ...) {
  cat(sprintf("<Identified: P(%s | do(%s)) by %s criterion, set {%s}>\n  %s\n", x$y, x$x, x$strategy,
              paste(x$set, collapse = ", "), x$estimand))
  invisible(x)
}

# Structural causal models ---------------------------------------------------------------

#' Discrete structural causal models
#'
#' An SCM specifies exogenous variables with independent prior
#' distributions and, for every endogenous variable, a deterministic
#' structural function of its parents (endogenous and exogenous). This
#' explicit specification is what licenses counterfactual reasoning; an
#' ordinary Bayesian network does not determine counterfactuals.
#'
#' `counterfactual()` evaluates P(variables under do(interventions) | evidence)
#' by abduction (the posterior of the exogenous variables given the factual
#' evidence), action (mutilation) and prediction (Pearl, 2009, Chapter 7).
#' The computation enumerates exogenous configurations with positive
#' posterior probability.
#'
#' @param exogenous Named list of probability vectors (named by state) for
#'   the exogenous variables.
#' @param equations Named list of functions; each takes the values of its
#'   parents as named character arguments and returns one state.
#' @param states Named list of states of the endogenous variables.
#' @param parents Named list of parent names for each endogenous variable.
#' @param model A `pgm_scm`.
#' @param variables Query variables.
#' @param evidence Factual evidence.
#' @param interventions Counterfactual interventions.
#' @return `scm()` returns an object of class `pgm_scm` (with element `bn`, the
#'   induced Bayesian network); `counterfactual()` returns a `pgm_factor`.
#' @examples
#' m <- scm(
#'   exogenous = list(U = c(a = 0.6, b = 0.4), V = c(a = 0.7, b = 0.3)),
#'   parents = list(X = "U", Y = c("X", "V")),
#'   states = list(X = c("0", "1"), Y = c("0", "1")),
#'   equations = list(
#'     X = function(U) if (U == "b") "1" else "0",
#'     Y = function(X, V) if (X == "1" || V == "b") "1" else "0"
#'   )
#' )
#' # Had X been 0, would Y have been 1, given that X = 1 and Y = 1?
#' counterfactual(m, "Y", evidence = list(X = "1", Y = "1"), interventions = list(X = "0"))
#' @export
scm <- function(exogenous, equations, states, parents) {
  if (!is.list(exogenous) || is.null(names(exogenous))) pgm_abort("`exogenous` must be a named list.")
  if (!setequal(names(equations), names(parents)) || !all(names(equations) %in% names(states))) {
    pgm_abort("`equations`, `parents` and `states` must name the same endogenous variables.")
  }
  endo <- names(equations)
  edges <- list()
  for (v in endo) for (p in parents[[v]]) edges[[length(edges) + 1L]] <- c(p, v)
  bn <- bayesian_network(edges, nodes = c(names(exogenous), endo))
  for (u in names(exogenous)) bn <- add_cpd(bn, cpd(u, exogenous[[u]]))
  allst <- c(lapply(exogenous, names), states)
  for (v in topological_sort(bn$graph)) {
    if (v %in% names(exogenous)) next
    pa <- parents(bn$graph, v)
    fun <- equations[[v]]
    bn <- add_cpd(bn, deterministic_cpd(v, pa, allst[c(v, pa)], function(p) do.call(fun, as.list(p))))
  }
  structure(list(bn = bn, exogenous = names(exogenous), endogenous = endo, equations = equations),
            class = c("pgm_scm", "pgm_model"))
}

#' @export
print.pgm_scm <- function(x, ...) {
  cat(sprintf("<Structural causal model: %d exogenous, %d endogenous variables>\n",
              length(x$exogenous), length(x$endogenous)))
  for (v in x$endogenous) cat(sprintf("  %s := f(%s)\n", v, paste(parents(x$bn$graph, v), collapse = ", ")))
  invisible(x)
}

#' @rdname scm
#' @export
counterfactual <- function(model, variables, evidence = list(), interventions = list()) {
  if (!inherits(model, "pgm_scm")) pgm_abort("Counterfactuals need an explicit structural causal model (`scm()`).")
  bn <- model$bn
  evidence <- normalize_evidence(evidence)
  u <- model$exogenous
  post <- query(bn, u, evidence = evidence, method = "ve")
  grid <- expand.grid(post$states, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  mut <- intervene(bn, normalize_evidence(interventions, "interventions"))
  out_states <- model_states(bn)[variables]
  acc <- setNames(numeric(prod(vapply(out_states, length, integer(1)))), NULL)
  for (r in which(post$values > 0)) {
    ev_u <- as.list(grid[r, , drop = FALSE])
    q <- query(mut, variables, evidence = ev_u, method = "ve")
    acc <- acc + post$values[r] * factor_reorder(q, variables)$values
  }
  f <- new_factor(variables, out_states, acc)
  attr(f, "method") <- "counterfactual"
  f
}
