# Exact inference ---------------------------------------------------------------

# Interaction graph adjacency for a set of factors, restricted to `vars`.
interaction_adjacency <- function(factors, vars) {
  scopes <- lapply(factors, function(f) intersect(f$vars, vars))
  owner <- rep(seq_along(scopes), lengths(scopes))
  flat <- unlist(scopes, use.names = FALSE)
  adj <- lapply(setNames(vars, vars), function(v) {
    u <- unique(flat[owner %in% owner[flat == v]])
    u[u != v]
  })
  adj
}

#' Variable elimination orderings
#'
#' Computes a greedy elimination ordering on the interaction graph of a set
#' of factors (two variables interact when they share a factor).
#' `"min_fill"` eliminates the variable whose elimination adds the fewest
#' fill-in edges; `"weighted_min_fill"` weights each fill edge by the product
#' of the cardinalities of its endpoints; `"min_degree"` picks the variable
#' with fewest neighbors; `"min_weight"` minimizes the size of the created
#' factor. Ties are broken by the order of `variables`.
#'
#' @param model A model or a list of factors.
#' @param variables Variables to eliminate (default: all).
#' @param heuristic One of `"min_fill"`, `"weighted_min_fill"`,
#'   `"min_degree"`, `"min_weight"`.
#' @return A character vector.
#' @examples
#' bn <- example_network("asia")
#' elimination_order(bn, heuristic = "min_degree")
#' @export
elimination_order <- function(model, variables = NULL, heuristic = c("min_fill", "weighted_min_fill", "min_degree", "min_weight")) {
  heuristic <- match.arg(heuristic)
  factors <- model_factors(model)
  all_vars <- unique(unlist(lapply(factors, `[[`, "vars"), use.names = FALSE))
  variables <- variables %||% all_vars
  card <- unlist(lapply(factors, function(f) f$card), use.names = TRUE)
  card <- card[!duplicated(names(card))]
  adj <- interaction_adjacency(factors, all_vars)
  greedy_order(adj, heuristic, weights = card, eliminate = intersect(all_vars, variables))
}

# Eliminate `eliminate` from `factors` in `order` by summation (or max).
eliminate_vars <- function(factors, order, op = c("sum", "max"), trace = FALSE) {
  op <- match.arg(op)
  pool <- factors
  traces <- list()
  for (v in order) {
    hit <- vapply(pool, function(f) v %in% f$vars, logical(1))
    if (!any(hit)) next
    prod_f <- factor_product(pool[hit])
    if (op == "sum") {
      new_f <- factor_marginalize(prod_f, v)
    } else {
      new_f <- factor_maximize(prod_f, v, argmax = FALSE)
      if (trace) traces[[v]] <- prod_f
    }
    pool <- c(pool[!hit], list(new_f))
  }
  list(factors = pool, traces = traces)
}

resolve_order <- function(factors, eliminate, elimination_order) {
  if (!length(eliminate)) return(character())
  if (is.character(elimination_order) && length(elimination_order) == 1L &&
      elimination_order %in% c("min_fill", "weighted_min_fill", "min_degree", "min_weight")) {
    all_vars <- unique(unlist(lapply(factors, `[[`, "vars"), use.names = FALSE))
    card <- unlist(lapply(factors, function(f) f$card), use.names = TRUE)
    card <- card[!duplicated(names(card))]
    adj <- interaction_adjacency(factors, all_vars)
    return(greedy_order(adj, elimination_order, weights = card, eliminate = intersect(all_vars, eliminate)))
  }
  ord <- as.character(elimination_order)
  if (anyDuplicated(ord)) pgm_abort("`elimination_order` has duplicated variables.")
  ord <- ord[ord %in% eliminate]
  if (!setequal(ord, eliminate)) {
    pgm_abort("A user supplied `elimination_order` must list each eliminated variable exactly once.")
  }
  ord
}

# Inference problem setup --------------------------------------------------------

# Returns list(factors, states, variables, evidence, model) with evidence
# validated, interventions applied, virtual evidence turned into factors,
# and (for Bayesian networks) barren nodes pruned.
prepare_discrete <- function(model, variables, evidence, virtual_evidence = NULL,
                             interventions = NULL, prune = TRUE, log = FALSE) {
  evidence <- normalize_evidence(evidence)
  interventions <- normalize_evidence(interventions, "interventions")
  if (length(interventions)) model <- intervene(model, interventions)
  states <- model_states(model)
  all_vars <- nodes(model)
  variables <- as.character(variables %||% character())
  check_names(variables, "variables")
  bad <- setdiff(c(variables, names(evidence), names(virtual_evidence)), all_vars)
  if (length(bad)) pgm_abort(paste0("Unknown variable(s): ", paste(bad, collapse = ", "), "."), class = "pgm_unknown_variable")
  for (v in names(evidence)) {
    if (!(as.character(evidence[[v]]) %in% states[[v]])) {
      pgm_abort(sprintf("'%s' is not a state of '%s' (states: %s).", evidence[[v]], v, paste(states[[v]], collapse = ", ")),
                class = "pgm_unknown_state")
    }
    evidence[[v]] <- as.character(evidence[[v]])
  }
  overlap <- intersect(variables, names(evidence))
  if (length(overlap)) {
    pgm_abort(paste0("Query variables cannot also be observed: ", paste(overlap, collapse = ", "), "."))
  }
  vfactors <- virtual_evidence_factors(virtual_evidence, states)
  if (inherits(model, "bayesian_network")) {
    check_model(model)
    if (prune) {
      keep_set <- unique(c(variables, names(evidence), names(virtual_evidence)))
      keep <- if (length(keep_set)) c(keep_set, ancestors(model$graph, keep_set)) else all_vars
      model$cpds <- model$cpds[intersect(all_vars, keep)]
    }
  }
  factors <- if (inherits(model, "bayesian_network")) unname(bn_factors(model)) else unname(model_factors(model))
  factors <- c(factors, vfactors)
  if (log) factors <- lapply(factors, factor_to_log)
  list(factors = factors, states = states, variables = variables, evidence = evidence, model = model)
}

virtual_evidence_factors <- function(virtual_evidence, states) {
  if (is.null(virtual_evidence) || !length(virtual_evidence)) return(list())
  if (!is.list(virtual_evidence) || is.null(names(virtual_evidence))) {
    pgm_abort("`virtual_evidence` must be a named list of likelihood vectors.")
  }
  lapply(names(virtual_evidence), function(v) {
    lik <- virtual_evidence[[v]]
    st <- states[[v]]
    if (is.null(st)) pgm_abort(sprintf("Virtual evidence on '%s' requires a discrete variable.", v))
    if (!is.null(names(lik))) {
      if (!setequal(names(lik), st)) pgm_abort(sprintf("Virtual evidence on '%s' must name every state.", v))
      lik <- lik[st]
    }
    if (length(lik) != length(st)) pgm_abort(sprintf("Virtual evidence on '%s' needs %d likelihoods.", v, length(st)))
    check_prob(lik, "virtual_evidence")
    if (all(lik == 0)) impossible_evidence()
    discrete_factor(v, setNames(list(st), v), as.numeric(lik))
  })
}

reduce_all <- function(factors, evidence) {
  lapply(factors, factor_reduce, evidence = evidence)
}

# Variable elimination ------------------------------------------------------------

ve_posterior <- function(prep, elimination_order = "min_fill") {
  fs <- reduce_all(prep$factors, prep$evidence)
  scope <- unique(unlist(lapply(fs, `[[`, "vars"), use.names = FALSE))
  missing <- setdiff(prep$variables, scope)
  # Query variables with no remaining factor (cannot happen for BNs) get uniform.
  if (length(missing)) pgm_abort(paste0("No factor mentions: ", paste(missing, collapse = ", "), "."))
  elim <- setdiff(scope, prep$variables)
  ord <- resolve_order(fs, elim, elimination_order)
  res <- eliminate_vars(fs, ord, "sum")
  joint <- factor_product(res$factors)
  joint <- factor_reorder(joint, prep$variables)
  joint
}

finish_posterior <- function(joint, prep, log, method, joint_out) {
  z <- if (joint$log) log_sum_exp(joint$values) else sum(joint$values)
  if ((joint$log && z == -Inf) || (!joint$log && !(z > 0))) impossible_evidence(prep$evidence)
  post <- factor_normalize(joint)
  if (!log) post <- factor_to_linear(post)
  attr(post, "method") <- method
  attr(post, "evidence") <- prep$evidence
  if (joint_out || length(prep$variables) == 1L) return(post)
  out <- lapply(setNames(prep$variables, prep$variables), function(v) {
    m <- factor_marginalize(post, setdiff(post$vars, v))
    attr(m, "method") <- method
    m
  })
  out
}

#' Probabilistic queries
#'
#' Computes the posterior distribution of `variables` given hard evidence,
#' soft (likelihood) evidence and interventions. Exact methods:
#'
#' * `"ve"`: variable elimination (default);
#' * `"enumeration"`: sums the full joint over all completions of the
#'   evidence (exponential; a reference method);
#' * `"jt"`: junction tree calibration;
#' * `"bp"`: belief propagation on a factor graph that is a tree (exact); on
#'   graphs with cycles it stops with an error pointing to `"loopy_bp"`.
#'
#' Approximate methods return estimates with diagnostics (see
#' [approximate_inference]): `"loopy_bp"`, `"rejection"`,
#' `"likelihood_weighting"`, `"importance"`, `"gibbs"`. Additional methods can
#' be registered with [register_inference_method()].
#'
#' Hard evidence `list(X = "x")` conditions on the observation X = x.
#' Interventions `list(X = "x")` set X by the do operator: the CPD of X is
#' replaced by a point mass and its incoming arcs are cut (graph mutilation)
#' before inference, so the result is P(variables | do(X = x), evidence).
#' Soft evidence `list(X = c(x1 = 0.8, x2 = 0.2))` multiplies the joint by a
#' likelihood vector over the states of X. Variables absent from all three
#' are marginalized (missing evidence). Evidence of probability zero raises an
#' error of class `pgm_impossible_evidence`.
#'
#' For Gaussian and conditional Gaussian Bayesian networks `query()` returns
#' exact Gaussian or mixture results (see [gaussian_query()]).
#'
#' @param model A Bayesian network, Markov network, factor graph or junction
#'   tree.
#' @param variables Query variables.
#' @param evidence Named list of observed states.
#' @param method Inference method name.
#' @param joint If `TRUE` return the joint posterior of `variables`;
#'   otherwise a list of marginal posteriors.
#' @param virtual_evidence Named list of likelihood vectors.
#' @param interventions Named list of states set by intervention.
#' @param elimination_order Heuristic name or explicit order for `"ve"`.
#' @param log If `TRUE` compute in log space and return a log scale factor.
#' @param ... Passed to the inference method (for example `n`, `seed`).
#' @return A normalized `pgm_factor` (attribute `method`), or a named list of
#'   them when `joint = FALSE`. Approximate methods add class
#'   `pgm_approximate` and a `diagnostics` attribute.
#' @examples
#' bn <- example_network("asia")
#' query(bn, "lung", evidence = list(xray = "yes", smoke = "yes"))
#' query(bn, "lung", evidence = list(xray = "yes"), method = "jt")
#' query(bn, "dysp", interventions = list(smoke = "no"))
#' @export
query <- function(model, variables, evidence = list(), method = "ve", joint = TRUE,
                  virtual_evidence = NULL, interventions = NULL, elimination_order = "min_fill",
                  log = FALSE, ...) {
  if (inherits(model, "bayesian_network") && !is_discrete_bn(model) && length(model$cpds)) {
    if (!method %in% c("ve", "exact", "enumeration", "jt")) {
      pgm_abort("Continuous and hybrid networks support exact inference only (method = 've').")
    }
    return(gaussian_query(model, variables, evidence, interventions = interventions))
  }
  fun <- get_inference_method(method)
  fun(model = model, variables = variables, evidence = evidence, joint = joint,
      virtual_evidence = virtual_evidence, interventions = interventions,
      elimination_order = elimination_order, log = log, ...)
}

infer_ve <- function(model, variables, evidence, joint = TRUE, virtual_evidence = NULL,
                     interventions = NULL, elimination_order = "min_fill", log = FALSE, ...) {
  prep <- prepare_discrete(model, variables, evidence, virtual_evidence, interventions, log = log)
  j <- ve_posterior(prep, elimination_order)
  finish_posterior(j, prep, log, "ve", joint)
}

# Enumeration: sum the joint over every complete assignment consistent with
# the evidence, evaluating each factor by direct table lookup.
infer_enumeration <- function(model, variables, evidence, joint = TRUE, virtual_evidence = NULL,
                              interventions = NULL, log = FALSE, max_cells = 1e7, ...) {
  prep <- prepare_discrete(model, variables, evidence, virtual_evidence, interventions, prune = FALSE)
  fs <- prep$factors
  vars <- unique(unlist(lapply(fs, `[[`, "vars"), use.names = FALSE))
  st <- lapply(setNames(vars, vars), function(v) {
    if (v %in% names(prep$evidence)) prep$evidence[[v]] else prep$states[[v]]
  })
  n <- prod(vapply(st, length, integer(1)))
  if (n > max_cells) pgm_abort(sprintf("Enumeration needs %g assignments; increase `max_cells` or use method = 've'.", n))
  grid <- expand.grid(st, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  logp <- numeric(nrow(grid))
  for (f in fs) {
    pos <- integer(nrow(grid))
    stride <- 1
    for (k in seq_along(f$vars)) {
      v <- f$vars[k]
      pos <- pos + (match(grid[[v]], f$states[[v]]) - 1L) * stride
      stride <- stride * f$card[[k]]
    }
    val <- f$values[pos + 1]
    logp <- logp + if (f$log) val else log(val)
  }
  key <- if (length(prep$variables)) do.call(paste, c(unname(as.list(grid[prep$variables])), sep = "\r")) else rep("", nrow(grid))
  out_grid <- expand.grid(prep$states[prep$variables], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out_key <- if (length(prep$variables)) do.call(paste, c(unname(as.list(out_grid)), sep = "\r")) else ""
  lv <- vapply(out_key, function(k) {
    sel <- logp[key == k]
    if (length(sel)) log_sum_exp(sel) else -Inf
  }, numeric(1))
  joint_f <- new_factor(prep$variables, prep$states[prep$variables], lv, TRUE)
  finish_posterior(joint_f, prep, log, "enumeration", joint)
}

#' Probability of evidence and joint probabilities
#'
#' `evidence_probability()` returns P(evidence) (the partition function of the
#' reduced model for Markov networks, divided by the unreduced partition
#' function). `joint_probability()` returns the probability of a full or
#' partial assignment. `marginal_probability()` returns the prior marginal
#' distribution of `variables`. All use exact variable elimination.
#'
#' @param model A discrete model.
#' @param evidence,assignment Named list of states.
#' @param variables Variables.
#' @param log Return the log probability.
#' @param ... Passed to [query()].
#' @return A numeric scalar, or a `pgm_factor` for `marginal_probability()`.
#' @examples
#' bn <- example_network("asia")
#' evidence_probability(bn, list(xray = "yes", dysp = "yes"))
#' marginal_probability(bn, "lung")
#' @export
evidence_probability <- function(model, evidence, log = FALSE) {
  prep <- prepare_discrete(model, character(), evidence, log = TRUE)
  fs <- reduce_all(prep$factors, prep$evidence)
  scope <- unique(unlist(lapply(fs, `[[`, "vars"), use.names = FALSE))
  ord <- resolve_order(fs, scope, "min_fill")
  res <- eliminate_vars(fs, ord, "sum")
  lz <- sum(vapply(res$factors, function(f) f$values, numeric(1)))
  if (!inherits(model, "bayesian_network")) {
    lz <- lz - partition_function(model, log = TRUE)
  }
  if (log) lz else exp(lz)
}

#' @rdname evidence_probability
#' @export
joint_probability <- function(model, assignment, log = FALSE) {
  evidence_probability(model, assignment, log = log)
}

#' @rdname evidence_probability
#' @export
marginal_probability <- function(model, variables, ...) {
  query(model, variables, evidence = list(), ...)
}

# MAP and MPE ---------------------------------------------------------------------

#' MAP and MPE queries
#'
#' `map_query()` returns the most probable joint assignment of `variables`
#' given the evidence, maximizing the marginal posterior
#' P(variables | evidence) with the remaining variables summed out.
#' `mpe_query()` returns the most probable explanation: the jointly most
#' probable assignment of all unobserved variables. These generally differ:
#' MAP of a subset is not the restriction of the MPE.
#'
#' Both use exact variable elimination: summed variables are eliminated
#' first, then the MAP variables are maximized out and the maximizing
#' assignment is recovered by traceback. Ties are broken toward the first
#' state.
#'
#' @param model A discrete model.
#' @param variables Variables to maximize over.
#' @param evidence Named list of observed states.
#' @param interventions Named list of interventions.
#' @param virtual_evidence Named list of likelihood vectors.
#' @return A list of class `pgm_map` with `assignment` (named character
#'   vector) and `probability` (posterior probability of the assignment given
#'   the evidence).
#' @examples
#' bn <- example_network("asia")
#' map_query(bn, c("lung", "tub"), evidence = list(xray = "yes"))
#' mpe_query(bn, evidence = list(xray = "yes", dysp = "yes"))
#' @export
map_query <- function(model, variables, evidence = list(), interventions = NULL, virtual_evidence = NULL) {
  prep <- prepare_discrete(model, variables, evidence, virtual_evidence, interventions, log = TRUE)
  fs <- reduce_all(prep$factors, prep$evidence)
  scope <- unique(unlist(lapply(fs, `[[`, "vars"), use.names = FALSE))
  sum_vars <- setdiff(scope, prep$variables)
  ord_sum <- resolve_order(fs, sum_vars, "min_fill")
  res <- eliminate_vars(fs, ord_sum, "sum")
  fs2 <- res$factors
  lz <- {
    tmp <- eliminate_vars(fs2, resolve_order(fs2, prep$variables, "min_fill"), "sum")$factors
    sum(vapply(tmp, function(f) f$values, numeric(1)))
  }
  if (lz == -Inf) impossible_evidence(prep$evidence)
  ord_max <- resolve_order(fs2, prep$variables, "min_fill")
  mx <- eliminate_vars(fs2, ord_max, "max", trace = TRUE)
  lmax <- sum(vapply(mx$factors, function(f) f$values, numeric(1)))
  assignment <- setNames(character(length(ord_max)), ord_max)
  for (v in rev(ord_max)) {
    tf <- mx$traces[[v]]
    if (is.null(tf)) {
      assignment[[v]] <- prep$states[[v]][1L]
      next
    }
    others <- setdiff(tf$vars, v)
    red <- if (length(others)) factor_reduce(tf, as.list(assignment[others])) else tf
    assignment[[v]] <- red$states[[v]][which.max(red$values)]
  }
  assignment <- assignment[prep$variables]
  structure(list(assignment = assignment, probability = exp(lmax - lz), log_probability = lmax - lz,
                 evidence = prep$evidence, type = "MAP"), class = "pgm_map")
}

#' @rdname map_query
#' @export
mpe_query <- function(model, evidence = list(), interventions = NULL, virtual_evidence = NULL) {
  evidence <- normalize_evidence(evidence)
  vars <- setdiff(nodes(model), names(evidence))
  out <- map_query(model, vars, evidence, interventions, virtual_evidence)
  out$type <- "MPE"
  out
}

#' @export
print.pgm_map <- function(x, ...) {
  cat(sprintf("<%s assignment> posterior probability %s\n", x$type, format(x$probability, digits = 6)))
  print(x$assignment, quote = FALSE)
  invisible(x)
}
