# Markov networks and factor graphs ----------------------------------------------

factor_registry <- function(factors) {
  reg <- list()
  for (f in factors) {
    for (v in f$vars) {
      if (is.null(reg[[v]])) reg[[v]] <- list(type = "discrete", states = f$states[[v]])
      else if (!identical(reg[[v]]$states, f$states[[v]])) {
        pgm_abort(sprintf("Factors disagree on the states of '%s'.", v), class = "pgm_state_mismatch")
      }
    }
  }
  reg
}

#' Markov networks
#'
#' A discrete Markov network (Markov random field) is a set of nonnegative
#' potentials over arbitrary scopes. The joint distribution is the
#' normalized product of the potentials; the graph connects every pair of
#' variables that share a potential. Potentials need not be normalized.
#'
#' @param factors A list of `pgm_factor` potentials.
#' @param nodes Optional variable names (variables without potentials need
#'   `states`).
#' @param states Optional named list of states for variables without a
#'   potential.
#' @return An object of class `markov_network`.
#' @examples
#' st <- list(A = c("0", "1"), B = c("0", "1"), C = c("0", "1"))
#' ising <- function(a, b) discrete_factor(c(a, b), st[c(a, b)], c(2, 1, 1, 2))
#' mn <- markov_network(list(ising("A", "B"), ising("B", "C"), ising("C", "A")))
#' partition_function(mn)
#' query(mn, "A", evidence = list(B = "1"))
#' @export
markov_network <- function(factors = list(), nodes = NULL, states = NULL) {
  if (inherits(factors, "pgm_factor")) factors <- list(factors)
  factors <- lapply(factors, function(f) {
    if (!inherits(f, "pgm_factor")) pgm_abort("All potentials must be `pgm_factor` objects.")
    as_plain_factor(f)
  })
  reg <- factor_registry(factors)
  if (!is.null(states)) {
    for (v in names(states)) {
      if (!is.null(reg[[v]]) && !identical(reg[[v]]$states, as.character(states[[v]]))) {
        pgm_abort(sprintf("`states` disagrees with the potentials for '%s'.", v), class = "pgm_state_mismatch")
      }
      reg[[v]] <- list(type = "discrete", states = as.character(states[[v]]))
    }
  }
  vars <- unique(c(as.character(nodes %||% character()), names(reg)))
  miss <- setdiff(vars, names(reg))
  if (length(miss)) pgm_abort(paste0("States are unknown for: ", paste(miss, collapse = ", "), "."))
  mn <- structure(list(graph = NULL, factors = factors, registry = reg[vars]), class = c("markov_network", "pgm_model"))
  mn$graph <- mn_graph(vars, factors)
  mn
}

mn_graph <- function(vars, factors) {
  from <- to <- character()
  for (f in factors) {
    if (length(f$vars) > 1L) {
      cmb <- combn(f$vars, 2L)
      from <- c(from, cmb[1, ])
      to <- c(to, cmb[2, ])
    }
  }
  undirected_from_pairs(vars, from, to)
}

#' Add a potential to a Markov network or factor graph
#'
#' @param model A `markov_network` or `factor_graph`.
#' @param factor A `pgm_factor`.
#' @return The modified model.
#' @export
add_factor <- function(model, factor) {
  if (!inherits(factor, "pgm_factor")) pgm_abort("`factor` must be a `pgm_factor`.")
  if (inherits(model, "markov_network")) {
    return(markov_network(c(model$factors, list(factor)), nodes = model$graph$nodes,
                          states = model_states(model)))
  }
  if (inherits(model, "factor_graph")) {
    return(factor_graph(c(model$factors, list(factor)), states = model_states(model)))
  }
  pgm_abort("`add_factor()` needs a markov_network or factor_graph.")
}

#' @export
model_factors.markov_network <- function(model) model$factors

#' @export
check_model.markov_network <- function(model, error = TRUE) {
  probs <- character()
  for (i in seq_along(model$factors)) {
    f <- model$factors[[i]]
    if (!f$log && any(f$values < 0)) probs <- c(probs, sprintf("Potential %d has negative values.", i))
    if (any(is.na(f$values))) probs <- c(probs, sprintf("Potential %d has missing values.", i))
    for (v in f$vars) if (!identical(model$registry[[v]]$states, f$states[[v]])) {
      probs <- c(probs, sprintf("Potential %d uses states for '%s' that differ from the registry.", i, v))
    }
  }
  if (error && length(probs)) pgm_abort(paste(c("The model is invalid:", paste0("  ", probs)), collapse = "\n"), class = "pgm_invalid_model")
  if (error) invisible(TRUE) else probs
}

#' @export
print.markov_network <- function(x, ...) {
  cat(sprintf("<Markov network: %d variables, %d potentials, %d edges>\n",
              length(x$graph$nodes), length(x$factors), nrow(x$graph$edges)))
  for (f in x$factors) cat("  phi(", paste(f$vars, collapse = ", "), ")\n", sep = "")
  invisible(x)
}

#' Partition function
#'
#' Computes the normalizing constant Z, the sum over all joint assignments of
#' the product of potentials, by exact variable elimination in log space. For
#' a Bayesian network Z equals one up to rounding.
#'
#' @param model A Markov network, factor graph, Bayesian network or list of
#'   factors.
#' @param log Return log Z.
#' @param elimination_order Ordering heuristic or explicit order.
#' @return A numeric scalar.
#' @export
partition_function <- function(model, log = FALSE, elimination_order = "min_fill") {
  fs <- lapply(model_factors(model), factor_to_log)
  scope <- unique(unlist(lapply(fs, `[[`, "vars"), use.names = FALSE))
  ord <- resolve_order(fs, scope, elimination_order)
  res <- eliminate_vars(fs, ord, "sum")
  lz <- sum(vapply(res$factors, function(f) f$values, numeric(1)))
  # Variables in the registry without any potential contribute a factor of
  # their cardinality.
  if (!is.null(model$registry)) {
    free <- setdiff(names(model$registry), scope)
    lz <- lz + sum(log(vapply(model$registry[free], function(r) length(r$states), integer(1))))
  }
  if (log) lz else exp(lz)
}

#' Factor graphs
#'
#' A factor graph is the bipartite graph between variables and factors in
#' which each factor is joined to the variables in its scope. `factor_graph()`
#' builds one from factors; `as_factor_graph()` converts a discrete Bayesian
#' network (one factor per CPD) or Markov network (one factor per potential).
#' Both conversions preserve the joint distribution exactly; converting a
#' Bayesian network loses the arc directions, which are kept in the
#' `source_graph` element.
#'
#' @param factors A list of `pgm_factor` objects.
#' @param states Optional named list of states for variables with no factor.
#' @param model A model to convert.
#' @return An object of class `factor_graph`.
#' @examples
#' fg <- as_factor_graph(example_network("sprinkler"))
#' fg
#' loopy_bp(fg, evidence = list(WetGrass = "T"))$converged
#' @export
factor_graph <- function(factors = list(), states = NULL) {
  if (inherits(factors, "pgm_factor")) factors <- list(factors)
  factors <- lapply(factors, as_plain_factor)
  reg <- factor_registry(factors)
  for (v in names(states)) if (is.null(reg[[v]])) reg[[v]] <- list(type = "discrete", states = as.character(states[[v]]))
  vars <- names(reg)
  fnames <- paste0("f", seq_along(factors))
  names(factors) <- fnames
  from <- unlist(lapply(seq_along(factors), function(i) rep(fnames[i], length(factors[[i]]$vars))), use.names = FALSE)
  to <- unlist(lapply(factors, `[[`, "vars"), use.names = FALSE)
  structure(list(factors = factors, registry = reg, variables = vars,
                 graph = new_pgm_graph(c(vars, fnames), edge_frame(from, to, rep("undirected", length(from))), "pgm_ugraph")),
            class = c("factor_graph", "pgm_model"))
}

#' @rdname factor_graph
#' @export
as_factor_graph <- function(model) {
  if (inherits(model, "factor_graph")) return(model)
  fg <- factor_graph(unname(model_factors(model)), states = model_states(model))
  if (inherits(model, "bayesian_network")) fg$source_graph <- model$graph
  fg
}

#' Convert to a Markov network
#'
#' Converts a discrete Bayesian network to the Markov network over its moral
#' graph whose potentials are the CPDs. The joint distribution is preserved
#' exactly; the arc directions (and hence the independencies implied by
#' v-structures) are lost.
#'
#' @param model A discrete `bayesian_network` or `factor_graph`.
#' @return A `markov_network`.
#' @export
as_markov_network <- function(model) {
  if (inherits(model, "markov_network")) return(model)
  markov_network(unname(model_factors(model)), states = model_states(model))
}

#' @export
model_factors.factor_graph <- function(model) unname(model$factors)

#' @export
nodes.factor_graph <- function(x) x$variables

#' @export
check_model.factor_graph <- function(model, error = TRUE) check_model.markov_network(model, error)

#' @export
print.factor_graph <- function(x, ...) {
  cat(sprintf("<Factor graph: %d variables, %d factors>\n", length(x$variables), length(x$factors)))
  for (n in names(x$factors)) cat("  ", n, ": ", paste(x$factors[[n]]$vars, collapse = ", "), "\n", sep = "")
  invisible(x)
}
