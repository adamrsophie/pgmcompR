# Bayesian networks -------------------------------------------------------------

#' Bayesian networks
#'
#' Creates a Bayesian network from a DAG and, optionally, its conditional
#' distributions. Conditional distributions may be tabular ([cpd()]), linear
#' Gaussian ([lg_cpd()]), conditional Gaussian ([cg_cpd()]) or any custom type
#' that follows the CPD contract (see [cpd_contract]).
#'
#' The variable registry records, for every node, whether it is discrete or
#' continuous and, for discrete nodes, its ordered state labels. States can be
#' declared up front with `states`, or are learned from the CPDs as they are
#' added. Every CPD must use the registered state order of its variables.
#'
#' @param edges Edges (see [dag()]) or a `pgm_dag`.
#' @param nodes Optional node names, for isolated nodes.
#' @param cpds Optional list of CPDs.
#' @param states Optional named list of state labels for discrete nodes.
#' @return An object of class `bayesian_network`.
#' @examples
#' bn <- bayesian_network(list(c("Smoking", "Cancer"), c("Cancer", "XRay")))
#' bn <- add_cpd(bn, "Smoking", c(no = 0.7, yes = 0.3))
#' bn <- add_cpd(bn, "Cancer", matrix(c(0.99, 0.01, 0.9, 0.1), 2,
#'               dimnames = list(c("no", "yes"), NULL)))
#' bn <- add_cpd(bn, "XRay", matrix(c(0.95, 0.05, 0.2, 0.8), 2,
#'               dimnames = list(c("negative", "positive"), NULL)))
#' check_model(bn)
#' bn
#' @export
bayesian_network <- function(edges = NULL, nodes = NULL, cpds = NULL, states = NULL) {
  g <- if (inherits(edges, "pgm_graph")) {
    if (!inherits(edges, "pgm_dag")) pgm_abort("A Bayesian network needs a DAG.")
    add_nodes(edges, as.character(nodes %||% character()))
  } else dag(edges, nodes)
  bn <- structure(list(graph = g, cpds = list(), registry = list()), class = c("bayesian_network", "pgm_model"))
  if (!is.null(states)) {
    if (!is.list(states) || is.null(names(states))) pgm_abort("`states` must be a named list.")
    check_nodes_in(g, names(states), "states")
    for (v in names(states)) {
      bn$registry[[v]] <- list(type = "discrete", states = normalize_states(v, states[v])[[1L]], declared = TRUE)
    }
  }
  for (c in cpds %||% list()) bn <- add_cpd(bn, c)
  bn
}

is_discrete_cpd <- function(x) inherits(x, "pgm_cpd")

cpd_kind <- function(x) {
  if (inherits(x, "pgm_cpd")) "discrete"
  else if (inherits(x, "pgm_lg_cpd")) "gaussian"
  else if (inherits(x, "pgm_cg_cpd")) "conditional_gaussian"
  else "custom"
}

cpd_registry_entries <- function(x) {
  if (inherits(x, "pgm_cpd")) {
    return(lapply(setNames(x$vars, x$vars), function(v) list(type = "discrete", states = x$states[[v]])))
  }
  if (!inherits(x, c("pgm_lg_cpd", "pgm_cg_cpd"))) return(x$registry_entries %||% list())
  ent <- list()
  ent[[cpd_variable(x)]] <- list(type = "continuous", states = NULL)
  if (!is.null(x$discrete_parents)) {
    for (v in x$discrete_parents) ent[[v]] <- list(type = "discrete", states = x$discrete_states[[v]])
  }
  if (!is.null(x$continuous_parents)) {
    for (v in x$continuous_parents) ent[[v]] <- list(type = "continuous", states = NULL)
  }
  if (!is.null(x$registry_entries)) ent <- utils::modifyList(ent, x$registry_entries)
  ent
}

#' Add, remove and extract CPDs
#'
#' `add_cpd()` attaches a conditional distribution to a node, replacing any
#' existing one. Either pass a CPD object, or pass `variable` and
#' `probabilities` to build a tabular CPD whose parents are the node's parents
#' in the graph (in graph order) and whose parent states come from the
#' registry. `remove_cpd()` detaches it and `get_cpd()` returns it. `cpds()`
#' returns all CPDs.
#'
#' @param model A `bayesian_network`.
#' @param cpd A CPD object, or the name of the child when `probabilities` is
#'   given.
#' @param probabilities Probabilities passed to [cpd()].
#' @param states Optional state labels passed to [cpd()].
#' @param ... Passed to [cpd()].
#' @param variable Node name.
#' @return `add_cpd()` and `remove_cpd()` return the modified model;
#'   `get_cpd()` returns a CPD; `cpds()` a named list.
#' @export
add_cpd <- function(model, cpd, probabilities = NULL, states = NULL, ...) {
  if (!inherits(model, "bayesian_network")) pgm_abort("`model` must be a bayesian_network.")
  if (is.character(cpd)) {
    variable <- cpd
    check_string(variable)
    check_nodes_in(model$graph, variable, "variable")
    pa <- parents(model$graph, variable)
    st <- states %||% list()
    for (p in pa) if (is.null(st[[p]]) && !is.null(model$registry[[p]]$states)) st[[p]] <- model$registry[[p]]$states
    named_array <- is.array(probabilities) && !is.matrix(probabilities) && !is.null(names(dimnames(probabilities)))
    child_labeled <- (named_array && !is.null(dimnames(probabilities)[[1L]])) ||
      (is.matrix(probabilities) && !is.null(rownames(probabilities))) ||
      (is.null(dim(probabilities)) && !is.null(names(probabilities)) && !length(pa))
    if (is.null(st[[variable]]) && !child_labeled && !is.null(model$registry[[variable]]$states)) {
      st[[variable]] <- model$registry[[variable]]$states
    }
    cpd <- cpd(variable, probabilities, parents = if (named_array) NULL else pa, states = st, ...)
  }
  v <- cpd_variable(cpd)
  if (is.null(v)) pgm_abort("The CPD does not name its variable.")
  check_nodes_in(model$graph, v, "variable")
  gp <- parents(model$graph, v)
  cp <- cpd_parents(cpd)
  if (!setequal(gp, cp) || length(gp) != length(cp)) {
    pgm_abort(sprintf("The CPD of '%s' has parents {%s} but the graph has parents {%s}.", v,
                      paste(cp, collapse = ", "), paste(gp, collapse = ", ")), class = "pgm_invalid_cpd")
  }
  ent <- cpd_registry_entries(cpd)
  for (nm in names(ent)) {
    old <- model$registry[[nm]]
    new <- ent[[nm]]
    if (!is.null(old)) {
      if (!identical(old$type, new$type)) {
        pgm_abort(sprintf("'%s' is registered as %s but the CPD treats it as %s.", nm, old$type, new$type),
                  class = "pgm_state_mismatch")
      }
      if (identical(new$type, "discrete") && !identical(old$states, new$states)) {
        used <- any(vapply(model$cpds[names(model$cpds) != v], function(c) nm %in% cpd_parents(c), logical(1)))
        if (nm == v && !used && !isTRUE(old$declared)) next

        pgm_abort(sprintf("The CPD uses states {%s} for '%s' but the model has {%s}.",
                          paste(new$states, collapse = ", "), nm, paste(old$states, collapse = ", ")),
                  class = "pgm_state_mismatch")
      }
    }
  }
  for (nm in names(ent)) {
    if (is.null(model$registry[[nm]])) model$registry[[nm]] <- ent[[nm]]
    else if (nm == v) model$registry[[nm]]$states <- ent[[nm]]$states
  }
  model$cpds[[v]] <- cpd
  model$cpds <- model$cpds[intersect(model$graph$nodes, names(model$cpds))]
  model
}

#' @rdname add_cpd
#' @export
remove_cpd <- function(model, variable) {
  check_string(variable)
  if (is.null(model$cpds[[variable]])) pgm_abort(sprintf("'%s' has no CPD.", variable))
  model$cpds[[variable]] <- NULL
  model
}

#' @rdname add_cpd
#' @export
get_cpd <- function(model, variable) {
  check_string(variable)
  out <- model$cpds[[variable]]
  if (is.null(out)) pgm_abort(sprintf("'%s' has no CPD.", variable))
  out
}

#' @rdname add_cpd
#' @export
cpds <- function(model) model$cpds

#' Variable registry
#'
#' Returns the registry of a model: for each variable its type
#' (`"discrete"` or `"continuous"`) and its state labels.
#'
#' @param model A model.
#' @return A data frame with columns `variable`, `type`, `cardinality` and
#'   `states`.
#' @export
variable_registry <- function(model) {
  reg <- model$registry
  vars <- intersect(nodes(model), names(reg))
  data.frame(
    variable = vars,
    type = vapply(reg[vars], function(r) r$type, ""),
    cardinality = vapply(reg[vars], function(r) if (is.null(r$states)) NA_integer_ else length(r$states), integer(1)),
    states = vapply(reg[vars], function(r) paste(r$states, collapse = ", "), ""),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

#' @rdname variable_registry
#' @return `model_states()` returns a named list of the state labels of the
#'   discrete variables.
#' @export
model_states <- function(model) {
  reg <- model$registry
  out <- lapply(reg, function(r) r$states)
  out[!vapply(out, is.null, logical(1))]
}

#' @export
print.bayesian_network <- function(x, ...) {
  kinds <- vapply(x$cpds, cpd_kind, "")
  type <- if (!length(kinds)) "empty" else if (all(kinds == "discrete")) "discrete"
  else if (all(kinds == "gaussian")) "linear Gaussian" else if (all(kinds %in% c("discrete", "conditional_gaussian", "gaussian"))) "hybrid" else "custom"
  cat(sprintf("<Bayesian network (%s): %d nodes, %d arcs, %d of %d CPDs>\n", type,
              length(x$graph$nodes), nrow(x$graph$edges), length(x$cpds), length(x$graph$nodes)))
  for (v in x$graph$nodes) {
    pa <- parents(x$graph, v)
    r <- x$registry[[v]]
    st <- if (!is.null(r$states)) paste0(" [", paste(r$states, collapse = ", "), "]") else if (identical(r$type, "continuous")) " [continuous]" else ""
    cat(sprintf("  %s%s%s%s\n", v, if (length(pa)) paste0(" | ", paste(pa, collapse = ", ")) else "", st,
                if (is.null(x$cpds[[v]])) "  (no CPD)" else ""))
  }
  invisible(x)
}

#' Validate a model
#'
#' Checks that every node has a CPD, that each CPD's parents match the graph,
#' that state labels agree across CPDs, and that every conditional
#' distribution is normalized and nonnegative. Gaussian CPDs must have
#' positive residual standard deviations. Markov networks and factor graphs
#' are checked for factor scopes and nonnegative potentials.
#'
#' @param model A model.
#' @param error If `TRUE` (default) problems raise an error; otherwise they
#'   are returned.
#' @return Invisibly `TRUE` when valid; with `error = FALSE`, a character
#'   vector of problems (empty when valid).
#' @export
check_model <- function(model, error = TRUE) UseMethod("check_model")

#' @export
check_model.bayesian_network <- function(model, error = TRUE) {
  probs <- character()
  g <- model$graph
  for (v in g$nodes) {
    c <- model$cpds[[v]]
    if (is.null(c)) {
      probs <- c(probs, sprintf("'%s' has no CPD.", v))
      next
    }
    pa <- parents(g, v)
    if (!setequal(pa, cpd_parents(c))) probs <- c(probs, sprintf("CPD parents of '%s' do not match the graph.", v))
    if (inherits(c, "pgm_cpd")) {
      for (u in c$vars) {
        r <- model$registry[[u]]
        if (is.null(r) || !identical(r$states, c$states[[u]])) {
          probs <- c(probs, sprintf("CPD of '%s' uses states for '%s' that differ from the registry.", v, u))
        }
        if (!is.null(model$cpds[[u]]) && !inherits(model$cpds[[u]], "pgm_cpd") && !inherits(model$cpds[[u]], "pgm_custom_discrete")) {
          if (u != v) probs <- c(probs, sprintf("Discrete CPD of '%s' has a continuous parent '%s'.", v, u))
        }
      }
      m <- matrix(c$values, nrow = c$card[[1L]])
      if (any(m < 0)) probs <- c(probs, sprintf("CPD of '%s' has negative entries.", v))
      if (any(abs(colSums(m) - 1) > cpd_tolerance())) probs <- c(probs, sprintf("CPD of '%s' is not normalized.", v))
    } else {
      chk <- tryCatch(validate_cpd(c, model), error = function(e) conditionMessage(e))
      if (is.character(chk)) probs <- c(probs, chk)
    }
  }
  if (error && length(probs)) pgm_abort(paste(c("The model is invalid:", paste0("  ", probs)), collapse = "\n"), class = "pgm_invalid_model")
  if (error) invisible(TRUE) else probs
}

#' Validate a CPD against a model
#'
#' Generic used by [check_model()] for non tabular CPDs. Methods should signal
#' an error describing the problem, or return `TRUE`.
#'
#' @param x A CPD.
#' @param model The model it belongs to.
#' @param ... Unused.
#' @return `TRUE` invisibly.
#' @export
validate_cpd <- function(x, model, ...) UseMethod("validate_cpd")

#' @export
validate_cpd.default <- function(x, model, ...) invisible(TRUE)

# Discrete factors of a model ---------------------------------------------------

#' Factors of a model
#'
#' Returns the discrete factors that define a model: the CPD tables of a
#' discrete Bayesian network (custom discrete CPDs are tabulated with
#' [cpd_to_factor()]), or the potentials of a Markov network or factor graph.
#'
#' @param model A model.
#' @return A list of `pgm_factor` objects.
#' @export
model_factors <- function(model) UseMethod("model_factors")

#' @export
model_factors.bayesian_network <- function(model) {
  check_model(model)
  bn_factors(model)
}

bn_factors <- function(model) {
  lapply(model$cpds, function(c) {
    if (inherits(c, "pgm_cpd")) return(c)
    f <- cpd_to_factor(c, model)
    if (!inherits(f, "pgm_factor")) pgm_abort(sprintf("The CPD of '%s' is not discrete.", cpd_variable(c)))
    f
  })
}

#' @export
model_factors.default <- function(model) {
  if (is.list(model) && all(vapply(model, inherits, logical(1), "pgm_factor"))) return(model)
  pgm_abort("Cannot extract factors from this object.")
}

#' Tabulate a CPD
#'
#' Generic that converts a discrete CPD into a `pgm_factor` over the child
#' and its parents, with the child first. Custom discrete CPD classes must
#' provide a method.
#'
#' @param x A CPD.
#' @param model The model (for state labels).
#' @param ... Unused.
#' @return A `pgm_factor`.
#' @export
cpd_to_factor <- function(x, model, ...) UseMethod("cpd_to_factor")

#' @export
cpd_to_factor.pgm_cpd <- function(x, model, ...) x

#' @export
cpd_to_factor.default <- function(x, model, ...) {
  pgm_abort(sprintf("No `cpd_to_factor()` method for class '%s'.", class(x)[1L]))
}

#' Model type predicates
#'
#' `is_discrete_bn()` is `TRUE` when every variable is discrete,
#' `is_gaussian_bn()` when every CPD is linear Gaussian and `is_hybrid_bn()`
#' when the network mixes discrete and continuous variables.
#'
#' @param model A Bayesian network.
#' @return A logical scalar.
#' @examples
#' is_discrete_bn(example_network("asia"))
#' is_gaussian_bn(example_network("gaussian"))
#' @export
is_discrete_bn <- function(model) {
  inherits(model, "bayesian_network") && length(model$cpds) &&
    all(vapply(model$registry[model$graph$nodes], function(r) identical(r$type, "discrete"), logical(1)))
}

#' @rdname is_discrete_bn
#' @export
is_gaussian_bn <- function(model) {
  inherits(model, "bayesian_network") && length(model$cpds) &&
    all(vapply(model$cpds, inherits, logical(1), "pgm_lg_cpd"))
}

#' @rdname is_discrete_bn
#' @export
is_hybrid_bn <- function(model) {
  inherits(model, "bayesian_network") && length(model$cpds) && !is_discrete_bn(model) && !is_gaussian_bn(model)
}
