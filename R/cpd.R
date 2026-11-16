# Tabular conditional probability distributions --------------------------------

cpd_tolerance <- function() getOption("pgmcompR.cpd_tolerance", 1e-8)

#' Tabular conditional probability distributions
#'
#' `cpd()` builds the conditional distribution P(variable | parents) as a table.
#' `probabilities` may be
#'
#' * a named numeric vector (no parents): names are the states;
#' * a matrix with one row per state of `variable` and one column per parent
#'   configuration, the first parent varying fastest;
#' * an array whose first dimension is `variable` and whose remaining
#'   dimensions are the parents, with named dimnames;
#' * a plain numeric vector of the right length (child fastest), with
#'   `states` supplied.
#'
#' Every column must sum to one within `getOption("pgmcompR.cpd_tolerance",
#' 1e-8)`; use `normalize = TRUE` to rescale columns instead. A column of
#' zeros is never silently replaced.
#'
#' @param variable Name of the child variable.
#' @param probabilities Probabilities in one of the formats above.
#' @param parents Character vector of parent names (taken from the array's
#'   dimnames when omitted).
#' @param states Named list of state labels for `variable` and the parents,
#'   used when they cannot be read from names or dimnames.
#' @param normalize Logical; rescale each column to sum to one.
#' @return An object of class `pgm_cpd` (also a `pgm_factor` whose first
#'   variable is the child).
#' @examples
#' cpd("Smoking", c(no = 0.7, yes = 0.3))
#' cpd("Cancer", matrix(c(0.99, 0.01, 0.9, 0.1), 2,
#'                      dimnames = list(Cancer = c("no", "yes"), NULL)),
#'     parents = "Smoking", states = list(Smoking = c("no", "yes")))
#' @export
cpd <- function(variable, probabilities, parents = NULL, states = NULL, normalize = FALSE) {
  check_string(variable)
  states <- states %||% list()
  if (is.array(probabilities) && length(dim(probabilities)) >= 2L && !is.matrix(probabilities)) {
    dn <- dimnames(probabilities)
    nm <- names(dn)
    if (is.null(nm) || any(!nzchar(nm))) {
      if (is.null(parents)) pgm_abort("An array of probabilities needs named dimnames or `parents`.")
      nm <- c(variable, parents)
    }
    if (nm[1L] != variable) pgm_abort("The first dimension of the array must be the child variable.")
    if (is.null(parents)) parents <- nm[-1L]
    if (!identical(nm[-1L], parents)) {
      pgm_abort("The array dimensions must be ordered as the child followed by `parents`.")
    }
    for (k in seq_along(nm)) if (!is.null(dn[[k]])) states[[nm[k]]] <- dn[[k]]
    vals <- as.numeric(probabilities)
    dims <- dim(probabilities)
  } else if (is.matrix(probabilities)) {
    parents <- parents %||% character()
    rn <- rownames(probabilities)
    # Ignore partial row names such as those created by rbind(1 - p, p).
    if (!is.null(rn) && all(nzchar(rn)) && is.null(states[[variable]])) states[[variable]] <- rn
    vals <- as.numeric(probabilities)
    dims <- dim(probabilities)
  } else if (is.numeric(probabilities)) {
    parents <- parents %||% character()
    if (!is.null(names(probabilities)) && !length(parents)) states[[variable]] <- names(probabilities)
    vals <- as.numeric(probabilities)
    dims <- NULL
  } else {
    pgm_abort("`probabilities` must be numeric (vector, matrix or array).")
  }
  check_names(c(variable, parents), "c(variable, parents)")
  vars <- c(variable, parents)
  missing <- setdiff(vars, names(states))
  if (length(missing)) {
    pgm_abort(paste0("State labels are unknown for: ", paste(missing, collapse = ", "),
                     ". Supply them through `states` or dimnames."))
  }
  st <- normalize_states(vars, states[vars])
  card <- vapply(st, length, integer(1))
  if (!is.null(dims) && length(dims) == 2L && length(parents) && dims[2L] != prod(card[-1L])) {
    pgm_abort(sprintf("The matrix has %d columns but the parents have %d configurations.", dims[2L], prod(card[-1L])))
  }
  if (!is.null(dims) && dims[1L] != card[1L]) {
    pgm_abort(sprintf("The table has %d rows but '%s' has %d states.", dims[1L], variable, card[1L]))
  }
  if (length(vals) != prod(card)) {
    pgm_abort(sprintf("`probabilities` has %d values; expected %d (%s).", length(vals), prod(card),
                      paste(card, collapse = " x ")))
  }
  check_prob(vals, "probabilities")
  m <- matrix(vals, nrow = card[1L])
  cs <- colSums(m)
  if (any(cs == 0)) {
    pgm_abort(sprintf("Column(s) %s of the CPD of '%s' are all zero: a conditional distribution must be specified for every parent configuration.",
                      paste(which(cs == 0), collapse = ", "), variable), class = "pgm_invalid_cpd")
  }
  if (normalize) {
    m <- sweep(m, 2L, cs, "/")
  } else if (any(abs(cs - 1) > cpd_tolerance())) {
    bad <- which(abs(cs - 1) > cpd_tolerance())
    pgm_abort(sprintf("The CPD of '%s' is not normalized: column(s) %s sum to %s. Use `normalize = TRUE` to rescale.",
                      variable, paste(head(bad, 5), collapse = ", "), paste(format(head(cs[bad], 5)), collapse = ", ")),
              class = "pgm_invalid_cpd")
  }
  new_cpd(variable, parents, st, as.numeric(m))
}

new_cpd <- function(variable, parents, states, values) {
  f <- new_factor(c(variable, parents), states, values, FALSE, class = "pgm_cpd")
  f$child <- variable
  f$parents <- parents
  f
}

#' Parametric discrete CPDs
#'
#' Convenience constructors that build a tabular CPD from a parametric rule.
#' `noisy_or_cpd()` gives a binary child that is "on" unless every active
#' parent is independently inhibited: P(off | parents) = (1 - leak) times the
#' product of (1 - p_i) over active parents. `deterministic_cpd()` puts
#' probability one on `fun(parent states)`. `custom_cpd()` calls `fun` on each
#' parent configuration (a named character vector) to obtain the child's
#' distribution.
#'
#' @param variable Child name.
#' @param parents Parent names.
#' @param states Named list of states for the child and parents.
#' @param probs Named numeric vector of activation probabilities per parent.
#' @param active Named list giving, for each parent, the state(s) that count as
#'   active. Defaults to the last state.
#' @param leak Leak probability.
#' @param fun A function of a named character vector of parent states.
#' @return A `pgm_cpd`.
#' @examples
#' st <- list(Fever = c("no", "yes"), Flu = c("no", "yes"), Cold = c("no", "yes"))
#' noisy_or_cpd("Fever", c("Flu", "Cold"), st, probs = c(Flu = 0.8, Cold = 0.4), leak = 0.05)
#' @export
noisy_or_cpd <- function(variable, parents, states, probs, active = NULL, leak = 0) {
  st <- normalize_states(c(variable, parents), states)
  if (length(st[[variable]]) != 2L) pgm_abort("A noisy-OR child must be binary.")
  if (is.null(names(probs)) || !setequal(names(probs), parents)) pgm_abort("`probs` must be named by `parents`.")
  if (any(probs < 0 | probs > 1) || leak < 0 || leak > 1) pgm_abort("Probabilities must lie in [0, 1].")
  active <- active %||% lapply(st[parents], function(s) s[length(s)])
  custom_cpd(variable, parents, st, function(pa) {
    off <- 1 - leak
    for (p in parents) if (pa[[p]] %in% active[[p]]) off <- off * (1 - probs[[p]])
    c(off, 1 - off)
  })
}

#' @rdname noisy_or_cpd
#' @export
deterministic_cpd <- function(variable, parents, states, fun) {
  st <- normalize_states(c(variable, parents), states)
  custom_cpd(variable, parents, st, function(pa) {
    val <- as.character(fun(pa))
    if (length(val) != 1L || !val %in% st[[variable]]) pgm_abort("`fun` must return one state of the child.")
    as.numeric(st[[variable]] == val)
  })
}

#' @rdname noisy_or_cpd
#' @export
custom_cpd <- function(variable, parents, states, fun) {
  st <- normalize_states(c(variable, parents), states)
  if (!length(parents)) {
    return(cpd(variable, fun(setNames(character(), character())), states = st))
  }
  grid <- expand.grid(st[parents], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  cols <- lapply(seq_len(nrow(grid)), function(i) {
    p <- fun(unlist(grid[i, , drop = FALSE]))
    if (!is.numeric(p) || length(p) != length(st[[variable]])) {
      pgm_abort("`fun` must return one probability per child state.")
    }
    p
  })
  cpd(variable, do.call(cbind, cols), parents = parents, states = st)
}

#' CPD components
#'
#' @param x A CPD.
#' @return `cpd_variable()` returns the child name and `cpd_parents()` the
#'   parent names.
#' @export
cpd_variable <- function(x) x$child %||% x$variable

#' @rdname cpd_variable
#' @export
cpd_parents <- function(x) x$parents %||% character()

#' @export
print.pgm_cpd <- function(x, ...) {
  cat(sprintf("<CPD: P(%s%s)>\n", x$child, if (length(x$parents)) paste0(" | ", paste(x$parents, collapse = ", ")) else ""))
  m <- matrix(x$values, nrow = x$card[[1L]], dimnames = list(x$states[[1L]], NULL))
  if (length(x$parents)) {
    grid <- expand.grid(x$states[x$parents], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    colnames(m) <- apply(grid, 1L, function(r) paste(paste0(x$parents, "=", r), collapse = ","))
  } else colnames(m) <- "P"
  names(dimnames(m)) <- c(x$child, if (length(x$parents)) "" else "")
  if (ncol(m) > 8L) {
    print(round(m[, 1:8, drop = FALSE], 6))
    cat(sprintf("... %d more parent configurations\n", ncol(m) - 8L))
  } else print(round(m, 6))
  invisible(x)
}
