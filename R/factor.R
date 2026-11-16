# Discrete factors --------------------------------------------------------------
#
# A factor stores values as a plain numeric vector in column major order: the
# first variable varies fastest, exactly like an R array with dim = card.

new_factor <- function(vars, states, values, log = FALSE, class = NULL) {
  card <- vapply(states, length, integer(1))
  names(card) <- vars
  names(states) <- vars
  structure(list(vars = vars, card = card, states = states,
                 values = as.numeric(values), log = isTRUE(log)),
            class = c(class, "pgm_factor"))
}

#' Discrete factors
#'
#' Creates a discrete factor (a nonnegative function of a set of discrete
#' variables). Values are given in column major order: the first variable
#' varies fastest, as in an R array. `values` may also be an array whose
#' dimnames name the variables and states; then `variables` and `states` are
#' taken from it. Lengths are checked exactly: values are never recycled.
#'
#' The name `discrete_factor()` is used instead of `factor()` so that base
#' R's `factor()` is not masked.
#'
#' @param variables Character vector of variable names.
#' @param states Named list of state labels per variable, or an integer vector
#'   of cardinalities (states are then labeled `"1"`, `"2"`, ...).
#' @param values Numeric vector or array of nonnegative values.
#' @param log If `TRUE`, `values` are log potentials (may contain `-Inf`).
#' @return An object of class `pgm_factor`.
#' @examples
#' f <- discrete_factor(c("A", "B"), list(A = c("a0", "a1"), B = c("b0", "b1")),
#'                      c(0.3, 0.7, 0.9, 0.1))
#' f
#' as.data.frame(f)
#' @export
discrete_factor <- function(variables = NULL, states = NULL, values, log = FALSE) {
  if (is.array(values) && !is.null(names(dimnames(values))) && is.null(variables)) {
    dn <- dimnames(values)
    variables <- names(dn)
    if (any(vapply(dn, is.null, logical(1)))) pgm_abort("Every dimension of `values` needs state labels.")
    states <- dn
  }
  if (is.null(variables)) variables <- character()
  check_names(variables, "variables")
  states <- normalize_states(variables, states)
  card <- vapply(states, length, integer(1))
  expected <- prod(card)
  if (length(values) != expected) {
    pgm_abort(sprintf("`values` has length %d but the variables' cardinalities imply %d.",
                      length(values), expected))
  }
  if (is.array(values) && length(dim(values)) == length(card) && !all(dim(values) == card)) {
    pgm_abort("The dimensions of `values` do not match the cardinalities.")
  }
  values <- as.numeric(values)
  if (isTRUE(log)) {
    if (anyNA(values) || any(values == Inf)) pgm_abort("Log values must be finite or -Inf.")
  } else {
    check_prob(values, "values")
  }
  new_factor(variables, states, values, log)
}

normalize_states <- function(variables, states) {
  if (is.null(states)) pgm_abort("`states` is required.")
  if (is.numeric(states)) {
    if (length(states) != length(variables)) pgm_abort("`states` cardinalities must match `variables`.")
    states <- lapply(states, function(k) as.character(seq_len(check_count(k, "states"))))
    names(states) <- variables
  }
  if (!is.list(states)) pgm_abort("`states` must be a named list or integer vector.")
  if (is.null(names(states))) {
    if (length(states) != length(variables)) pgm_abort("`states` must have one element per variable.")
    names(states) <- variables
  }
  missing <- setdiff(variables, names(states))
  if (length(missing)) pgm_abort(paste0("Missing states for: ", paste(missing, collapse = ", "), "."))
  states <- lapply(states[variables], function(s) {
    s <- as.character(s)
    if (!length(s)) pgm_abort("Each variable needs at least one state.")
    check_names(s, "states")
    s
  })
  states
}

#' @export
print.pgm_factor <- function(x, n = 20L, ...) {
  kind <- if (inherits(x, "pgm_cpd")) "CPD" else "Factor"
  if (length(x$vars) == 0L) {
    cat(sprintf("<%s: scalar> %s\n", kind, format(x$values)))
    return(invisible(x))
  }
  label <- if (inherits(x, "pgm_cpd")) {
    paste0("P(", x$child, if (length(x$parents)) paste0(" | ", paste(x$parents, collapse = ", ")), ")")
  } else paste(x$vars, collapse = ", ")
  cat(sprintf("<%s: %s>%s\n", kind, label, if (x$log) " (log scale)" else ""))
  df <- as.data.frame(x)
  if (nrow(df) > n) {
    print(utils::head(df, n), row.names = FALSE)
    cat(sprintf("... %d more rows\n", nrow(df) - n))
  } else {
    print(df, row.names = FALSE)
  }
  invisible(x)
}

#' @export
as.data.frame.pgm_factor <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (!length(x$vars)) return(data.frame(value = x$values))
  grid <- expand.grid(x$states, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid[[if (x$log) "log_value" else "value"]] <- x$values
  grid
}

#' @export
as.array.pgm_factor <- function(x, ...) {
  if (!length(x$vars)) return(x$values)
  array(x$values, dim = unname(x$card), dimnames = x$states)
}

#' Factor accessors
#'
#' `factor_scope()` returns the variables of a factor in storage order,
#' `factor_states()` the state labels and `factor_values()` the values in
#' column major order.
#'
#' @param f A `pgm_factor`.
#' @return Character vector, named list, or numeric vector.
#' @export
factor_scope <- function(f) f$vars

#' @rdname factor_scope
#' @export
factor_states <- function(f) f$states

#' @rdname factor_scope
#' @export
factor_values <- function(f) f$values

#' Scalar factor
#'
#' @param value A single nonnegative number (or log value when `log = TRUE`).
#' @param log Logical.
#' @return A `pgm_factor` with empty scope.
#' @export
scalar_factor <- function(value = 1, log = FALSE) {
  new_factor(character(), list(), value, log)
}

#' Change factor scale
#'
#' `factor_to_log()` converts a factor to log scale and `factor_to_linear()`
#' back to the linear scale.
#'
#' @param f A `pgm_factor`.
#' @return A `pgm_factor`.
#' @export
factor_to_log <- function(f) {
  if (f$log) return(f)
  f$values <- log(f$values)
  f$log <- TRUE
  f
}

#' @rdname factor_to_log
#' @export
factor_to_linear <- function(f) {
  if (!f$log) return(f)
  f$values <- exp(f$values)
  f$log <- FALSE
  f
}

# Strip CPD class so algebra results are plain factors.
as_plain_factor <- function(f) {
  new_factor(f$vars, f$states, f$values, f$log)
}

check_compatible <- function(a, b) {
  shared <- intersect(a$vars, b$vars)
  for (v in shared) {
    if (!identical(a$states[[v]], b$states[[v]])) {
      if (setequal(a$states[[v]], b$states[[v]])) {
        pgm_abort(sprintf("Variable '%s' has the same states in a different order in the two factors; reorder states first.", v),
                  class = "pgm_state_mismatch")
      }
      pgm_abort(sprintf("Variable '%s' has incompatible states: {%s} versus {%s}.", v,
                        paste(a$states[[v]], collapse = ", "), paste(b$states[[v]], collapse = ", ")),
                class = "pgm_state_mismatch")
    }
  }
  invisible(TRUE)
}

# Zero-based index into a factor with variables `fvars` (cardinalities
# `fcard`) for every cell of the full grid over `uvars` (cardinalities `ucard`).
expand_index <- function(ucard, uvars, fvars, fcard) {
  n <- prod(ucard)
  if (!length(fvars)) return(integer(n))
  backend_expand_index(as.integer(ucard), match(fvars, uvars) - 1L, as.integer(fcard))
}

expand_index_r <- function(ucard, pos, fcard) {
  n <- prod(ucard)
  ustride <- c(1, cumprod(ucard))[seq_along(ucard)]
  fstride <- c(1, cumprod(fcard))[seq_along(fcard)]
  idx <- numeric(n)
  cell <- 0:(n - 1)
  for (k in seq_along(pos)) {
    j <- pos[k] + 1L
    idx <- idx + ((cell %/% ustride[j]) %% ucard[j]) * fstride[k]
  }
  idx
}

#' Factor product
#'
#' Multiplies factors pointwise over the union of their scopes. Shared
#' variables must have identical state labels in the same order. The scope of
#' the result lists the variables of the first factor, then new variables in
#' order of appearance. Log scale factors are added instead.
#'
#' @param ... `pgm_factor` objects, or a single list of them.
#' @return A `pgm_factor`.
#' @examples
#' a <- discrete_factor("A", list(A = c("0", "1")), c(0.4, 0.6))
#' b <- discrete_factor(c("B", "A"), list(B = c("0", "1"), A = c("0", "1")), c(0.9, 0.1, 0.2, 0.8))
#' factor_product(a, b)
#' @export
factor_product <- function(...) {
  fs <- list(...)
  if (length(fs) == 1L && is.list(fs[[1L]]) && !inherits(fs[[1L]], "pgm_factor")) fs <- fs[[1L]]
  if (!length(fs)) return(scalar_factor(1))
  for (f in fs) if (!inherits(f, "pgm_factor")) pgm_abort("All arguments must be `pgm_factor` objects.")
  use_log <- any(vapply(fs, function(f) f$log, logical(1)))
  if (use_log) fs <- lapply(fs, factor_to_log)
  out <- as_plain_factor(fs[[1L]])
  for (f in fs[-1L]) out <- product2(out, f)
  out
}

product2 <- function(a, b) {
  check_compatible(a, b)
  uvars <- c(a$vars, setdiff(b$vars, a$vars))
  ustates <- c(a$states, b$states[setdiff(b$vars, a$vars)])
  ucard <- vapply(ustates, length, integer(1))
  n <- prod(ucard)
  va <- rep_len(a$values, n) # a's scope is a prefix of the union: exact tiling
  vb <- b$values[expand_index(ucard, uvars, b$vars, b$card) + 1]
  vals <- if (a$log) va + vb else va * vb
  if (a$log) vals[is.nan(vals)] <- -Inf
  new_factor(uvars, ustates, vals, a$log)
}

#' Reorder factor variables
#'
#' Permutes the storage order of a factor's variables without changing the
#' function it represents.
#'
#' @param f A `pgm_factor`.
#' @param variables The new variable order (a permutation of the scope).
#' @return A `pgm_factor`.
#' @export
factor_reorder <- function(f, variables) {
  if (!setequal(variables, f$vars) || length(variables) != length(f$vars)) {
    pgm_abort("`variables` must be a permutation of the factor scope.")
  }
  if (identical(variables, f$vars)) return(f)
  perm <- match(variables, f$vars)
  arr <- array(f$values, dim = unname(f$card))
  out <- f
  out$values <- as.numeric(aperm(arr, perm))
  out$vars <- variables
  out$card <- f$card[variables]
  out$states <- f$states[variables]
  if (inherits(f, "pgm_cpd")) out <- as_plain_factor(out)
  out
}

#' Reorder the states of a variable in a factor
#'
#' @param f A `pgm_factor`.
#' @param variable Variable name.
#' @param states The new state order (a permutation of the current states).
#' @return A `pgm_factor`.
#' @export
factor_reorder_states <- function(f, variable, states) {
  check_string(variable)
  if (!variable %in% f$vars) pgm_abort(sprintf("'%s' is not in the factor scope.", variable))
  old <- f$states[[variable]]
  if (!setequal(old, states) || length(old) != length(states)) pgm_abort("`states` must be a permutation of the current states.")
  arr <- array(f$values, dim = unname(f$card))
  k <- match(variable, f$vars)
  idx <- lapply(f$card, seq_len)
  idx[[k]] <- match(states, old)
  f$values <- as.numeric(do.call(`[`, c(list(arr), unname(idx), drop = FALSE)))
  f$states[[variable]] <- states
  f
}

reduce_matrix <- function(f, drop_vars, fun = c("sum", "max")) {
  fun <- match.arg(fun)
  keep <- setdiff(f$vars, drop_vars)
  g <- factor_reorder(as_plain_factor(f), c(keep, intersect(f$vars, drop_vars)))
  nk <- prod(g$card[keep])
  m <- matrix(g$values, nrow = nk)
  list(keep = keep, m = m, g = g)
}

#' Marginalize and maximize factors
#'
#' `factor_marginalize()` sums the listed variables out of a factor (log-sum-exp
#' on the log scale). `factor_maximize()` maximizes them out; with
#' `argmax = TRUE` the maximizing states are attached as attribute `argmax`, a
#' data frame with one row per cell of the result.
#'
#' @param f A `pgm_factor`.
#' @param variables Variables to remove.
#' @param argmax Logical; record maximizing states.
#' @return A `pgm_factor`.
#' @examples
#' f <- discrete_factor(c("A", "B"), list(A = c("0", "1"), B = c("0", "1")), c(0.1, 0.2, 0.3, 0.4))
#' factor_marginalize(f, "B")
#' factor_maximize(f, "B")
#' @export
factor_marginalize <- function(f, variables) {
  bad <- setdiff(variables, f$vars)
  if (length(bad)) pgm_abort(paste0("Not in the factor scope: ", paste(bad, collapse = ", "), "."))
  if (!length(variables)) return(as_plain_factor(f))
  backend_marginalize(f, variables)
}

marginalize_r <- function(f, variables) {
  r <- reduce_matrix(f, variables, "sum")
  vals <- if (f$log) row_log_sum_exp(r$m) else rowSums(r$m)
  new_factor(r$keep, f$states[r$keep], vals, f$log)
}

#' @rdname factor_marginalize
#' @export
factor_maximize <- function(f, variables, argmax = FALSE) {
  bad <- setdiff(variables, f$vars)
  if (length(bad)) pgm_abort(paste0("Not in the factor scope: ", paste(bad, collapse = ", "), "."))
  if (!length(variables)) return(as_plain_factor(f))
  r <- reduce_matrix(f, variables, "max")
  which_max <- max.col(r$m, ties.method = "first")
  vals <- r$m[cbind(seq_len(nrow(r$m)), which_max)]
  out <- new_factor(r$keep, f$states[r$keep], vals, f$log)
  if (argmax) {
    dv <- intersect(r$g$vars, variables)
    dstates <- f$states[dv]
    grid <- expand.grid(dstates, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    attr(out, "argmax") <- grid[which_max, , drop = FALSE]
    rownames(attr(out, "argmax")) <- NULL
  }
  out
}

#' Reduce a factor by evidence
#'
#' Conditions on hard evidence by selecting the slice of the factor at the
#' observed states and dropping those variables. Evidence on variables outside
#' the scope is ignored.
#'
#' @param f A `pgm_factor`.
#' @param evidence Named list of observed states.
#' @return A `pgm_factor`.
#' @examples
#' f <- discrete_factor(c("A", "B"), list(A = c("0", "1"), B = c("0", "1")), c(0.1, 0.2, 0.3, 0.4))
#' factor_reduce(f, list(B = "1"))
#' @export
factor_reduce <- function(f, evidence) {
  evidence <- normalize_evidence(evidence)
  ev <- evidence[names(evidence) %in% f$vars]
  if (!length(ev)) return(as_plain_factor(f))
  idx <- lapply(f$vars, function(v) {
    if (v %in% names(ev)) {
      s <- as.character(ev[[v]])
      k <- match(s, f$states[[v]])
      if (is.na(k)) pgm_abort(sprintf("'%s' is not a state of '%s' (states: %s).", s, v,
                                      paste(f$states[[v]], collapse = ", ")), class = "pgm_unknown_state")
      k
    } else seq_len(f$card[[v]])
  })
  arr <- array(f$values, dim = unname(f$card))
  vals <- as.numeric(do.call(`[`, c(list(arr), idx, drop = FALSE)))
  keep <- setdiff(f$vars, names(ev))
  new_factor(keep, f$states[keep], vals, f$log)
}

#' Normalize a factor
#'
#' Scales the factor to sum to one. A factor whose values are all zero cannot
#' be normalized; this signals impossible evidence (an error of class
#' `pgm_impossible_evidence`).
#'
#' @param f A `pgm_factor`.
#' @return A normalized `pgm_factor` on the same scale as the input.
#' @export
factor_normalize <- function(f) {
  if (f$log) {
    z <- log_sum_exp(f$values)
    if (z == -Inf) impossible_evidence()
    f$values <- f$values - z
  } else {
    z <- sum(f$values)
    if (!(z > 0)) impossible_evidence()
    f$values <- f$values / z
  }
  as_plain_factor(f)
}

#' Divide factors
#'
#' Divides `numerator` by `denominator`, whose scope must be contained in that
#' of the numerator. Uses the convention 0/0 = 0 required by junction tree
#' message updates; a nonzero value divided by zero is an error.
#'
#' @param numerator,denominator `pgm_factor` objects.
#' @return A `pgm_factor` over the numerator's scope.
#' @export
factor_divide <- function(numerator, denominator) {
  if (!all(denominator$vars %in% numerator$vars)) {
    pgm_abort("The scope of the denominator must be contained in the scope of the numerator.")
  }
  check_compatible(numerator, denominator)
  lg <- numerator$log || denominator$log
  a <- if (lg) factor_to_log(numerator) else numerator
  b <- if (lg) factor_to_log(denominator) else denominator
  vb <- b$values[expand_index(a$card, a$vars, b$vars, b$card) + 1]
  va <- a$values
  if (lg) {
    zero_b <- vb == -Inf
    if (any(zero_b & va > -Inf)) pgm_abort("Division of a nonzero value by zero.", class = "pgm_division_by_zero")
    vals <- va - vb
    vals[zero_b] <- -Inf
  } else {
    zero_b <- vb == 0
    if (any(zero_b & va != 0)) pgm_abort("Division of a nonzero value by zero.", class = "pgm_division_by_zero")
    vals <- va / vb
    vals[zero_b] <- 0
  }
  new_factor(a$vars, a$states, vals, lg)
}

# Align values of `f` onto the scope and order of `ref` (same variable set).
aligned_values <- function(f, ref) {
  if (!setequal(f$vars, ref$vars)) pgm_abort("The factors must have the same scope.")
  check_compatible(f, ref)
  factor_reorder(f, ref$vars)$values
}

#' Information measures for discrete factors
#'
#' `factor_entropy()` computes the Shannon entropy of the normalized factor;
#' `kl_divergence()` computes KL(p || q) for factors with the same scope (Inf
#' when q is zero where p is positive); `mutual_information()` computes the
#' mutual information between two sets of variables of a joint distribution
#' factor. Natural logarithms are used unless `base` is given.
#'
#' @param f,p,q `pgm_factor` objects.
#' @param x,y Disjoint character vectors of variables in the scope of `f`.
#' @param base Logarithm base.
#' @return A numeric scalar.
#' @examples
#' f <- discrete_factor(c("A", "B"), list(A = c("0", "1"), B = c("0", "1")), c(0.4, 0.1, 0.1, 0.4))
#' factor_entropy(f, base = 2)
#' mutual_information(f, "A", "B")
#' @export
factor_entropy <- function(f, base = exp(1)) {
  p <- factor_to_linear(factor_normalize(f))$values
  p <- p[p > 0]
  -sum(p * log(p)) / log(base)
}

#' @rdname factor_entropy
#' @export
kl_divergence <- function(p, q, base = exp(1)) {
  pv <- factor_to_linear(factor_normalize(p))$values
  qv <- aligned_values(factor_to_linear(factor_normalize(q)), factor_to_linear(p))
  pos <- pv > 0
  if (any(qv[pos] == 0)) return(Inf)
  sum(pv[pos] * log(pv[pos] / qv[pos])) / log(base)
}

#' @rdname factor_entropy
#' @export
mutual_information <- function(f, x, y, base = exp(1)) {
  if (length(intersect(x, y))) pgm_abort("`x` and `y` must be disjoint.")
  if (!all(c(x, y) %in% f$vars)) pgm_abort("`x` and `y` must be in the factor scope.")
  j <- factor_normalize(factor_marginalize(f, setdiff(f$vars, c(x, y))))
  hx <- factor_entropy(factor_marginalize(j, y))
  hy <- factor_entropy(factor_marginalize(j, x))
  hxy <- factor_entropy(j)
  max(0, hx + hy - hxy) / log(base)
}

#' Value of a factor at an assignment
#'
#' @param f A `pgm_factor`.
#' @param assignment Named list or character vector giving a state for every
#'   variable in the scope.
#' @return Numeric scalar (on the factor's scale).
#' @export
factor_value <- function(f, assignment) {
  assignment <- normalize_evidence(assignment, "assignment")
  if (!all(f$vars %in% names(assignment))) pgm_abort("`assignment` must give a state for every variable.")
  factor_reduce(f, assignment[f$vars])$values
}

#' Arithmetic on factors
#'
#' `*` multiplies and `/` divides factors (see [factor_product()] and
#' [factor_divide()]).
#'
#' @param e1,e2 `pgm_factor` objects.
#' @return A `pgm_factor`.
#' @name factor-arithmetic
#' @export
Ops.pgm_factor <- function(e1, e2) {
  switch(.Generic,
    "*" = factor_product(e1, e2),
    "/" = factor_divide(e1, e2),
    pgm_abort(sprintf("Operator `%s` is not defined for factors.", .Generic))
  )
}

#' @export
all.equal.pgm_factor <- function(target, current, ...) {
  if (!inherits(current, "pgm_factor")) return("`current` is not a pgm_factor")
  if (!setequal(target$vars, current$vars)) return("Scopes differ")
  for (v in target$vars) if (!identical(target$states[[v]], current$states[[v]])) return(sprintf("States of '%s' differ", v))
  cur <- factor_reorder(as_plain_factor(current), target$vars)
  tv <- if (target$log) exp(target$values) else target$values
  cv <- if (cur$log) exp(cur$values) else cur$values
  all.equal(tv, cv, ...)
}
