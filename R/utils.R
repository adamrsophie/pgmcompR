#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib pgmcompR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rnorm runif setNames median var sd qnorm pnorm pchisq dnorm cor lm.fit
#' @importFrom utils combn head tail
## usethis namespace: end
NULL

# Condition helpers ----------------------------------------------------------

pgm_abort <- function(message, class = NULL, call = NULL, ...) {
  cond <- structure(
    class = c(class, "pgm_error", "error", "condition"),
    list(message = message, call = call, ...)
  )
  stop(cond)
}

pgm_warn <- function(message, class = NULL, ...) {
  cond <- structure(
    class = c(class, "pgm_warning", "warning", "condition"),
    list(message = message, call = NULL, ...)
  )
  warning(cond)
}

impossible_evidence <- function(evidence = NULL) {
  txt <- if (length(evidence)) {
    paste0(" (", paste(names(evidence), vapply(evidence, function(e) paste(format(e), collapse = "/"), ""),
                       sep = " = ", collapse = ", "), ")")
  } else ""
  pgm_abort(paste0("The evidence has probability zero under the model", txt, "."),
            class = "pgm_impossible_evidence", evidence = evidence)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

check_string <- function(x, arg = deparse(substitute(x))) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    pgm_abort(sprintf("`%s` must be a single non-empty string.", arg))
  }
  invisible(x)
}

check_names <- function(x, arg = deparse(substitute(x))) {
  if (!is.character(x) || anyNA(x) || any(!nzchar(x))) {
    pgm_abort(sprintf("`%s` must be a character vector of non-empty names.", arg))
  }
  if (anyDuplicated(x)) {
    pgm_abort(sprintf("`%s` contains duplicated names: %s.", arg,
                      paste(unique(x[duplicated(x)]), collapse = ", ")))
  }
  invisible(x)
}

check_count <- function(x, arg = deparse(substitute(x)), min = 1) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < min || x != round(x)) {
    pgm_abort(sprintf("`%s` must be a whole number >= %d.", arg, min))
  }
  invisible(as.integer(x))
}

check_prob <- function(x, arg = deparse(substitute(x))) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x < 0)) {
    pgm_abort(sprintf("`%s` must contain finite non-negative numbers.", arg))
  }
  invisible(x)
}

# Numerically stable log(sum(exp(x))).
log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

# Row-wise log-sum-exp of a matrix.
row_log_sum_exp <- function(m) {
  mx <- apply(m, 1L, max)
  mx_safe <- ifelse(is.finite(mx), mx, 0)
  out <- mx_safe + log(rowSums(exp(m - mx_safe)))
  out[mx == -Inf] <- -Inf
  out
}

with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  expr
}

# Coerce evidence to a named list of character state labels.
normalize_evidence <- function(evidence, arg = "evidence") {
  if (is.null(evidence) || length(evidence) == 0L) return(list())
  if (is.data.frame(evidence)) {
    if (nrow(evidence) != 1L) pgm_abort(sprintf("`%s` data frame must have one row.", arg))
    evidence <- lapply(evidence, function(v) if (is.factor(v)) as.character(v) else v)
  }
  if (is.atomic(evidence)) evidence <- as.list(evidence)
  if (!is.list(evidence) || is.null(names(evidence)) || any(!nzchar(names(evidence)))) {
    pgm_abort(sprintf("`%s` must be a named list.", arg))
  }
  check_names(names(evidence), arg)
  lapply(evidence, function(v) {
    if (is.factor(v)) v <- as.character(v)
    if (length(v) != 1L || is.na(v)) {
      pgm_abort(sprintf("Each element of `%s` must be a single non-missing value.", arg))
    }
    v
  })
}
