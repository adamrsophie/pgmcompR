#' The CPD contract
#'
#' Any S3 class can serve as a conditional distribution in a
#' [bayesian_network()] if objects carry `child` (or `variable`) and
#' `parents` elements and the class provides methods for these generics:
#'
#' * [sample_cpd()] (required for simulation and sampling based inference);
#' * [cpd_log_density()] (required for likelihood weighting, scoring and
#'   [log_likelihood()]);
#' * [cpd_to_factor()] (required for exact discrete inference): a
#'   `pgm_factor` with the child first and every column normalized;
#' * [validate_cpd()] (optional): signals an error if the CPD is invalid for
#'   the model.
#'
#' Discrete custom CPDs should also set `registry_entries`, a named list of
#' `list(type = "discrete", states = ...)` for the child and parents, so the
#' model registry learns their states. Contracts: the factor's dimension order
#' is (child, parents in the CPD's order); states follow the registry order;
#' probabilities are finite and nonnegative and every column sums to one;
#' `cpd_log_density()` returns `-Inf` (never `NaN`) for impossible values.
#'
#' @name cpd_contract
#' @examples
#' # A logistic CPD for a binary child with one binary parent
#' logistic_cpd <- function(child, parent, a, b) {
#'   structure(list(child = child, parents = parent, a = a, b = b,
#'                  registry_entries = list()), class = "logistic_cpd")
#' }
#' cpd_to_factor.logistic_cpd <- function(x, model, ...) {
#'   p1 <- plogis(x$a + x$b * c(0, 1))
#'   cpd(x$child, rbind(1 - p1, p1), parents = x$parents,
#'       states = model$registry[c(x$child, x$parents)] |> lapply(`[[`, "states"))
#' }
#' registerS3method("cpd_to_factor", "logistic_cpd", cpd_to_factor.logistic_cpd)
#' bn <- bayesian_network(list(c("A", "B")), states = list(A = c("0", "1"), B = c("0", "1")))
#' bn <- add_cpd(bn, cpd("A", c(`0` = 0.5, `1` = 0.5)))
#' bn <- add_cpd(bn, logistic_cpd("B", "A", -1, 2))
#' query(bn, "B")
NULL
