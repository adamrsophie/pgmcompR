# Extension registries ------------------------------------------------------------

.pgm_registry <- new.env(parent = emptyenv())
.pgm_registry$inference <- list()
.pgm_registry$score <- list()
.pgm_registry$estimator <- list()
.pgm_registry$sampler <- list()

registry_set <- function(kind, name, fun, overwrite) {
  check_string(name)
  if (!is.function(fun)) pgm_abort("`fun` must be a function.")
  tab <- .pgm_registry[[kind]]
  if (!is.null(tab[[name]]) && !overwrite) {
    pgm_abort(sprintf("A %s method named '%s' already exists; use `overwrite = TRUE`.", kind, name))
  }
  tab[[name]] <- fun
  .pgm_registry[[kind]] <- tab
  invisible(name)
}

registry_get <- function(kind, name) {
  fun <- .pgm_registry[[kind]][[name]]
  if (is.null(fun)) {
    pgm_abort(sprintf("Unknown %s method '%s'. Available: %s.", kind, name,
                      paste(names(.pgm_registry[[kind]]), collapse = ", ")))
  }
  fun
}

#' Register extensions
#'
#' Adds user defined algorithms without modifying the package.
#'
#' * Inference methods are called by [query()] as
#'   `fun(model, variables, evidence, joint, virtual_evidence, interventions,
#'   elimination_order, log, ...)` and must return a normalized `pgm_factor`
#'   over `variables` (or a list of marginals when `joint = FALSE`).
#'   Approximate methods should add class `pgm_approximate` and a
#'   `diagnostics` attribute.
#' * Structure scores are called by the score based learners as
#'   `fun(node, parents, data, ...)` and must return a decomposable local
#'   score (larger is better).
#' * Parameter estimators are called by [fit_parameters()] as
#'   `fun(model, data, ...)` and must return a model with CPDs.
#' * Samplers are called by [simulate_bn()] style code as
#'   `fun(model, n, ...)` and must return a data frame.
#'
#' `list_extensions()` lists registered names.
#'
#' @param name Method name.
#' @param fun Implementation.
#' @param overwrite Replace an existing method.
#' @return The name, invisibly. `list_extensions()` returns a named list of
#'   character vectors.
#' @examples
#' register_inference_method("my_ve", function(model, variables, evidence, ...) {
#'   query(model, variables, evidence, method = "ve")
#' })
#' query(example_network("asia"), "lung", list(smoke = "yes"), method = "my_ve")
#' unregister_extension("inference", "my_ve")
#' @export
register_inference_method <- function(name, fun, overwrite = FALSE) registry_set("inference", name, fun, overwrite)

#' @rdname register_inference_method
#' @export
register_score <- function(name, fun, overwrite = FALSE) registry_set("score", name, fun, overwrite)

#' @rdname register_inference_method
#' @export
register_estimator <- function(name, fun, overwrite = FALSE) registry_set("estimator", name, fun, overwrite)

#' @rdname register_inference_method
#' @export
register_sampler <- function(name, fun, overwrite = FALSE) registry_set("sampler", name, fun, overwrite)

#' @rdname register_inference_method
#' @param kind One of `"inference"`, `"score"`, `"estimator"`, `"sampler"`.
#' @export
unregister_extension <- function(kind, name) {
  kind <- match.arg(kind, c("inference", "score", "estimator", "sampler"))
  if (name %in% .pgm_registry$builtin[[kind]]) pgm_abort("Built in methods cannot be removed.")
  tab <- .pgm_registry[[kind]]
  tab[[name]] <- NULL
  .pgm_registry[[kind]] <- tab
  invisible(name)
}

#' @rdname register_inference_method
#' @export
list_extensions <- function() {
  lapply(list(inference = "inference", score = "score", estimator = "estimator", sampler = "sampler"),
         function(k) names(.pgm_registry[[k]]))
}

get_inference_method <- function(name) registry_get("inference", name)

register_builtins <- function() {
  inf <- list(
    ve = infer_ve, enumeration = infer_enumeration, jt = infer_jt, bp = infer_bp,
    loopy_bp = infer_loopy_bp, rejection = infer_rejection, likelihood_weighting = infer_lw,
    importance = infer_importance, gibbs = infer_gibbs
  )
  .pgm_registry$inference <- inf
  .pgm_registry$score <- list(bic = local_score_bic, aic = local_score_aic, loglik = local_score_loglik,
                              bdeu = local_score_bdeu, k2 = local_score_k2, bic_g = local_score_bic_g)
  .pgm_registry$estimator <- list(mle = estimate_mle, bayes = estimate_bayes, em = estimate_em)
  .pgm_registry$sampler <- list(forward = function(model, n, ...) simulate_bn(model, n, ...))
  .pgm_registry$builtin <- lapply(list(inference = "inference", score = "score", estimator = "estimator", sampler = "sampler"),
                           function(k) names(.pgm_registry[[k]]))
  invisible(NULL)
}

.onLoad <- function(libname, pkgname) {
  register_builtins()
  invisible(NULL)
}
