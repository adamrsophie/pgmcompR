# Backend selection -------------------------------------------------------------

.pgm_state <- new.env(parent = emptyenv())
.pgm_state$backend <- "R"

#' Select the numerical backend
#'
#' pgmcompR implements every algorithm in R. The `"cpp"` backend replaces the
#' factor algebra kernels (index expansion for products and division, and
#' summation for marginalization) with compiled C++ code. Both backends compute
#' the same quantities; results agree to floating point rounding (the package
#' tests require a relative difference below 1e-12). Switching backends never
#' changes the user-facing API.
#'
#' @param backend `"R"` or `"cpp"`. If missing, the current backend is
#'   returned.
#' @return The previous backend, invisibly when setting.
#' @examples
#' old <- pgm_backend("cpp")
#' pgm_backend()
#' pgm_backend(old)
#' @export
pgm_backend <- function(backend) {
  if (missing(backend)) return(.pgm_state$backend)
  backend <- match.arg(backend, c("R", "cpp"))
  old <- .pgm_state$backend
  .pgm_state$backend <- backend
  invisible(old)
}

#' Run code with a temporary backend
#'
#' @param backend `"R"` or `"cpp"`.
#' @param code Code to evaluate.
#' @return The value of `code`.
#' @export
with_backend <- function(backend, code) {
  old <- pgm_backend(backend)
  on.exit(pgm_backend(old), add = TRUE)
  force(code)
}

backend_expand_index <- function(ucard, pos, fcard) {
  if (.pgm_state$backend == "cpp") cpp_expand_index(ucard, pos, fcard)
  else expand_index_r(ucard, pos, fcard)
}

backend_marginalize <- function(f, variables) {
  if (.pgm_state$backend == "cpp" && !f$log) {
    keep <- setdiff(f$vars, variables)
    kcard <- f$card[keep]
    outidx <- if (length(keep)) cpp_expand_index(as.integer(f$card), match(keep, f$vars) - 1L, as.integer(kcard)) else numeric(length(f$values))
    vals <- cpp_scatter_sum(f$values, outidx, prod(kcard))
    return(new_factor(keep, f$states[keep], vals, FALSE))
  }
  if (.pgm_state$backend == "cpp") {
    r <- reduce_matrix(f, variables, "sum")
    vals <- cpp_row_reduce(as.numeric(r$m), nrow(r$m), TRUE)
    return(new_factor(r$keep, f$states[r$keep], vals, TRUE))
  }
  marginalize_r(f, variables)
}
