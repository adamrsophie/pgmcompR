# Parameter estimation -------------------------------------------------------------------

data_states <- function(model, data, v) {
  st <- model$registry[[v]]$states
  if (!is.null(st)) return(st)
  x <- data[[v]]
  if (is.factor(x)) levels(x) else sort(unique(as.character(x[!is.na(x)])))
}

is_continuous_column <- function(model, data, v) {
  r <- model$registry[[v]]
  if (!is.null(r)) return(identical(r$type, "continuous"))
  is.numeric(data[[v]])
}

prepare_data <- function(model, data, allow_missing = FALSE) {
  if (!is.data.frame(data)) pgm_abort("`data` must be a data frame.")
  miss <- setdiff(nodes(model), names(data))
  if (length(miss)) pgm_abort(paste0("`data` lacks columns for: ", paste(miss, collapse = ", "), "."))
  data <- data[nodes(model)]
  if (!allow_missing && anyNA(data)) {
    pgm_abort("`data` has missing values; use method = 'em' or remove incomplete rows.", class = "pgm_missing_data")
  }
  for (v in names(data)) {
    if (!is_continuous_column(model, data, v)) {
      st <- data_states(model, data, v)
      x <- as.character(data[[v]])
      bad <- setdiff(unique(x[!is.na(x)]), st)
      if (length(bad)) pgm_abort(sprintf("Column '%s' has values outside the states {%s}: %s.", v,
                                         paste(st, collapse = ", "), paste(bad, collapse = ", ")), class = "pgm_unknown_state")
      data[[v]] <- factor(x, levels = st)
    } else if (!is.numeric(data[[v]])) {
      pgm_abort(sprintf("Column '%s' must be numeric for a continuous node.", v))
    }
  }
  data
}

# Count table (child states x parent configurations) with optional weights.
family_counts <- function(data, v, pa, states, weights = NULL) {
  vars <- c(v, pa)
  card <- vapply(states[vars], length, integer(1))
  idx <- rep(1, nrow(data))
  stride <- 1
  for (k in seq_along(vars)) {
    idx <- idx + (as.integer(data[[vars[k]]]) - 1) * stride
    stride <- stride * card[[k]]
  }
  w <- weights %||% rep(1, nrow(data))
  cnt <- numeric(prod(card))
  tab <- tapply(w, idx, sum)
  cnt[as.integer(names(tab))] <- tab
  matrix(cnt, nrow = card[[1L]])
}

#' Estimate parameters
#'
#' Fits the CPDs of a Bayesian network with known structure.
#'
#' * `"mle"`: maximum likelihood. Tabular CPDs use relative frequencies;
#'   linear Gaussian CPDs use least squares with the maximum likelihood
#'   residual variance (divisor n); conditional Gaussian CPDs fit one
#'   regression per discrete parent configuration.
#' * `"bayes"`: posterior mean under a Dirichlet prior, which is also the
#'   posterior predictive probability of the next observation. `prior` is
#'   `"bdeu"` (pseudo counts `equivalent_sample_size / (q r)` per cell, with q
#'   parent configurations and r states), `"k2"` (one per cell) or a number
#'   (the same pseudo count everywhere). `pseudo_counts` overrides the prior
#'   per variable with a number or a matrix shaped like the CPD.
#' * `"em"`: expectation maximization for data with missing values
#'   (discrete networks). Expected counts are computed exactly with a
#'   calibrated junction tree per distinct pattern of observed values. With
#'   `prior`, the M step is the posterior mode (MAP) rather than the MLE.
#'
#' Parent configurations that never occur cannot be estimated by maximum
#' likelihood. `unobserved` makes the policy explicit: `"error"` stops,
#' `"uniform"` assigns a uniform distribution and warns, `"na"` leaves `NA`
#' (the model then fails [check_model()]). Bayesian estimation always uses the
#' prior mean for such configurations.
#'
#' `structural_zeros` is a named list of logical matrices (child states by
#' parent configurations) marking impossible cells; they are fixed at zero
#' and the rest of each column is renormalized.
#'
#' @param model A `bayesian_network` (CPDs, if any, are replaced).
#' @param data A data frame with one column per node. Discrete columns may
#'   be factors or character.
#' @param method `"mle"`, `"bayes"`, `"em"` or a registered estimator.
#' @param prior Prior for `"bayes"` (and optionally `"em"`).
#' @param equivalent_sample_size BDeu equivalent sample size.
#' @param pseudo_counts Optional named list of per variable pseudo counts.
#' @param unobserved Policy for unobserved parent configurations.
#' @param structural_zeros Optional named list of logical matrices.
#' @param max_iter,tol EM controls (relative change of the log likelihood).
#' @param seed Seed for EM initialization.
#' @param ... Passed to a registered estimator.
#' @return A fitted `bayesian_network`. EM fits carry attributes
#'   `loglik` (trace), `converged` and `iterations`.
#' @examples
#' bn <- example_network("asia")
#' d <- simulate_bn(bn, 2000, seed = 1)
#' fit <- fit_parameters(bayesian_network(edges(bn)), d, method = "bayes")
#' get_cpd(fit, "lung")
#' @export
fit_parameters <- function(model, data, method = "mle", prior = "bdeu", equivalent_sample_size = 1,
                           pseudo_counts = NULL, unobserved = c("uniform", "error", "na"),
                           structural_zeros = NULL, max_iter = 100L, tol = 1e-8, seed = NULL, ...) {
  if (!inherits(model, "bayesian_network")) pgm_abort("`model` must be a bayesian_network.")
  unobserved <- match.arg(unobserved)
  est <- registry_get("estimator", method)
  est(model, data, prior = prior, equivalent_sample_size = equivalent_sample_size, pseudo_counts = pseudo_counts,
      unobserved = unobserved, structural_zeros = structural_zeros, max_iter = max_iter, tol = tol, seed = seed, ...)
}

reset_model <- function(model) {
  out <- bayesian_network(model$graph)
  declared <- Filter(function(r) isTRUE(r$declared), model$registry)
  out$registry <- declared
  out
}

estimate_mle <- function(model, data, ...) estimate_tabular(model, data, bayes = FALSE, ...)

estimate_bayes <- function(model, data, ...) estimate_tabular(model, data, bayes = TRUE, ...)

dirichlet_alpha <- function(v, card, q, prior, equivalent_sample_size, pseudo_counts) {
  pc <- pseudo_counts[[v]]
  if (!is.null(pc)) {
    if (length(pc) == 1L) return(matrix(pc, card, q))
    pc <- as.matrix(pc)
    if (!all(dim(pc) == c(card, q))) pgm_abort(sprintf("`pseudo_counts` for '%s' must be %d x %d.", v, card, q))
    return(pc)
  }
  if (is.numeric(prior)) return(matrix(prior, card, q))
  switch(prior,
    bdeu = matrix(equivalent_sample_size / (card * q), card, q),
    k2 = matrix(1, card, q),
    uniform = matrix(1, card, q),
    pgm_abort("`prior` must be 'bdeu', 'k2' or a number.")
  )
}

estimate_tabular <- function(model, data, bayes, prior = "bdeu", equivalent_sample_size = 1, pseudo_counts = NULL,
                             unobserved = "uniform", structural_zeros = NULL, weights = NULL, ...) {
  data <- prepare_data(model, data)
  out <- reset_model(model)
  states <- lapply(setNames(nodes(model), nodes(model)), function(v) if (is.factor(data[[v]])) levels(data[[v]]) else NULL)
  for (v in topological_sort(model$graph)) {
    pa <- parents(model$graph, v)
    if (is.null(states[[v]])) {
      out <- add_cpd(out, fit_continuous_cpd(v, pa, data, states, weights))
      next
    }
    if (any(vapply(states[pa], is.null, logical(1)))) {
      pgm_abort(sprintf("Discrete node '%s' has a continuous parent; this is not a conditional Gaussian network.", v))
    }
    cnt <- family_counts(data, v, pa, states, weights)
    card <- nrow(cnt); q <- ncol(cnt)
    zeros <- structural_zeros[[v]]
    if (!is.null(zeros) && !all(dim(as.matrix(zeros)) == c(card, q))) pgm_abort(sprintf("`structural_zeros` for '%s' must be %d x %d.", v, card, q))
    if (bayes) {
      alpha <- dirichlet_alpha(v, card, q, prior, equivalent_sample_size, pseudo_counts)
      num <- cnt + alpha
    } else num <- cnt
    if (!is.null(zeros)) num[as.matrix(zeros)] <- 0
    cs <- colSums(num)
    empty <- which(cs == 0)
    if (length(empty)) {
      if (unobserved == "error") {
        pgm_abort(sprintf("%d parent configuration(s) of '%s' never occur in the data; maximum likelihood is undefined. Use method = 'bayes' or `unobserved = 'uniform'`.",
                          length(empty), v), class = "pgm_unobserved_configuration")
      }
      if (unobserved == "uniform") {
        pgm_warn(sprintf("%d unobserved parent configuration(s) of '%s' were given uniform distributions.", length(empty), v),
                 class = "pgm_unobserved_configuration")
        fill <- matrix(1, card, length(empty))
        if (!is.null(zeros)) fill[as.matrix(zeros)[, empty, drop = FALSE]] <- 0
        num[, empty] <- fill
      }
    }
    probs <- sweep(num, 2L, colSums(num), "/")
    st <- states[c(v, pa)]
    cp <- new_cpd(v, pa, st, as.numeric(probs))
    attr(cp, "counts") <- cnt
    out <- add_cpd(out, cp)
  }
  out
}

fit_continuous_cpd <- function(v, pa, data, states, weights = NULL) {
  dpa <- pa[!vapply(states[pa], is.null, logical(1))]
  cpa <- setdiff(pa, dpa)
  w <- weights %||% rep(1, nrow(data))
  fit_one <- function(rows) {
    y <- data[[v]][rows]
    X <- cbind(1, as.matrix(data[rows, cpa, drop = FALSE]))
    ww <- w[rows]
    if (sum(ww > 0) <= ncol(X)) return(NULL)
    sw <- sqrt(ww)
    fit <- lm.fit(X * sw, y * sw)
    if (any(is.na(fit$coefficients))) pgm_abort(sprintf("The regression for '%s' is rank deficient (collinear parents).", v), class = "pgm_singular_matrix")
    res <- y - as.numeric(X %*% fit$coefficients)
    s <- sqrt(sum(ww * res^2) / sum(ww))
    if (!(s > 0)) pgm_abort(sprintf("'%s' is a deterministic function of its parents in the data (zero residual variance).", v), class = "pgm_singular_matrix")
    list(beta = fit$coefficients, sd = s)
  }
  if (!length(dpa)) {
    f <- fit_one(seq_len(nrow(data)))
    if (is.null(f)) pgm_abort(sprintf("Too few observations to fit '%s'.", v))
    return(lg_cpd(v, intercept = f$beta[[1L]], coefficients = setNames(f$beta[-1L], cpa), sd = f$sd))
  }
  st <- states[dpa]
  grid <- expand.grid(st, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  conf <- rep(1, nrow(data)); stride <- 1
  for (p in dpa) {
    conf <- conf + (as.integer(data[[p]]) - 1) * stride
    stride <- stride * length(st[[p]])
  }
  fits <- lapply(seq_len(nrow(grid)), function(i) fit_one(which(conf == i)))
  if (any(vapply(fits, is.null, logical(1)))) {
    pgm_abort(sprintf("Some discrete parent configurations of '%s' have too few observations for a regression.", v),
              class = "pgm_unobserved_configuration")
  }
  coef <- do.call(rbind, lapply(fits, function(f) f$beta[-1L]))
  if (is.null(dim(coef))) coef <- matrix(coef, nrow(grid), length(cpa))
  colnames(coef) <- cpa
  cg_cpd(v, dpa, st, intercept = vapply(fits, function(f) f$beta[[1L]], numeric(1)), coefficients = coef,
         sd = vapply(fits, `[[`, numeric(1), "sd"), continuous_parents = cpa)
}

# Expectation maximization -------------------------------------------------------------

estimate_em <- function(model, data, prior = NULL, equivalent_sample_size = 1, pseudo_counts = NULL,
                        unobserved = "uniform", structural_zeros = NULL, max_iter = 100L, tol = 1e-8, seed = NULL, ...) {
  data <- prepare_data(model, data, allow_missing = TRUE)
  if (any(vapply(data, is.numeric, logical(1)))) pgm_abort("EM is implemented for discrete networks.")
  states <- lapply(data, levels)
  if (!anyNA(data)) {
    return(estimate_tabular(model, data, bayes = !is.null(prior) && !identical(prior, "none"), prior = prior %||% "bdeu",
                            equivalent_sample_size = equivalent_sample_size, pseudo_counts = pseudo_counts,
                            unobserved = unobserved, structural_zeros = structural_zeros))
  }
  use_prior <- !is.null(prior) && !identical(prior, "none")
  nodes_v <- nodes(model)
  # Initialize from complete cases when possible, otherwise random CPDs.
  cur <- reset_model(model)
  for (v in nodes_v) cur$registry[[v]] <- list(type = "discrete", states = states[[v]])
  with_seed(seed, {
    for (v in topological_sort(model$graph)) {
      pa <- parents(model$graph, v)
      q <- prod(vapply(states[pa], length, integer(1)))
      m <- matrix(stats::rgamma(length(states[[v]]) * q, 5), length(states[[v]]))
      m <- sweep(m, 2L, colSums(m), "/")
      cur <- add_cpd(cur, new_cpd(v, pa, states[c(v, pa)], as.numeric(m)))
    }
  })
  key <- do.call(paste, c(lapply(data, function(x) ifelse(is.na(x), "\r", as.character(x))), sep = "\t"))
  pat <- !duplicated(key)
  upat <- data[pat, , drop = FALSE]
  wts <- as.numeric(table(factor(key, levels = key[pat])))
  ll_trace <- numeric()
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    exp_cnt <- lapply(setNames(nodes_v, nodes_v), function(v) {
      pa <- parents(model$graph, v)
      matrix(0, length(states[[v]]), prod(vapply(states[pa], length, integer(1))))
    })
    jt <- junction_tree(cur)
    ll <- 0
    for (r in seq_len(nrow(upat))) {
      row <- upat[r, , drop = FALSE]
      obs <- names(row)[!is.na(unlist(row))]
      ev <- lapply(row[obs], as.character)
      cal <- tryCatch(calibrate(jt, ev), pgm_impossible_evidence = function(e) NULL)
      if (is.null(cal)) {
        ll <- -Inf
        next
      }
      ll <- ll + wts[r] * evidence_probability(cur, ev, log = TRUE)
      for (v in nodes_v) {
        fam <- c(v, parents(model$graph, v))
        hid <- setdiff(fam, obs)
        i <- which(vapply(cal$cliques, function(cl) all(fam %in% cl), logical(1)))[1L]
        b <- cal$beliefs[[i]]
        fb <- factor_reorder(factor_marginalize(b, setdiff(b$vars, fam)), fam)
        exp_cnt[[v]] <- exp_cnt[[v]] + wts[r] * matrix(fb$values, nrow = length(states[[v]]))
        invisible(hid)
      }
    }
    ll_trace <- c(ll_trace, ll)
    new <- reset_model(model)
    for (v in nodes_v) new$registry[[v]] <- list(type = "discrete", states = states[[v]])
    for (v in topological_sort(model$graph)) {
      pa <- parents(model$graph, v)
      num <- exp_cnt[[v]]
      if (use_prior) num <- num + dirichlet_alpha(v, nrow(num), ncol(num), prior, equivalent_sample_size, pseudo_counts)
      if (!is.null(structural_zeros[[v]])) num[as.matrix(structural_zeros[[v]])] <- 0
      cs <- colSums(num)
      if (any(cs == 0)) {
        if (unobserved == "error") pgm_abort(sprintf("A parent configuration of '%s' has zero expected count.", v), class = "pgm_unobserved_configuration")
        num[, cs == 0] <- 1
      }
      new <- add_cpd(new, new_cpd(v, pa, states[c(v, pa)], as.numeric(sweep(num, 2L, colSums(num), "/"))))
    }
    cur <- new
    if (it > 1L && abs(ll_trace[it] - ll_trace[it - 1L]) < tol * (abs(ll_trace[it - 1L]) + tol)) {
      converged <- TRUE
      break
    }
  }
  if (!converged) pgm_warn(sprintf("EM did not converge in %d iterations.", max_iter), class = "pgm_convergence_warning")
  attr(cur, "loglik") <- ll_trace
  attr(cur, "converged") <- converged
  attr(cur, "iterations") <- length(ll_trace)
  cur
}

#' Log likelihood of data under a Bayesian network
#'
#' @param model A fitted `bayesian_network`.
#' @param data A complete data frame.
#' @param by_row Return per row contributions.
#' @return A numeric scalar (or vector).
#' @examples
#' bn <- example_network("asia")
#' log_likelihood(bn, simulate_bn(bn, 100, seed = 2))
#' @export
log_likelihood <- function(model, data, by_row = FALSE) {
  check_model(model)
  data <- prepare_data(model, data)
  ll <- numeric(nrow(data))
  for (v in nodes(model)) ll <- ll + cpd_log_density(model$cpds[[v]], data, model)
  if (by_row) ll else sum(ll)
}

#' Predict unobserved variables
#'
#' For each row of `newdata`, computes the posterior of `variable` given the
#' other observed columns (missing values are marginalized) by exact
#' inference.
#'
#' @param model A discrete `bayesian_network`.
#' @param newdata A data frame of evidence.
#' @param variable The variable to predict.
#' @param type `"class"` (posterior mode) or `"prob"` (posterior
#'   probabilities).
#' @param method Exact inference method.
#' @return A factor (`type = "class"`) or a matrix of probabilities.
#' @examples
#' bn <- example_network("asia")
#' nd <- data.frame(smoke = c("yes", "no"), xray = c("yes", "no"))
#' predict_bn(bn, nd, "lung", type = "prob")
#' @export
predict_bn <- function(model, newdata, variable, type = c("class", "prob"), method = "ve") {
  type <- match.arg(type)
  check_string(variable)
  st <- model$registry[[variable]]$states
  if (is.null(st)) pgm_abort("`predict_bn()` predicts discrete variables.")
  cols <- setdiff(intersect(names(newdata), nodes(model)), variable)
  key <- do.call(paste, c(lapply(newdata[cols], function(x) ifelse(is.na(x), "\r", as.character(x))), sep = "\t"))
  uk <- unique(key)
  probs <- matrix(NA_real_, length(uk), length(st), dimnames = list(NULL, st))
  for (i in seq_along(uk)) {
    r <- which(key == uk[i])[1L]
    ev <- lapply(newdata[r, cols, drop = FALSE], as.character)
    ev <- ev[!vapply(ev, is.na, logical(1))]
    probs[i, ] <- factor_values(query(model, variable, ev, method = method))
  }
  probs <- probs[match(key, uk), , drop = FALSE]
  if (type == "prob") return(probs)
  factor(st[max.col(probs, ties.method = "first")], levels = st)
}
