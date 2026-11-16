# Simulation and sampling based inference -------------------------------------------

#' Sample from a CPD
#'
#' Generic used by [simulate_bn()] and the sampling based inference methods.
#' Methods receive the parent values as columns of `parent_data` (factors for
#' discrete parents, numeric for continuous ones) and return one draw per row.
#'
#' @param x A CPD.
#' @param parent_data A data frame of parent values with `n` rows.
#' @param n Number of draws.
#' @param model The model, for state labels.
#' @param ... Unused.
#' @return A vector of length `n`: state labels (character) for discrete CPDs,
#'   numbers for continuous ones.
#' @export
sample_cpd <- function(x, parent_data, n, model, ...) UseMethod("sample_cpd")

#' @export
sample_cpd.default <- function(x, parent_data, n, model, ...) {
  f <- cpd_to_factor(x, model)
  sample_cpd.pgm_cpd(f, parent_data, n, model)
}

cpd_column_index <- function(x, parent_data) {
  col <- rep(1, nrow(parent_data))
  stride <- 1
  for (k in seq_along(x$parents)) {
    p <- x$parents[k]
    lv <- match(as.character(parent_data[[p]]), x$states[[p]])
    if (anyNA(lv)) pgm_abort(sprintf("Parent '%s' has values outside its states.", p))
    col <- col + (lv - 1) * stride
    stride <- stride * x$card[[k + 1L]]
  }
  col
}

#' @export
sample_cpd.pgm_cpd <- function(x, parent_data, n, model, ...) {
  m <- matrix(x$values, nrow = x$card[[1L]])
  col <- if (length(x$parents)) cpd_column_index(x, parent_data) else rep(1, n)
  draw_categorical(m, col, x$states[[1L]])
}

draw_categorical <- function(m, col, labels) {
  n <- length(col)
  cm <- apply(m, 2L, cumsum)
  if (!is.matrix(cm)) cm <- matrix(cm, nrow = nrow(m))
  u <- runif(n)
  k <- colSums(sweep(cm[, col, drop = FALSE], 2L, u, "<")) + 1L
  k[k > nrow(m)] <- nrow(m)
  labels[k]
}

# Probability of observed child values given parents, one per row.
cpd_density <- function(x, data, model) {
  if (inherits(x, "pgm_cpd")) {
    m <- matrix(x$values, nrow = x$card[[1L]])
    col <- if (length(x$parents)) cpd_column_index(x, data) else rep(1, nrow(data))
    row <- match(as.character(data[[x$child]]), x$states[[1L]])
    return(m[cbind(row, col)])
  }
  exp(cpd_log_density(x, data, model))
}

#' Log density of a CPD
#'
#' Generic giving log P(child | parents) for each row of `data`. Used by
#' likelihood weighting, scoring and custom CPD types.
#'
#' @param x A CPD.
#' @param data Data frame with the child and parent columns.
#' @param model The model.
#' @param ... Unused.
#' @return Numeric vector.
#' @export
cpd_log_density <- function(x, data, model, ...) UseMethod("cpd_log_density")

#' @export
cpd_log_density.pgm_cpd <- function(x, data, model, ...) log(cpd_density(x, data, model))

#' @export
cpd_log_density.default <- function(x, data, model, ...) {
  f <- cpd_to_factor(x, model)
  log(cpd_density(f, data, model))
}

as_output_frame <- function(cols, model) {
  out <- as.data.frame(cols, stringsAsFactors = FALSE, optional = TRUE)
  for (v in names(out)) {
    st <- model$registry[[v]]$states
    if (!is.null(st)) out[[v]] <- factor(out[[v]], levels = st)
  }
  out
}

#' Simulate from a Bayesian network
#'
#' Draws independent samples by ancestral (forward) sampling in topological
#' order. Interventions replace the CPDs of the intervened variables by point
#' masses before sampling, giving draws from the interventional distribution.
#'
#' @param model A Bayesian network with all CPDs.
#' @param n Number of samples.
#' @param interventions Named list of states (or values) set by intervention.
#' @param seed Optional seed; the global random number generator state is
#'   restored afterwards.
#' @return A data frame with one column per node (factors for discrete
#'   nodes, with levels in registry order).
#' @examples
#' d <- simulate_bn(example_network("sprinkler"), 5, seed = 1)
#' d
#' @export
simulate_bn <- function(model, n, interventions = NULL, seed = NULL) {
  if (!inherits(model, "bayesian_network")) pgm_abort("`model` must be a bayesian_network.")
  check_count(n, "n")
  check_model(model)
  interventions <- normalize_evidence(interventions, "interventions")
  if (length(interventions)) model <- intervene(model, interventions)
  with_seed(seed, {
    cols <- list()
    for (v in topological_sort(model$graph)) {
      x <- model$cpds[[v]]
      pa <- cpd_parents(x)
      pd <- if (length(pa)) as.data.frame(cols[pa], stringsAsFactors = FALSE, optional = TRUE) else data.frame(row.names = seq_len(n))
      cols[[v]] <- sample_cpd(x, pd, n, model)
    }
    as_output_frame(cols[model$graph$nodes], model)
  })
}

# Posterior estimate from (weighted) draws over the query variables.
weighted_estimate <- function(draws, weights, variables, states) {
  grid <- expand.grid(states[variables], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  key_grid <- do.call(paste, c(unname(as.list(grid)), sep = "\r"))
  key <- do.call(paste, c(unname(lapply(draws[variables], as.character)), sep = "\r"))
  cell <- match(key, key_grid)
  sw <- sum(weights)
  est <- vapply(seq_along(key_grid), function(k) sum(weights[cell == k]), numeric(1)) / sw
  # Delta method standard error of a self-normalized ratio estimator.
  mcse <- vapply(seq_along(key_grid), function(k) {
    ind <- as.numeric(cell == k)
    sqrt(sum(weights^2 * (ind - est[k])^2)) / sw
  }, numeric(1))
  list(est = est, mcse = mcse)
}

approx_factor <- function(variables, states, est, diagnostics, method, joint, evidence) {
  f <- new_factor(variables, states[variables], est)
  class(f) <- c("pgm_approximate", class(f))
  attr(f, "method") <- method
  attr(f, "evidence") <- evidence
  attr(f, "diagnostics") <- diagnostics
  if (joint || length(variables) == 1L) return(f)
  lapply(setNames(variables, variables), function(v) {
    m <- factor_marginalize(f, setdiff(variables, v))
    class(m) <- c("pgm_approximate", class(m))
    attr(m, "method") <- method
    attr(m, "diagnostics") <- diagnostics[setdiff(names(diagnostics), "mcse")]
    m
  })
}

#' Approximate inference
#'
#' Sampling based query methods for discrete models, used through [query()]
#' with `method` set to one of the names below. All accept `n` (samples) and
#' `seed`. Each returns an estimate of the posterior with class
#' `pgm_approximate` and a `diagnostics` attribute.
#'
#' * `"rejection"`: forward samples, keeping those that match the evidence.
#'   Diagnostics: accepted draws, acceptance rate, Monte Carlo standard error
#'   (MCSE) per cell. Soft evidence is not supported.
#' * `"likelihood_weighting"`: evidence variables are clamped and each draw is
#'   weighted by the likelihood of the evidence. Diagnostics: effective sample
#'   size (ESS, Kish), MCSE (delta method).
#' * `"importance"`: self normalized importance sampling with a user supplied
#'   `proposal` Bayesian network over the same variables (likelihood weighting
#'   is the special case `proposal = model`). Diagnostics as for likelihood
#'   weighting.
#' * `"gibbs"`: Gibbs sampling over the unobserved variables, conditioning
#'   each variable on its Markov blanket; works for Bayesian and Markov
#'   networks. Arguments `chains` (default 4), `burn_in` (default `n / 2`),
#'   `thin`. Diagnostics: split R-hat and bulk ESS per cell (Vehtari et al.,
#'   2021) and MCSE; a `pgm_convergence_warning` is raised when R-hat exceeds
#'   1.01. Deterministic relations can make the chain reducible; this is
#'   flagged in the diagnostics.
#' * `"loopy_bp"`: see [loopy_bp()].
#'
#' Approximate methods should be validated against exact inference on small
#' networks; the package tests do so.
#'
#' @name approximate_inference
#' @examples
#' bn <- example_network("sprinkler")
#' q <- query(bn, "Rain", list(WetGrass = "T"), method = "likelihood_weighting", n = 2000, seed = 1)
#' q
#' attr(q, "diagnostics")$ess
NULL

infer_rejection <- function(model, variables, evidence, joint = TRUE, virtual_evidence = NULL,
                            interventions = NULL, n = 5000, seed = NULL, max_draws = 1e6, ...) {
  if (!inherits(model, "bayesian_network")) pgm_abort("Rejection sampling needs a Bayesian network.")
  if (length(virtual_evidence)) pgm_abort("Rejection sampling does not support soft evidence; use likelihood weighting.")
  prep <- prepare_discrete(model, variables, evidence, NULL, interventions, prune = FALSE)
  ev <- prep$evidence
  with_seed(seed, {
    acc <- NULL
    total <- 0
    while ((is.null(acc) || nrow(acc) < n) && total < max_draws) {
      batch <- min(max(1000, 2 * n), max_draws - total)
      d <- simulate_bn(prep$model, batch)
      total <- total + batch
      keep <- rep(TRUE, batch)
      for (v in names(ev)) keep <- keep & as.character(d[[v]]) == ev[[v]]
      acc <- rbind(acc, d[keep, , drop = FALSE])
    }
    if (is.null(acc) || !nrow(acc)) {
      pgm_abort(sprintf("No sample matched the evidence in %g draws; the evidence is very unlikely or impossible. Use likelihood weighting.", total),
                class = "pgm_sampling_failure")
    }
    acc <- utils::head(acc, n)
    est <- weighted_estimate(acc, rep(1, nrow(acc)), prep$variables, prep$states)
    diag <- list(n = nrow(acc), draws = total, acceptance_rate = nrow(acc) / total, ess = nrow(acc), mcse = est$mcse)
    if (nrow(acc) < n) pgm_warn(sprintf("Only %d of %d requested samples were accepted.", nrow(acc), n), class = "pgm_sampling_warning")
    approx_factor(prep$variables, prep$states, est$est, diag, "rejection", joint, ev)
  })
}

# Likelihood weighting draws, optionally sampling from a proposal network.
lw_draws <- function(model, evidence, virtual_evidence, n, proposal = NULL) {
  q <- proposal %||% model
  ord <- topological_sort(q$graph)
  cols <- list()
  logw <- numeric(n)
  for (v in ord) {
    if (v %in% names(evidence)) {
      cols[[v]] <- rep(evidence[[v]], n)
    } else {
      xq <- q$cpds[[v]]
      pd <- as.data.frame(cols[cpd_parents(xq)], stringsAsFactors = FALSE, optional = TRUE)
      if (!length(cpd_parents(xq))) pd <- data.frame(row.names = seq_len(n))
      cols[[v]] <- sample_cpd(xq, pd, n, q)
    }
  }
  d <- as.data.frame(cols, stringsAsFactors = FALSE, optional = TRUE)
  for (v in model$graph$nodes) {
    xp <- model$cpds[[v]]
    lp <- log(cpd_density(xp, d, model))
    if (v %in% names(evidence)) {
      logw <- logw + lp
    } else if (!is.null(proposal)) {
      logw <- logw + lp - log(cpd_density(q$cpds[[v]], d, q))
    }
  }
  for (f in virtual_evidence_factors(virtual_evidence, model_states(model))) {
    v <- f$vars
    logw <- logw + log(f$values[match(d[[v]], f$states[[v]])])
  }
  list(draws = d, logw = logw)
}

lw_finish <- function(res, prep, joint, method, n) {
  lw <- res$logw
  if (all(lw == -Inf)) {
    pgm_abort("Every sample has zero weight: the evidence is impossible or the proposal does not cover the posterior.",
              class = "pgm_sampling_failure")
  }
  w <- exp(lw - max(lw))
  est <- weighted_estimate(res$draws, w, prep$variables, prep$states)
  ess <- sum(w)^2 / sum(w^2)
  diag <- list(n = n, ess = ess, mcse = est$mcse, log_mean_weight = max(lw) + log(mean(w)))
  if (ess < 0.01 * n) pgm_warn(sprintf("Effective sample size is only %.1f of %d draws.", ess, n), class = "pgm_sampling_warning")
  approx_factor(prep$variables, prep$states, est$est, diag, method, joint, prep$evidence)
}

infer_lw <- function(model, variables, evidence, joint = TRUE, virtual_evidence = NULL,
                     interventions = NULL, n = 5000, seed = NULL, ...) {
  if (!inherits(model, "bayesian_network")) pgm_abort("Likelihood weighting needs a Bayesian network.")
  prep <- prepare_discrete(model, variables, evidence, virtual_evidence, interventions, prune = FALSE)
  with_seed(seed, {
    res <- lw_draws(prep$model, prep$evidence, virtual_evidence, n)
    lw_finish(res, prep, joint, "likelihood_weighting", n)
  })
}

infer_importance <- function(model, variables, evidence, joint = TRUE, virtual_evidence = NULL,
                             interventions = NULL, n = 5000, seed = NULL, proposal = NULL, ...) {
  if (!inherits(model, "bayesian_network")) pgm_abort("Importance sampling needs a Bayesian network.")
  prep <- prepare_discrete(model, variables, evidence, virtual_evidence, interventions, prune = FALSE)
  proposal <- proposal %||% prep$model
  if (!inherits(proposal, "bayesian_network") || !setequal(nodes(proposal), nodes(model))) {
    pgm_abort("`proposal` must be a Bayesian network over the same variables.")
  }
  for (v in nodes(model)) {
    if (!identical(proposal$registry[[v]]$states, model$registry[[v]]$states)) {
      pgm_abort(sprintf("The proposal uses different states for '%s'.", v))
    }
  }
  with_seed(seed, {
    res <- lw_draws(prep$model, prep$evidence, virtual_evidence, n, proposal = proposal)
    lw_finish(res, prep, joint, "importance", n)
  })
}

# Gibbs sampling -----------------------------------------------------------------

infer_gibbs <- function(model, variables, evidence, joint = TRUE, virtual_evidence = NULL,
                        interventions = NULL, n = 5000, seed = NULL, chains = 4L, burn_in = NULL,
                        thin = 1L, ...) {
  prep <- prepare_discrete(model, variables, evidence, virtual_evidence, interventions, prune = FALSE)
  chains <- check_count(chains, "chains")
  thin <- check_count(thin, "thin")
  burn_in <- burn_in %||% max(100L, n %/% 2L)
  fs <- reduce_all(prep$factors, prep$evidence)
  fs <- fs[vapply(fs, function(f) length(f$vars) > 0L, logical(1))]
  free <- setdiff(nodes(prep$model), names(prep$evidence))
  free <- intersect(free, unique(unlist(lapply(fs, `[[`, "vars"), use.names = FALSE)))
  st <- prep$states
  card <- vapply(st[free], length, integer(1))
  deterministic <- any(vapply(fs, function(f) any(f$values == 0), logical(1)))
  with_seed(seed, {
    S <- gibbs_init(prep, free, chains)
    fac_of <- lapply(setNames(free, free), function(v) which(vapply(fs, function(f) v %in% f$vars, logical(1))))
    fidx <- lapply(fs, function(f) {
      list(pos = match(f$vars, free), stride = c(1, cumprod(f$card))[seq_along(f$vars)])
    })
    total <- burn_in + n * thin
    keep_rows <- burn_in + seq(thin, n * thin, by = thin)
    trace <- array(NA_integer_, c(n, chains, length(prep$variables)))
    qpos <- match(prep$variables, free)
    rec <- 0L
    for (t in seq_len(total)) {
      for (j in seq_along(free)) {
        v <- free[j]
        lp <- matrix(0, chains, card[[j]])
        for (fi in fac_of[[v]]) {
          f <- fs[[fi]]
          ix <- fidx[[fi]]
          base <- rep(1, chains)
          for (k in seq_along(ix$pos)) {
            if (ix$pos[k] != j) base <- base + (S[, ix$pos[k]] - 1) * ix$stride[k]
          }
          sj <- ix$stride[match(j, ix$pos)]
          for (s in seq_len(card[[j]])) lp[, s] <- lp[, s] + log(f$values[base + (s - 1) * sj])
        }
        mx <- apply(lp, 1L, max)
        p <- exp(lp - mx)
        p <- p / rowSums(p)
        cum <- p
        if (card[[j]] > 1L) for (s in 2:card[[j]]) cum[, s] <- cum[, s - 1L] + p[, s]
        u <- runif(chains)
        S[, j] <- pmin(rowSums(u > cum) + 1L, card[[j]])
      }
      if (t %in% keep_rows) {
        rec <- rec + 1L
        trace[rec, , ] <- S[, qpos, drop = FALSE]
      }
    }
    grid <- expand.grid(lapply(st[prep$variables], seq_along), KEEP.OUT.ATTRS = FALSE)
    qcard <- vapply(st[prep$variables], length, integer(1))
    qstride <- c(1, cumprod(qcard))[seq_along(qcard)]
    cellm <- matrix(1, n, chains)
    for (k in seq_along(prep$variables)) cellm <- cellm + (trace[, , k] - 1) * qstride[k]
    ncell <- prod(qcard)
    est <- vapply(seq_len(ncell), function(k) mean(cellm == k), numeric(1))
    diag_cells <- lapply(seq_len(ncell), function(k) mcmc_diagnostics(matrix(as.numeric(cellm == k), n, chains)))
    rhat <- vapply(diag_cells, `[[`, numeric(1), "rhat")
    ess <- vapply(diag_cells, `[[`, numeric(1), "ess")
    mcse <- vapply(seq_len(ncell), function(k) {
      s <- sqrt(est[k] * (1 - est[k]))
      if (s == 0) 0 else s / sqrt(max(ess[k], 1))
    }, numeric(1))
    diag <- list(n = n * chains, chains = chains, burn_in = burn_in, thin = thin, rhat = rhat, ess = ess,
                 mcse = mcse, deterministic_factors = deterministic)
    if (any(rhat > 1.01, na.rm = TRUE)) {
      pgm_warn(sprintf("Gibbs chains have not mixed: max R-hat %.3f > 1.01. Increase `n` or `burn_in`.", max(rhat, na.rm = TRUE)),
               class = "pgm_convergence_warning")
    }
    approx_factor(prep$variables, prep$states, est, diag, "gibbs", joint, prep$evidence)
  })
}

gibbs_init <- function(prep, free, chains) {
  st <- prep$states
  S <- matrix(1L, chains, length(free))
  model <- prep$model
  if (inherits(model, "bayesian_network")) {
    res <- lw_draws(model, prep$evidence, NULL, max(50L * chains, 200L))
    ok <- which(res$logw > -Inf)
    if (!length(ok)) impossible_evidence(prep$evidence)
    pick <- ok[sample.int(length(ok), chains, replace = length(ok) < chains)]
    for (j in seq_along(free)) S[, j] <- match(as.character(res$draws[[free[j]]][pick]), st[[free[j]]])
    return(S)
  }
  # Markov networks: start from the MPE assignment for every chain, then
  # randomize where the potentials allow it.
  mpe <- mpe_query(model, prep$evidence)
  for (j in seq_along(free)) S[, j] <- match(mpe$assignment[[free[j]]], st[[free[j]]])
  S
}

# Split R-hat and bulk ESS (rank normalized) for a draws x chains matrix.
mcmc_diagnostics <- function(x) {
  n <- nrow(x)
  if (n < 4L) return(list(rhat = NA_real_, ess = NA_real_))
  if (stats::var(as.numeric(x)) == 0) return(list(rhat = NA_real_, ess = length(x)))
  half <- n %/% 2L
  split <- cbind(x[seq_len(half), , drop = FALSE], x[(n - half + 1L):n, , drop = FALSE])
  r <- rank(as.numeric(split), ties.method = "average")
  z <- matrix(qnorm((r - 3 / 8) / (length(r) + 1 / 4)), nrow = half)
  list(rhat = rhat_basic(z), ess = ess_basic(z))
}

rhat_basic <- function(z) {
  m <- ncol(z); n <- nrow(z)
  cm <- colMeans(z)
  w <- mean(apply(z, 2L, stats::var))
  b <- n * stats::var(cm)
  if (w == 0) return(if (b == 0) 1 else Inf)
  sqrt(((n - 1) / n * w + b / n) / w)
}

ess_basic <- function(z) {
  m <- ncol(z); n <- nrow(z)
  cm <- colMeans(z)
  w <- mean(apply(z, 2L, stats::var))
  if (w == 0) return(m * n)
  vplus <- (n - 1) / n * w + if (m > 1) stats::var(cm) else 0
  acov <- vapply(seq_len(m), function(j) {
    y <- z[, j] - mean(z[, j])
    ac <- stats::acf(y, lag.max = n - 1L, type = "covariance", plot = FALSE, demean = FALSE)$acf[, 1, 1]
    ac
  }, numeric(n))
  if (!is.matrix(acov)) acov <- matrix(acov, ncol = m)
  rho <- 1 - (w - rowMeans(acov)) / vplus
  rho[1] <- 1
  # Geyer's initial monotone sequence on paired sums.
  s <- 0
  prev <- Inf
  k <- 1L
  while (k + 1L <= n) {
    pr <- rho[k] + rho[k + 1L]
    if (pr < 0) break
    pr <- min(pr, prev)
    s <- s + pr
    prev <- pr
    k <- k + 2L
  }
  tau <- -1 + 2 * s
  max(1, min(m * n * log10(m * n), m * n / max(tau, 1e-8)))
}

#' @export
print.pgm_approximate <- function(x, ...) {
  d <- attr(x, "diagnostics")
  cat(sprintf("<Approximate posterior (%s)>", attr(x, "method") %||% "approximate"))
  if (!is.null(d$ess)) cat(sprintf(" ESS %s", paste(format(round(range(d$ess), 1)), collapse = " to ")))
  if (!is.null(d$rhat) && any(is.finite(d$rhat))) cat(sprintf(", max R-hat %.3f", max(d$rhat, na.rm = TRUE)))
  if (!is.null(d$converged)) cat(sprintf(", %s", if (d$converged) "converged" else "NOT converged"))
  cat("\n")
  df <- as.data.frame(x)
  if (!is.null(d$mcse) && length(d$mcse) == nrow(df)) df$mcse <- d$mcse
  print(df, row.names = FALSE)
  invisible(x)
}
