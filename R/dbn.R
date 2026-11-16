# Dynamic Bayesian networks -------------------------------------------------------------

prev_name <- function(v) paste0(v, "_prev")

#' Dynamic Bayesian networks
#'
#' A discrete two slice temporal Bayesian network (2TBN). `intra` arcs join
#' variables within a slice; `inter` arcs `c("A", "B")` join A at time t - 1
#' to B at time t (added later with `add_temporal_edge()`). The initial
#' slice (t = 1) uses the intra slice graph; later slices use the transition
#' graph, whose previous slice copies are named `"<variable>_prev"`.
#' `add_dbn_cpd()` attaches a CPD to the initial or transition model; with
#' `time = t` it sets a time specific transition CPD (time varying model).
#'
#' Exact inference ([filter_dbn()], [smooth_dbn()], [viterbi_dbn()],
#' [forecast_dbn()]) works on the joint state of one slice, so it is limited
#' to slices whose joint state space has at most `max_states` configurations
#' (default 4096). [particle_filter()] is an approximate alternative.
#'
#' @param states Named list of states of the slice variables.
#' @param intra,inter Arcs within and between slices.
#' @param dbn A `dynamic_bn`.
#' @param from,to Variable names for a temporal arc from `from` at t - 1 to
#'   `to` at t.
#' @param variable Variable name.
#' @param probabilities Probabilities passed to [cpd()], with parents in the
#'   order of the relevant graph (see `dbn_parents()`).
#' @param slice `"initial"` or `"transition"`.
#' @param time Optional time index (2 or larger) for a time specific
#'   transition CPD.
#' @return A `dynamic_bn`; `dbn_parents()` returns a character vector.
#' @examples
#' d <- hmm(hidden_states = c("rain", "dry"), observation_states = c("umbrella", "none"),
#'          initial = c(0.5, 0.5), transition = matrix(c(0.7, 0.3, 0.3, 0.7), 2),
#'          emission = matrix(c(0.9, 0.1, 0.2, 0.8), 2))
#' ev <- data.frame(Observation = c("umbrella", "umbrella", "none"))
#' filter_dbn(d, ev)$marginals[[3]]$Hidden
#' @export
dynamic_bn <- function(states, intra = NULL, inter = NULL) {
  if (!is.list(states) || is.null(names(states))) pgm_abort("`states` must be a named list.")
  vars <- names(states)
  states <- normalize_states(vars, states)
  g0 <- dag(intra, nodes = vars)
  ie <- parse_edges(inter)
  check_nodes_in(g0, c(ie$from, ie$to), "inter")
  dbn <- structure(list(states = states, intra = g0, inter = ie[, c("from", "to")], initial = list(),
                        transition = list(), time_cpds = list()), class = c("dynamic_bn", "pgm_model"))
  dbn
}

transition_graph <- function(dbn) {
  vars <- names(dbn$states)
  e <- rbind(dbn$intra$edges[, c("from", "to")],
             data.frame(from = prev_name(dbn$inter$from), to = dbn$inter$to, stringsAsFactors = FALSE))
  dag(e, nodes = c(prev_name(vars), vars))
}

#' @rdname dynamic_bn
#' @export
add_temporal_edge <- function(dbn, from, to) {
  check_nodes_in(dbn$intra, c(from, to))
  dbn$inter <- rbind(dbn$inter, data.frame(from = from, to = to, stringsAsFactors = FALSE))
  dbn$transition <- dbn$transition[setdiff(names(dbn$transition), to)]
  dbn
}

#' @rdname dynamic_bn
#' @export
dbn_parents <- function(dbn, variable, slice = c("initial", "transition")) {
  slice <- match.arg(slice)
  if (slice == "initial") parents(dbn$intra, variable) else parents(transition_graph(dbn), variable)
}

#' @rdname dynamic_bn
#' @export
add_dbn_cpd <- function(dbn, variable, probabilities, slice = c("initial", "transition"), time = NULL) {
  slice <- match.arg(slice)
  check_string(variable)
  vars <- names(dbn$states)
  if (!variable %in% vars) pgm_abort(sprintf("Unknown variable '%s'.", variable))
  pa <- dbn_parents(dbn, variable, slice)
  st <- dbn$states
  for (v in vars) st[[prev_name(v)]] <- dbn$states[[v]]
  x <- if (inherits(probabilities, "pgm_cpd")) probabilities else cpd(variable, probabilities, parents = pa, states = st[c(variable, pa)])
  if (!setequal(cpd_parents(x), pa)) pgm_abort(sprintf("The CPD of '%s' must have parents {%s}.", variable, paste(pa, collapse = ", ")))
  if (slice == "initial") dbn$initial[[variable]] <- x
  else if (is.null(time)) dbn$transition[[variable]] <- x
  else {
    key <- as.character(check_count(time, "time", min = 2))
    dbn$time_cpds[[key]][[variable]] <- x
  }
  dbn
}

#' @export
print.dynamic_bn <- function(x, ...) {
  cat(sprintf("<Dynamic Bayesian network: %d variables per slice, %d intra and %d temporal arcs>\n",
              length(x$states), nrow(x$intra$edges), nrow(x$inter)))
  cat(sprintf("  initial CPDs: %d of %d; transition CPDs: %d of %d%s\n", length(x$initial), length(x$states),
              length(x$transition), length(x$states), if (length(x$time_cpds)) "; time varying" else ""))
  invisible(x)
}

check_dbn <- function(dbn) {
  vars <- names(dbn$states)
  miss <- c(setdiff(vars, names(dbn$initial)), setdiff(vars, names(dbn$transition)))
  if (length(miss)) pgm_abort(paste0("The DBN lacks CPDs for: ", paste(unique(miss), collapse = ", "), "."), class = "pgm_invalid_model")
  invisible(TRUE)
}

#' @export
check_model.dynamic_bn <- function(model, error = TRUE) {
  vars <- names(model$states)
  miss <- unique(c(setdiff(vars, names(model$initial)), setdiff(vars, names(model$transition))))
  problems <- if (length(miss)) paste0("No CPD for: ", paste(miss, collapse = ", "), ".") else character()
  if (length(problems) && error) pgm_abort(paste(problems, collapse = "\n"), class = "pgm_invalid_model")
  if (length(problems)) problems else TRUE
}

# Joint state machinery -------------------------------------------------------------

dbn_initial_vector <- function(dbn) {
  vars <- names(dbn$states)
  f <- factor_reorder(factor_product(unname(dbn$initial)), vars)
  f$values
}

dbn_transition_matrix <- function(dbn, t = NULL) {
  vars <- names(dbn$states)
  cp <- dbn$transition
  if (!is.null(t) && !is.null(dbn$time_cpds[[as.character(t)]])) cp <- utils::modifyList(cp, dbn$time_cpds[[as.character(t)]])
  f <- factor_product(unname(cp))
  pv <- prev_name(vars)
  # Previous slice variables without children are absent: add them uniformly.
  missing_prev <- setdiff(pv, f$vars)
  if (length(missing_prev)) {
    st <- setNames(dbn$states[sub("_prev$", "", missing_prev)], missing_prev)
    f <- factor_product(f, new_factor(missing_prev, st, rep(1, prod(vapply(st, length, integer(1))))))
  }
  f <- factor_reorder(f, c(pv, vars))
  N <- prod(vapply(dbn$states, length, integer(1)))
  matrix(f$values, N, N)
}

dbn_state_grid <- function(dbn) expand.grid(dbn$states, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

evidence_rows <- function(evidence, vars) {
  if (is.data.frame(evidence)) {
    lapply(seq_len(nrow(evidence)), function(i) {
      r <- lapply(evidence[i, intersect(names(evidence), vars), drop = FALSE], function(x) as.character(x))
      r[!vapply(r, is.na, logical(1))]
    })
  } else if (is.list(evidence)) lapply(evidence, function(e) normalize_evidence(e))
  else pgm_abort("`evidence` must be a data frame (one row per time) or a list of named lists.")
}

evidence_likelihood <- function(grid, ev) {
  L <- rep(1, nrow(grid))
  for (v in names(ev)) {
    if (!v %in% names(grid)) pgm_abort(sprintf("Unknown variable '%s' in evidence.", v))
    L <- L * (grid[[v]] == ev[[v]])
  }
  L
}

dbn_setup <- function(dbn, evidence, max_states) {
  check_dbn(dbn)
  N <- prod(vapply(dbn$states, length, integer(1)))
  if (N > max_states) pgm_abort(sprintf("The slice has %d joint states (limit %d); use particle_filter().", N, max_states))
  grid <- dbn_state_grid(dbn)
  ev <- evidence_rows(evidence, names(dbn$states))
  list(N = N, grid = grid, ev = ev, p0 = dbn_initial_vector(dbn))
}

joint_marginals <- function(dbn, grid, p) {
  lapply(setNames(names(dbn$states), names(dbn$states)), function(v) {
    vals <- vapply(dbn$states[[v]], function(s) sum(p[grid[[v]] == s]), numeric(1))
    new_factor(v, dbn$states[v], vals)
  })
}

#' Inference in dynamic Bayesian networks
#'
#' Exact inference on the joint slice state (the forward backward algorithm
#' of the equivalent hidden Markov model). `filter_dbn()` returns filtered
#' marginals P(X_t | e_1..t) and the log likelihood of the evidence;
#' `smooth_dbn()` returns smoothed marginals P(X_t | e_1..T);
#' `viterbi_dbn()` returns the most probable joint state sequence;
#' `forecast_dbn()` returns predictive marginals for `horizon` steps after
#' the last evidence. Evidence is a data frame with one row per time (missing
#' values allowed) or a list of named lists.
#'
#' @param dbn A `dynamic_bn` with all CPDs.
#' @param evidence Evidence per time step.
#' @param horizon Number of steps to forecast.
#' @param max_states Limit on the joint slice state space.
#' @return `filter_dbn()` and `smooth_dbn()` return a list with `marginals`
#'   (a list over time of named lists of factors), `joint` (matrix of joint
#'   state probabilities, time by state) and `loglik`. `viterbi_dbn()`
#'   returns a data frame of states with attribute `log_probability`.
#'   `forecast_dbn()` returns a list over the forecast steps.
#' @export
filter_dbn <- function(dbn, evidence, max_states = 4096L) {
  s <- dbn_setup(dbn, evidence, max_states)
  Tn <- length(s$ev)
  alpha <- matrix(0, Tn, s$N)
  ll <- 0
  for (t in seq_len(Tn)) {
    prior <- if (t == 1L) s$p0 else as.numeric(alpha[t - 1L, ] %*% dbn_transition_matrix(dbn, t))
    a <- prior * evidence_likelihood(s$grid, s$ev[[t]])
    z <- sum(a)
    if (!(z > 0)) impossible_evidence(s$ev[[t]])
    ll <- ll + log(z)
    alpha[t, ] <- a / z
  }
  list(marginals = lapply(seq_len(Tn), function(t) joint_marginals(dbn, s$grid, alpha[t, ])),
       joint = alpha, loglik = ll, states = s$grid)
}

#' @rdname filter_dbn
#' @export
smooth_dbn <- function(dbn, evidence, max_states = 4096L) {
  s <- dbn_setup(dbn, evidence, max_states)
  f <- filter_dbn(dbn, evidence, max_states)
  Tn <- length(s$ev)
  beta <- matrix(1, Tn, s$N)
  gam <- f$joint
  if (Tn > 1L) for (t in (Tn - 1L):1L) {
    Tm <- dbn_transition_matrix(dbn, t + 1L)
    b <- as.numeric(Tm %*% (evidence_likelihood(s$grid, s$ev[[t + 1L]]) * beta[t + 1L, ]))
    beta[t, ] <- b / sum(b)
    g <- f$joint[t, ] * beta[t, ]
    gam[t, ] <- g / sum(g)
  }
  list(marginals = lapply(seq_len(Tn), function(t) joint_marginals(dbn, s$grid, gam[t, ])),
       joint = gam, loglik = f$loglik, states = s$grid)
}

#' @rdname filter_dbn
#' @export
viterbi_dbn <- function(dbn, evidence, max_states = 4096L) {
  s <- dbn_setup(dbn, evidence, max_states)
  Tn <- length(s$ev)
  delta <- matrix(-Inf, Tn, s$N)
  psi <- matrix(0L, Tn, s$N)
  delta[1, ] <- log(s$p0) + log(evidence_likelihood(s$grid, s$ev[[1]]))
  for (t in seq_len(Tn)[-1]) {
    lT <- log(dbn_transition_matrix(dbn, t))
    cand <- delta[t - 1L, ] + lT
    psi[t, ] <- apply(cand, 2L, which.max)
    delta[t, ] <- cand[cbind(psi[t, ], seq_len(s$N))] + log(evidence_likelihood(s$grid, s$ev[[t]]))
  }
  if (all(delta[Tn, ] == -Inf)) impossible_evidence()
  path <- integer(Tn)
  path[Tn] <- which.max(delta[Tn, ])
  if (Tn > 1L) for (t in (Tn - 1L):1L) path[t] <- psi[t + 1L, path[t + 1L]]
  out <- s$grid[path, , drop = FALSE]
  rownames(out) <- NULL
  out$time <- seq_len(Tn)
  attr(out, "log_probability") <- max(delta[Tn, ]) - filter_dbn(dbn, evidence, max_states)$loglik
  out
}

#' @rdname filter_dbn
#' @export
forecast_dbn <- function(dbn, evidence, horizon = 1L, max_states = 4096L) {
  f <- filter_dbn(dbn, evidence, max_states)
  Tn <- nrow(f$joint)
  p <- f$joint[Tn, ]
  out <- list()
  for (h in seq_len(horizon)) {
    p <- as.numeric(p %*% dbn_transition_matrix(dbn, Tn + h))
    out[[h]] <- joint_marginals(dbn, f$states, p)
  }
  out
}

#' Simulate and fit dynamic Bayesian networks
#'
#' `simulate_dbn()` draws one sequence of length `n_steps`.
#' `fit_dbn()` estimates initial and time homogeneous transition CPDs from
#' complete sequences by maximum likelihood or Dirichlet posterior means.
#'
#' @param dbn A `dynamic_bn`.
#' @param n_steps Sequence length.
#' @param seed Optional seed.
#' @param data A data frame of sequences, or a list of data frames.
#' @param id Optional name of the column identifying sequences in `data`.
#' @param method `"mle"` or `"bayes"`.
#' @param equivalent_sample_size BDeu equivalent sample size for `"bayes"`.
#' @return `simulate_dbn()` returns a data frame with a `time` column;
#'   `fit_dbn()` a fitted `dynamic_bn`.
#' @export
simulate_dbn <- function(dbn, n_steps, seed = NULL) {
  check_dbn(dbn)
  vars <- names(dbn$states)
  with_seed(seed, {
    rows <- list()
    bn0 <- bayesian_network(dbn$intra, states = dbn$states)
    for (v in vars) bn0 <- add_cpd(bn0, dbn$initial[[v]])
    cur <- simulate_bn(bn0, 1)
    rows[[1]] <- cur
    tg <- transition_graph(dbn)
    ord <- setdiff(topological_sort(tg), prev_name(vars))
    for (t in seq_len(n_steps)[-1]) {
      cp <- dbn$transition
      if (!is.null(dbn$time_cpds[[as.character(t)]])) cp <- utils::modifyList(cp, dbn$time_cpds[[as.character(t)]])
      vals <- setNames(lapply(vars, function(v) as.character(cur[[v]])), prev_name(vars))
      for (v in ord) {
        x <- cp[[v]]
        pd <- as.data.frame(vals[cpd_parents(x)], stringsAsFactors = FALSE, optional = TRUE)
        if (!length(cpd_parents(x))) pd <- data.frame(row.names = 1)
        vals[[v]] <- sample_cpd(x, pd, 1, NULL)
      }
      cur <- as.data.frame(vals[vars], stringsAsFactors = FALSE)
      rows[[t]] <- cur
    }
    out <- do.call(rbind, lapply(rows, function(r) { r[] <- lapply(r, as.character); r }))
    for (v in vars) out[[v]] <- factor(out[[v]], levels = dbn$states[[v]])
    out$time <- seq_len(n_steps)
    rownames(out) <- NULL
    out
  })
}

#' @rdname simulate_dbn
#' @export
fit_dbn <- function(dbn, data, id = NULL, method = c("mle", "bayes"), equivalent_sample_size = 1) {
  method <- match.arg(method)
  vars <- names(dbn$states)
  seqs <- if (is.data.frame(data)) {
    if (is.null(id)) list(data) else split(data, data[[id]])
  } else data
  first <- do.call(rbind, lapply(seqs, function(s) s[1, vars, drop = FALSE]))
  pairs <- do.call(rbind, lapply(seqs, function(s) {
    if (nrow(s) < 2) return(NULL)
    cur <- s[-1, vars, drop = FALSE]
    prev <- s[-nrow(s), vars, drop = FALSE]
    names(prev) <- prev_name(vars)
    cbind(prev, cur)
  }))
  bn0 <- bayesian_network(dbn$intra, states = dbn$states)
  f0 <- fit_parameters(bn0, first, method = method, equivalent_sample_size = equivalent_sample_size)
  st <- dbn$states
  for (v in vars) st[[prev_name(v)]] <- dbn$states[[v]]
  tg <- transition_graph(dbn)
  bnt <- bayesian_network(tg, states = st)
  sub <- bnt
  ft <- fit_parameters(sub, pairs[tg$nodes], method = method, equivalent_sample_size = equivalent_sample_size)
  for (v in vars) {
    dbn$initial[[v]] <- f0$cpds[[v]]
    dbn$transition[[v]] <- ft$cpds[[v]]
  }
  dbn
}

#' Particle filtering
#'
#' Approximate filtering for a discrete DBN by sequential importance
#' resampling: particles are propagated through the transition CPDs,
#' weighted by the likelihood of the evidence at each step and resampled
#' (systematic resampling) when the effective sample size falls below
#' `resample_threshold * n_particles`.
#'
#' @param dbn A `dynamic_bn`.
#' @param evidence Evidence per time step (as for [filter_dbn()]).
#' @param n_particles Number of particles.
#' @param resample_threshold ESS fraction triggering resampling.
#' @param seed Optional seed.
#' @return A list with `marginals` (as for [filter_dbn()]), `ess` per step,
#'   `loglik` (estimate) and `approximate = TRUE`.
#' @export
particle_filter <- function(dbn, evidence, n_particles = 1000L, resample_threshold = 0.5, seed = NULL) {
  check_dbn(dbn)
  vars <- names(dbn$states)
  ev <- evidence_rows(evidence, vars)
  tg <- transition_graph(dbn)
  ord <- setdiff(topological_sort(tg), prev_name(vars))
  n <- n_particles
  with_seed(seed, {
    bn0 <- bayesian_network(dbn$intra, states = dbn$states)
    for (v in vars) bn0 <- add_cpd(bn0, dbn$initial[[v]])
    P <- simulate_bn(bn0, n)
    P[] <- lapply(P, as.character)
    logw <- numeric(n)
    marg <- list(); ess <- numeric(); ll <- 0
    for (t in seq_along(ev)) {
      if (t > 1L) {
        vals <- setNames(lapply(vars, function(v) P[[v]]), prev_name(vars))
        cp <- dbn$transition
        if (!is.null(dbn$time_cpds[[as.character(t)]])) cp <- utils::modifyList(cp, dbn$time_cpds[[as.character(t)]])
        for (v in ord) {
          x <- cp[[v]]
          pd <- as.data.frame(vals[cpd_parents(x)], stringsAsFactors = FALSE, optional = TRUE)
          if (!length(cpd_parents(x))) pd <- data.frame(row.names = seq_len(n))
          vals[[v]] <- sample_cpd(x, pd, n, NULL)
        }
        P <- as.data.frame(vals[vars], stringsAsFactors = FALSE)
      }
      lik <- rep(1, n)
      for (v in names(ev[[t]])) lik <- lik * (P[[v]] == ev[[t]][[v]])
      w_prev <- exp(logw - max(logw)); w_prev <- w_prev / sum(w_prev)
      inc <- sum(w_prev * lik)
      if (!(inc > 0)) pgm_abort(sprintf("All particles have zero weight at time %d.", t), class = "pgm_sampling_failure")
      ll <- ll + log(inc)
      logw <- logw + log(lik)
      w <- exp(logw - max(logw)); w <- w / sum(w)
      ess[t] <- 1 / sum(w^2)
      marg[[t]] <- lapply(setNames(vars, vars), function(v) {
        vals <- vapply(dbn$states[[v]], function(s) sum(w[P[[v]] == s]), numeric(1))
        new_factor(v, dbn$states[v], vals)
      })
      if (ess[t] < resample_threshold * n) {
        u <- (runif(1) + 0:(n - 1)) / n
        idx <- findInterval(u, cumsum(w), left.open = TRUE) + 1L
        idx[idx > n] <- n
        P <- P[idx, , drop = FALSE]
        logw <- numeric(n)
      }
    }
    list(marginals = marg, ess = ess, loglik = ll, approximate = TRUE)
  })
}

#' Hidden Markov models
#'
#' Builds a hidden Markov model as a special case of a dynamic Bayesian
#' network with variables `Hidden` (intra arc to `Observation`, temporal arc
#' `Hidden -> Hidden`). Matrices are column stochastic: column j of
#' `transition` is P(Hidden_t | Hidden_t-1 = j), column j of `emission` is
#' P(Observation | Hidden = j).
#'
#' @param hidden_states,observation_states State labels.
#' @param initial Initial distribution of the hidden state.
#' @param transition Transition matrix (columns sum to one).
#' @param emission Emission matrix (columns sum to one).
#' @return A `dynamic_bn`.
#' @export
hmm <- function(hidden_states, observation_states, initial, transition, emission) {
  d <- dynamic_bn(list(Hidden = hidden_states, Observation = observation_states),
                  intra = list(c("Hidden", "Observation")), inter = list(c("Hidden", "Hidden")))
  d <- add_dbn_cpd(d, "Hidden", stats::setNames(initial, hidden_states), "initial")
  d <- add_dbn_cpd(d, "Observation", emission, "initial")
  d <- add_dbn_cpd(d, "Hidden", transition, "transition")
  d <- add_dbn_cpd(d, "Observation", emission, "transition")
  d
}

# Kalman filtering --------------------------------------------------------------------

#' Kalman filter and smoother
#'
#' Exact filtering and Rauch, Tung and Striebel smoothing for the linear
#' Gaussian state space model x_t = A x_t-1 + w_t, w_t ~ N(0, Q), and
#' y_t = C x_t + v_t, v_t ~ N(0, R), with x_1 ~ N(m0, P0). Rows of `y` with
#' `NA` are treated as missing (prediction only). Innovations covariances are
#' factored by Cholesky decomposition.
#'
#' @param y Observation matrix (time by observation dimension) or vector.
#' @param A,C,Q,R State space matrices.
#' @param m0,P0 Initial mean and covariance.
#' @param filtered Output of `kalman_filter()`.
#' @return `kalman_filter()` returns a list with filtered means `m`
#'   (time by state), covariances `P` (list), one step predictions `m_pred`,
#'   `P_pred` and the log likelihood `loglik`. `kalman_smoother()` adds
#'   smoothed `m_smooth` and `P_smooth`.
#' @examples
#' set.seed(1)
#' x <- cumsum(rnorm(50)); y <- x + rnorm(50, sd = 0.5)
#' kf <- kalman_filter(y, A = 1, C = 1, Q = 1, R = 0.25, m0 = 0, P0 = 10)
#' ks <- kalman_smoother(kf)
#' head(ks$m_smooth)
#' @export
kalman_filter <- function(y, A, C, Q, R, m0, P0) {
  y <- as.matrix(y)
  A <- as.matrix(A); C <- as.matrix(C); Q <- as.matrix(Q); R <- as.matrix(R); P0 <- as.matrix(P0)
  k <- nrow(A); Tn <- nrow(y)
  if (ncol(C) != k || nrow(C) != ncol(y)) pgm_abort("Dimensions of `C` do not match the state and observations.")
  m <- matrix(NA_real_, Tn, k); mp <- m
  P <- Pp <- vector("list", Tn)
  ll <- 0
  mprev <- as.numeric(m0); Pprev <- P0
  for (t in seq_len(Tn)) {
    if (t == 1L) { mpt <- mprev; Ppt <- Pprev } else {
      mpt <- as.numeric(A %*% mprev); Ppt <- A %*% Pprev %*% t(A) + Q
    }
    Ppt <- (Ppt + t(Ppt)) / 2
    mp[t, ] <- mpt; Pp[[t]] <- Ppt
    obs <- !is.na(y[t, ])
    if (any(obs)) {
      Ct <- C[obs, , drop = FALSE]
      S <- Ct %*% Ppt %*% t(Ct) + R[obs, obs, drop = FALSE]
      Rs <- safe_chol((S + t(S)) / 2, "innovation covariance")
      e <- y[t, obs] - as.numeric(Ct %*% mpt)
      K <- t(chol_solve(Rs, Ct %*% Ppt))
      mnew <- mpt + as.numeric(K %*% e)
      Pnew <- Ppt - K %*% Ct %*% Ppt
      z <- forwardsolve(t(Rs), e)
      ll <- ll - 0.5 * sum(z^2) - sum(log(diag(Rs))) - sum(obs) / 2 * log(2 * pi)
    } else { mnew <- mpt; Pnew <- Ppt }
    m[t, ] <- mnew; P[[t]] <- (Pnew + t(Pnew)) / 2
    mprev <- mnew; Pprev <- P[[t]]
  }
  list(m = m, P = P, m_pred = mp, P_pred = Pp, loglik = ll, A = A)
}

#' @rdname kalman_filter
#' @export
kalman_smoother <- function(filtered) {
  Tn <- nrow(filtered$m)
  A <- filtered$A
  ms <- filtered$m; Ps <- filtered$P
  if (Tn > 1L) for (t in (Tn - 1L):1L) {
    Rp <- safe_chol(filtered$P_pred[[t + 1L]], "predicted covariance")
    G <- t(chol_solve(Rp, A %*% filtered$P[[t]]))
    ms[t, ] <- filtered$m[t, ] + as.numeric(G %*% (ms[t + 1L, ] - filtered$m_pred[t + 1L, ]))
    Pt <- filtered$P[[t]] + G %*% (Ps[[t + 1L]] - filtered$P_pred[[t + 1L]]) %*% t(G)
    Ps[[t]] <- (Pt + t(Pt)) / 2
  }
  c(filtered, list(m_smooth = ms, P_smooth = Ps))
}
