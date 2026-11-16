# Gaussian factors, linear Gaussian and conditional Gaussian CPDs --------------------

safe_chol <- function(S, what = "covariance matrix") {
  if (!length(S)) return(matrix(0, 0, 0))
  if (any(!is.finite(S))) pgm_abort(sprintf("The %s has non finite entries.", what), class = "pgm_singular_matrix")
  if (max(abs(S - t(S))) > 1e-8 * max(1, max(abs(S)))) pgm_abort(sprintf("The %s is not symmetric.", what), class = "pgm_singular_matrix")
  R <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(R)) {
    pgm_abort(sprintf("The %s is not positive definite (it is singular or indefinite).", what), class = "pgm_singular_matrix")
  }
  R
}

chol_solve <- function(R, b) backsolve(R, forwardsolve(t(R), b))

#' Gaussian factors
#'
#' `gaussian_factor()` represents a multivariate normal density (or, more
#' generally, a Gaussian potential) in canonical form
#' exp(-x'Kx/2 + h'x + g). Construct it from a mean and covariance (moment
#' form), or from `precision`, `h` and `g` (canonical form, which may be
#' degenerate, as for linear Gaussian CPDs). Products add canonical
#' parameters; marginalization and conditioning use Cholesky factorizations
#' and never form explicit inverses of the full matrix. Non positive definite
#' blocks raise an error of class `pgm_singular_matrix`.
#'
#' @param variables Character vector of variable names.
#' @param mean,covariance Moment parameters.
#' @param precision,h,g Canonical parameters (used when `mean` is `NULL`).
#' @param f,a,b Gaussian factors.
#' @param evidence Named numeric vector or list of observed values.
#' @return A `pgm_gaussian` object; `gaussian_moments()` returns a list with
#'   `mean` and `covariance`.
#' @examples
#' f <- gaussian_factor(c("X", "Y"), mean = c(0, 1), covariance = matrix(c(1, 0.5, 0.5, 2), 2))
#' gaussian_moments(gaussian_reduce(f, c(Y = 2)))
#' gaussian_moments(gaussian_marginalize(f, "Y"))
#' @export
gaussian_factor <- function(variables, mean = NULL, covariance = NULL, precision = NULL, h = NULL, g = 0) {
  check_names(variables, "variables")
  k <- length(variables)
  if (!is.null(mean)) {
    if (length(mean) != k) pgm_abort("`mean` must have one entry per variable.")
    covariance <- as.matrix(covariance)
    if (!all(dim(covariance) == c(k, k))) pgm_abort("`covariance` must be a square matrix matching `variables`.")
    R <- safe_chol(covariance)
    K <- chol2inv(R)
    hh <- as.numeric(K %*% mean)
    g <- -0.5 * sum(mean * hh) - sum(log(diag(R))) - k / 2 * log(2 * pi)
    out <- new_gaussian(variables, K, hh, g)
    out$moments <- list(mean = setNames(as.numeric(mean), variables), covariance = unname_dims(covariance, variables))
    return(out)
  }
  precision <- as.matrix(precision %||% matrix(0, k, k))
  h <- as.numeric(h %||% numeric(k))
  if (!all(dim(precision) == c(k, k)) || length(h) != k) pgm_abort("Canonical parameters do not match `variables`.")
  new_gaussian(variables, precision, h, g)
}

unname_dims <- function(m, vars) {
  m <- as.matrix(m)
  dimnames(m) <- list(vars, vars)
  m
}

new_gaussian <- function(vars, K, h, g) {
  K <- unname_dims((K + t(K)) / 2, vars)
  structure(list(vars = vars, K = K, h = setNames(as.numeric(h), vars), g = g), class = "pgm_gaussian")
}

#' @rdname gaussian_factor
#' @export
gaussian_moments <- function(f) {
  if (!is.null(f$moments)) return(f$moments)
  R <- safe_chol(f$K, "precision matrix")
  S <- chol2inv(R)
  list(mean = setNames(as.numeric(S %*% f$h), f$vars), covariance = unname_dims(S, f$vars))
}

#' @export
print.pgm_gaussian <- function(x, ...) {
  cat(sprintf("<Gaussian factor: %s>\n", paste(x$vars, collapse = ", ")))
  m <- tryCatch(gaussian_moments(x), error = function(e) NULL)
  if (is.null(m)) {
    cat("  canonical form (precision not positive definite)\n")
    print(x$K)
  } else {
    cat("  mean:\n")
    print(round(m$mean, 6))
    cat("  covariance:\n")
    print(round(m$covariance, 6))
  }
  invisible(x)
}

expand_gaussian <- function(f, vars) {
  K <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  h <- setNames(numeric(length(vars)), vars)
  K[f$vars, f$vars] <- f$K
  h[f$vars] <- f$h
  list(K = K, h = h)
}

#' @rdname gaussian_factor
#' @export
gaussian_product <- function(a, b) {
  vars <- c(a$vars, setdiff(b$vars, a$vars))
  ea <- expand_gaussian(a, vars)
  eb <- expand_gaussian(b, vars)
  new_gaussian(vars, ea$K + eb$K, ea$h + eb$h, a$g + b$g)
}

#' @rdname gaussian_factor
#' @export
gaussian_marginalize <- function(f, variables) {
  bad <- setdiff(variables, f$vars)
  if (length(bad)) pgm_abort(paste0("Not in the factor scope: ", paste(bad, collapse = ", "), "."))
  keep <- setdiff(f$vars, variables)
  if (!is.null(f$moments)) {
    out <- gaussian_factor(keep, f$moments$mean[keep], f$moments$covariance[keep, keep, drop = FALSE])
    return(out)
  }
  Kyy <- f$K[variables, variables, drop = FALSE]
  R <- safe_chol(Kyy, "precision block being integrated out")
  Kxy <- f$K[keep, variables, drop = FALSE]
  hy <- f$h[variables]
  A <- chol_solve(R, t(Kxy))
  a <- chol_solve(R, hy)
  K <- f$K[keep, keep, drop = FALSE] - Kxy %*% A
  h <- f$h[keep] - as.numeric(Kxy %*% a)
  g <- f$g + 0.5 * (length(variables) * log(2 * pi) - 2 * sum(log(diag(R))) + sum(hy * a))
  new_gaussian(keep, K, h, g)
}

#' @rdname gaussian_factor
#' @export
gaussian_reduce <- function(f, evidence) {
  ev <- unlist(evidence)
  if (is.null(names(ev))) pgm_abort("`evidence` must be named.")
  ev <- ev[names(ev) %in% f$vars]
  if (!length(ev)) return(f)
  y <- names(ev)
  keep <- setdiff(f$vars, y)
  if (!is.null(f$moments)) {
    m <- f$moments
    if (!length(keep)) return(gaussian_factor(character(), precision = matrix(0, 0, 0), h = numeric(), g = f$g))
    R <- safe_chol(m$covariance[y, y, drop = FALSE])
    Sxy <- m$covariance[keep, y, drop = FALSE]
    d <- as.numeric(ev) - m$mean[y]
    mu <- m$mean[keep] + as.numeric(Sxy %*% chol_solve(R, d))
    S <- m$covariance[keep, keep, drop = FALSE] - Sxy %*% chol_solve(R, t(Sxy))
    return(gaussian_factor(keep, mu, (S + t(S)) / 2))
  }
  yv <- as.numeric(ev)
  K <- f$K[keep, keep, drop = FALSE]
  h <- f$h[keep] - as.numeric(f$K[keep, y, drop = FALSE] %*% yv)
  g <- f$g + sum(f$h[y] * yv) - 0.5 * as.numeric(t(yv) %*% f$K[y, y, drop = FALSE] %*% yv)
  new_gaussian(keep, K, h, g)
}

# Linear Gaussian CPDs ---------------------------------------------------------------

#' Linear Gaussian and conditional Gaussian CPDs
#'
#' `lg_cpd()` defines X | U ~ N(intercept + sum(coefficients * U), sd^2) for
#' continuous parents U. `cg_cpd()` defines a conditional Gaussian CPD for a
#' continuous child with discrete parents D and continuous parents U: one
#' linear Gaussian regression per configuration of D (first discrete parent
#' varying fastest). Discrete children may not have continuous parents
#' (conditional Gaussian networks in the sense of Lauritzen and Wermuth).
#'
#' @param variable Child name.
#' @param intercept Intercept (for `cg_cpd()`, one per discrete configuration).
#' @param coefficients Named numeric vector of regression coefficients (for
#'   `cg_cpd()`, a matrix with one row per discrete configuration and one
#'   named column per continuous parent).
#' @param sd Residual standard deviation (positive; for `cg_cpd()` one per
#'   configuration).
#' @param discrete_parents Discrete parent names.
#' @param states Named list of states of the discrete parents.
#' @return A `pgm_lg_cpd` or `pgm_cg_cpd`.
#' @examples
#' lg_cpd("Y", intercept = 1, coefficients = c(X = 2), sd = 0.5)
#' cg_cpd("Weight", discrete_parents = "Sex", states = list(Sex = c("F", "M")),
#'        intercept = c(62, 78), sd = c(8, 10))
#' @export
lg_cpd <- function(variable, intercept = 0, coefficients = NULL, sd = 1) {
  check_string(variable)
  coefficients <- coefficients %||% setNames(numeric(), character())
  if (length(coefficients) && (is.null(names(coefficients)) || any(!nzchar(names(coefficients))))) {
    pgm_abort("`coefficients` must be named by the parents.")
  }
  if (!is.numeric(sd) || length(sd) != 1L || !(sd > 0) || !is.finite(sd)) pgm_abort("`sd` must be a positive number.")
  if (!is.numeric(intercept) || length(intercept) != 1L || !is.finite(intercept)) pgm_abort("`intercept` must be a finite number.")
  structure(list(variable = variable, child = variable, parents = names(coefficients),
                 continuous_parents = names(coefficients), intercept = intercept,
                 coefficients = coefficients, sd = sd), class = "pgm_lg_cpd")
}

#' @export
print.pgm_lg_cpd <- function(x, ...) {
  rhs <- c(format(x$intercept, digits = 6),
           if (length(x$coefficients)) paste0(format(x$coefficients, digits = 6), " * ", names(x$coefficients)))
  cat(sprintf("<Linear Gaussian CPD> %s ~ N(%s, %s^2)\n", x$variable, paste(rhs, collapse = " + "), format(x$sd, digits = 6)))
  invisible(x)
}

#' @rdname lg_cpd
#' @param continuous_parents Continuous parent names.
#' @export
cg_cpd <- function(variable, discrete_parents, states, intercept, coefficients = NULL, sd,
                   continuous_parents = NULL) {
  check_string(variable)
  st <- normalize_states(discrete_parents, states)
  nconf <- prod(vapply(st, length, integer(1)))
  cp <- continuous_parents %||% colnames(coefficients) %||% character()
  coefficients <- coefficients %||% matrix(0, nconf, 0)
  coefficients <- as.matrix(coefficients)
  if (length(cp) && is.null(colnames(coefficients))) colnames(coefficients) <- cp
  if (length(intercept) != nconf || length(sd) != nconf || nrow(coefficients) != nconf) {
    pgm_abort(sprintf("`intercept`, `sd` and the rows of `coefficients` need one entry per discrete configuration (%d).", nconf))
  }
  if (any(!(sd > 0))) pgm_abort("`sd` must be positive.")
  if (!setequal(colnames(coefficients), cp)) pgm_abort("The columns of `coefficients` must name the continuous parents.")
  if (length(cp)) coefficients <- coefficients[, cp, drop = FALSE]
  structure(list(variable = variable, child = variable, parents = c(discrete_parents, cp),
                 discrete_parents = discrete_parents, continuous_parents = cp, discrete_states = st,
                 intercept = as.numeric(intercept), coefficients = coefficients, sd = as.numeric(sd)),
            class = "pgm_cg_cpd")
}

#' @export
print.pgm_cg_cpd <- function(x, ...) {
  cat(sprintf("<Conditional Gaussian CPD> %s | %s\n", x$variable, paste(x$parents, collapse = ", ")))
  grid <- expand.grid(x$discrete_states, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  tab <- cbind(grid, intercept = x$intercept, x$coefficients, sd = x$sd)
  print(tab, row.names = FALSE)
  invisible(x)
}

#' @export
validate_cpd.pgm_lg_cpd <- function(x, model, ...) {
  for (p in x$continuous_parents) {
    if (!identical(model$registry[[p]]$type, "continuous")) pgm_abort(sprintf("Parent '%s' of '%s' must be continuous.", p, x$variable))
  }
  invisible(TRUE)
}

#' @export
validate_cpd.pgm_cg_cpd <- function(x, model, ...) {
  for (p in x$discrete_parents) {
    r <- model$registry[[p]]
    if (!identical(r$type, "discrete") || !identical(r$states, x$discrete_states[[p]])) {
      pgm_abort(sprintf("Discrete parent '%s' of '%s' does not match the registry.", p, x$variable))
    }
  }
  validate_cpd.pgm_lg_cpd(x, model)
}

cg_config_index <- function(x, data) {
  idx <- rep(1, nrow(data))
  stride <- 1
  for (p in x$discrete_parents) {
    lv <- match(as.character(data[[p]]), x$discrete_states[[p]])
    idx <- idx + (lv - 1) * stride
    stride <- stride * length(x$discrete_states[[p]])
  }
  idx
}

lg_mean <- function(x, data, conf = NULL) {
  if (inherits(x, "pgm_cg_cpd")) {
    conf <- conf %||% cg_config_index(x, data)
    mu <- x$intercept[conf]
    for (p in x$continuous_parents) mu <- mu + x$coefficients[conf, p] * as.numeric(data[[p]])
    return(list(mean = mu, sd = x$sd[conf]))
  }
  mu <- rep(x$intercept, nrow(data))
  for (p in names(x$coefficients)) mu <- mu + x$coefficients[[p]] * as.numeric(data[[p]])
  list(mean = mu, sd = rep(x$sd, nrow(data)))
}

#' @export
sample_cpd.pgm_lg_cpd <- function(x, parent_data, n, model, ...) {
  if (!nrow(parent_data)) parent_data <- data.frame(row.names = seq_len(n))
  m <- lg_mean(x, parent_data)
  rnorm(n, m$mean, m$sd)
}

#' @export
sample_cpd.pgm_cg_cpd <- function(x, parent_data, n, model, ...) {
  m <- lg_mean(x, parent_data)
  rnorm(n, m$mean, m$sd)
}

#' @export
cpd_log_density.pgm_lg_cpd <- function(x, data, model, ...) {
  m <- lg_mean(x, data)
  dnorm(as.numeric(data[[x$variable]]), m$mean, m$sd, log = TRUE)
}

#' @export
cpd_log_density.pgm_cg_cpd <- function(x, data, model, ...) {
  m <- lg_mean(x, data)
  dnorm(as.numeric(data[[x$variable]]), m$mean, m$sd, log = TRUE)
}

#' Canonical factor of a linear Gaussian CPD
#'
#' @param x A `pgm_lg_cpd`.
#' @return A `pgm_gaussian` over the child and its parents.
#' @export
lg_cpd_factor <- function(x) {
  b <- c(1, -as.numeric(x$coefficients))
  vars <- c(x$variable, names(x$coefficients))
  s2 <- x$sd^2
  K <- outer(b, b) / s2
  h <- b * x$intercept / s2
  g <- -x$intercept^2 / (2 * s2) - 0.5 * log(2 * pi * s2)
  new_gaussian(vars, K, h, g)
}

# Joint normal distribution of a linear Gaussian network (or of the continuous
# part of a CG network for a fixed discrete configuration).
lg_joint <- function(cpds, order, conf_of = NULL) {
  k <- length(order)
  B <- matrix(0, k, k, dimnames = list(order, order))
  b0 <- setNames(numeric(k), order)
  s2 <- setNames(numeric(k), order)
  for (v in order) {
    x <- cpds[[v]]
    if (inherits(x, "pgm_cg_cpd")) {
      i <- conf_of[[v]]
      b0[v] <- x$intercept[i]
      s2[v] <- x$sd[i]^2
      for (p in x$continuous_parents) B[v, p] <- x$coefficients[i, p]
    } else {
      b0[v] <- x$intercept
      s2[v] <- x$sd^2
      for (p in names(x$coefficients)) B[v, p] <- x$coefficients[[p]]
    }
  }
  # (I - B) is unit lower triangular in topological order.
  A <- forwardsolve(diag(k) - B, diag(k))
  dimnames(A) <- list(order, order)
  mu <- as.numeric(A %*% b0)
  S <- A %*% (s2 * t(A))
  list(mean = setNames(mu, order), covariance = unname_dims((S + t(S)) / 2, order))
}

#' Joint distribution of a linear Gaussian network
#'
#' @param model A linear Gaussian `bayesian_network`.
#' @return A `pgm_gaussian` in moment form over all nodes.
#' @examples
#' gaussian_moments(gaussian_joint(example_network("gaussian")))
#' @export
gaussian_joint <- function(model) {
  if (!is_gaussian_bn(model)) pgm_abort("`model` must be a linear Gaussian Bayesian network.")
  check_model(model)
  ord <- topological_sort(model$graph)
  j <- lg_joint(model$cpds, ord)
  gaussian_factor(model$graph$nodes, j$mean[model$graph$nodes], j$covariance[model$graph$nodes, model$graph$nodes])
}

#' Exact queries in Gaussian and conditional Gaussian networks
#'
#' For a linear Gaussian network the posterior of the query variables given
#' continuous evidence is multivariate normal and is computed from the joint
#' normal distribution by Cholesky based conditioning. For a conditional
#' Gaussian network (discrete and continuous nodes, no discrete node with a
#' continuous parent) the posterior is computed exactly by enumerating the
#' joint configurations of the discrete variables that are consistent with
#' the discrete evidence: the posterior of discrete query variables is a
#' table and the posterior of continuous query variables is a finite mixture
#' of normals. The enumeration is exponential in the number of discrete
#' variables and is limited by `max_configurations`.
#'
#' Interventions on continuous variables cut their incoming arcs and fix
#' their values (graph mutilation), then condition on the fixed values.
#'
#' @param model A Gaussian or conditional Gaussian `bayesian_network`.
#' @param variables Query variables.
#' @param evidence Named list of observed values (numbers for continuous
#'   variables, states for discrete ones).
#' @param interventions Named list of values set by intervention.
#' @param max_configurations Upper bound on discrete configurations.
#' @return For a linear Gaussian network, a `pgm_gaussian` (moment form). For
#'   a conditional Gaussian network, an object of class `pgm_cg_posterior`
#'   with elements `discrete` (a `pgm_factor` or `NULL`), `mixture` (a data
#'   frame of component weights and the discrete configuration),
#'   `components` (list of means and covariances), `mean` and `covariance`
#'   (overall moments of the continuous query variables) and
#'   `evidence_density`.
#' @examples
#' g <- example_network("gaussian")
#' gaussian_query(g, "X1", evidence = list(X3 = 2))
#' @export
gaussian_query <- function(model, variables, evidence = list(), interventions = NULL, max_configurations = 1e5) {
  evidence <- if (length(evidence)) as.list(evidence) else list()
  interventions <- if (length(interventions)) as.list(interventions) else list()
  check_model(model)
  bad <- setdiff(c(variables, names(evidence), names(interventions)), nodes(model))
  if (length(bad)) pgm_abort(paste0("Unknown variable(s): ", paste(bad, collapse = ", "), "."), class = "pgm_unknown_variable")
  if (length(intersect(variables, c(names(evidence), names(interventions))))) {
    pgm_abort("Query variables cannot also be observed or intervened on.")
  }
  if (length(interventions)) {
    model <- intervene(model, interventions)
    evidence <- c(evidence[setdiff(names(evidence), names(interventions))], interventions)
  }
  reg <- model$registry
  is_cont <- vapply(reg[nodes(model)], function(r) identical(r$type, "continuous"), logical(1))
  cvars <- names(is_cont)[is_cont]
  dvars <- names(is_cont)[!is_cont]
  ev_c <- evidence[intersect(names(evidence), cvars)]
  ev_d <- evidence[intersect(names(evidence), dvars)]
  for (v in names(ev_c)) {
    if (!is.numeric(ev_c[[v]]) || length(ev_c[[v]]) != 1L || !is.finite(ev_c[[v]])) {
      pgm_abort(sprintf("Evidence on continuous '%s' must be a finite number.", v))
    }
  }
  if (!length(dvars)) {
    j <- lg_joint(model$cpds, topological_sort(model$graph))
    f <- gaussian_factor(names(j$mean), j$mean, j$covariance)
    if (length(ev_c)) f <- gaussian_reduce(f, unlist(ev_c))
    out <- gaussian_marginalize(f, setdiff(f$vars, variables))
    out <- gaussian_factor(variables, gaussian_moments(out)$mean[variables],
                           gaussian_moments(out)$covariance[variables, variables, drop = FALSE])
    return(out)
  }
  # Conditional Gaussian: enumerate discrete configurations.
  ev_d <- normalize_evidence(ev_d)
  dstates <- lapply(setNames(dvars, dvars), function(v) if (v %in% names(ev_d)) as.character(ev_d[[v]]) else reg[[v]]$states)
  for (v in names(ev_d)) if (!ev_d[[v]] %in% reg[[v]]$states) pgm_abort(sprintf("'%s' is not a state of '%s'.", ev_d[[v]], v), class = "pgm_unknown_state")
  ncfg <- prod(vapply(dstates, length, integer(1)))
  if (ncfg > max_configurations) pgm_abort(sprintf("Exact conditional Gaussian inference needs %g discrete configurations (limit %g).", ncfg, max_configurations))
  grid <- expand.grid(dstates, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  ord <- topological_sort(model$graph)
  cord <- ord[ord %in% cvars]
  qc <- intersect(variables, cvars)
  qd <- intersect(variables, dvars)
  logw <- numeric(nrow(grid))
  comps <- vector("list", nrow(grid))
  for (r in seq_len(nrow(grid))) {
    row <- grid[r, , drop = FALSE]
    lp <- 0
    for (v in dvars) lp <- lp + log(cpd_density(model$cpds[[v]], row, model))
    conf_of <- list()
    for (v in cord) if (inherits(model$cpds[[v]], "pgm_cg_cpd")) conf_of[[v]] <- cg_config_index(model$cpds[[v]], row)
    comp <- NULL
    if (length(cord)) {
      j <- lg_joint(model$cpds, cord, conf_of)
      if (length(ev_c)) {
        y <- names(ev_c)
        R <- safe_chol(j$covariance[y, y, drop = FALSE])
        d <- unlist(ev_c)[y] - j$mean[y]
        z <- forwardsolve(t(R), d)
        lp <- lp - 0.5 * sum(z^2) - sum(log(diag(R))) - length(y) / 2 * log(2 * pi)
      }
      if (length(qc)) {
        f <- gaussian_factor(cord, j$mean, j$covariance)
        if (length(ev_c)) f <- gaussian_reduce(f, unlist(ev_c))
        m <- gaussian_moments(gaussian_marginalize(f, setdiff(f$vars, qc)))
        comp <- list(mean = m$mean[qc], covariance = m$covariance[qc, qc, drop = FALSE])
      }
    }
    logw[r] <- lp
    comps[[r]] <- comp
  }
  if (all(logw == -Inf)) impossible_evidence(evidence)
  lz <- log_sum_exp(logw)
  w <- exp(logw - lz)
  discrete <- NULL
  if (length(qd)) {
    key <- do.call(paste, c(unname(as.list(grid[qd])), sep = "\r"))
    qgrid <- expand.grid(reg_states(reg, qd), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    qkey <- do.call(paste, c(unname(as.list(qgrid)), sep = "\r"))
    vals <- vapply(qkey, function(k) sum(w[key == k]), numeric(1))
    discrete <- new_factor(qd, reg_states(reg, qd), vals)
  }
  mean_c <- cov_c <- NULL
  if (length(qc)) {
    keep <- w > 0
    mean_c <- Reduce(`+`, Map(function(c, wi) wi * c$mean, comps[keep], w[keep]))
    cov_c <- Reduce(`+`, Map(function(c, wi) wi * (c$covariance + tcrossprod(c$mean - mean_c)), comps[keep], w[keep]))
  }
  structure(list(discrete = discrete, mixture = cbind(grid, weight = w), components = comps,
                 mean = mean_c, covariance = cov_c, evidence_density = exp(lz), variables = variables),
            class = "pgm_cg_posterior")
}

reg_states <- function(reg, vars) lapply(setNames(vars, vars), function(v) reg[[v]]$states)

#' @export
print.pgm_cg_posterior <- function(x, ...) {
  cat("<Conditional Gaussian posterior>\n")
  if (!is.null(x$discrete)) print(x$discrete)
  if (!is.null(x$mean)) {
    cat(sprintf("Continuous: mixture of %d normal components\n  mean:\n", sum(x$mixture$weight > 0)))
    print(round(x$mean, 6))
    cat("  covariance:\n")
    print(round(x$covariance, 6))
  }
  invisible(x)
}

# Gaussian Markov networks --------------------------------------------------------

#' Gaussian graphical models
#'
#' `gaussian_markov_network()` defines a Gaussian Markov random field by its
#' precision matrix J and potential vector h (density proportional to
#' exp(-x'Jx/2 + h'x)); zeros in J are missing edges, that is conditional
#' independencies. `gaussian_bp()` runs Gaussian belief propagation
#' (scalar Gaussian messages). On tree structured graphs it is exact; on
#' loopy graphs, when it converges, the means are exact and the variances
#' approximate (Weiss and Freeman, 2001). `gaussian_markov_network()`
#' requires a positive definite J.
#'
#' @param precision Symmetric positive definite precision matrix with
#'   dimnames.
#' @param h Potential vector (default zero, so the mean is zero).
#' @param model A `gaussian_mn`.
#' @param max_iter,tol,damping Iteration controls for `gaussian_bp()`.
#' @return `gaussian_markov_network()` returns a `gaussian_mn`;
#'   `gaussian_bp()` returns a list with `mean`, `variance`, `converged`,
#'   `iterations`.
#' @examples
#' J <- matrix(c(2, -1, 0, -1, 2, -1, 0, -1, 2), 3,
#'             dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
#' gmn <- gaussian_markov_network(J, h = c(1, 0, 0))
#' gaussian_bp(gmn)$mean
#' gaussian_moments(gmn$factor)$mean
#' @export
gaussian_markov_network <- function(precision, h = NULL) {
  precision <- as.matrix(precision)
  vars <- rownames(precision) %||% colnames(precision) %||% paste0("X", seq_len(nrow(precision)))
  h <- h %||% numeric(length(vars))
  safe_chol(precision, "precision matrix")
  f <- new_gaussian(vars, precision, h, 0)
  from <- to <- character()
  for (i in seq_along(vars)) for (j in seq_len(i - 1L)) if (precision[i, j] != 0) {
    from <- c(from, vars[j]); to <- c(to, vars[i])
  }
  structure(list(factor = f, graph = undirected_from_pairs(vars, from, to),
                 registry = lapply(setNames(vars, vars), function(v) list(type = "continuous", states = NULL))),
            class = c("gaussian_mn", "pgm_model"))
}

#' @export
print.gaussian_mn <- function(x, ...) {
  cat(sprintf("<Gaussian Markov network: %d variables, %d edges>\n", length(x$graph$nodes), nrow(x$graph$edges)))
  invisible(x)
}

#' @rdname gaussian_markov_network
#' @export
gaussian_bp <- function(model, max_iter = 500L, tol = 1e-12, damping = 0) {
  J <- model$factor$K
  h <- model$factor$h
  n <- nrow(J)
  vars <- rownames(J)
  P <- matrix(0, n, n) # message precisions P[i, j] from i to j
  M <- matrix(0, n, n) # message means
  nbrs <- lapply(seq_len(n), function(i) setdiff(which(J[i, ] != 0), i))
  converged <- FALSE
  it <- 0L
  for (it in seq_len(max_iter)) {
    P_new <- P; M_new <- M
    for (i in seq_len(n)) for (j in nbrs[[i]]) {
      others <- setdiff(nbrs[[i]], j)
      Pi <- J[i, i] + sum(P[others, i])
      hi <- h[i] + sum(P[others, i] * M[others, i])
      if (Pi <= 0) pgm_abort("Gaussian BP produced a nonpositive precision; the model is not walk summable.", class = "pgm_convergence_failure")
      P_new[i, j] <- -J[i, j]^2 / Pi
      M_new[i, j] <- if (P_new[i, j] != 0) (-J[i, j] * hi / Pi) / P_new[i, j] else 0
    }
    if (damping > 0) {
      P_new <- (1 - damping) * P_new + damping * P
      M_new <- (1 - damping) * M_new + damping * M
    }
    delta <- max(abs(P_new - P), abs(M_new - M))
    P <- P_new; M <- M_new
    if (delta < tol) { converged <- TRUE; break }
  }
  prec <- diag(J) + colSums(P)
  hm <- h + colSums(P * M)
  if (!converged) pgm_warn("Gaussian belief propagation did not converge.", class = "pgm_convergence_warning")
  list(mean = setNames(hm / prec, vars), variance = setNames(1 / prec, vars), converged = converged, iterations = it,
       exact = converged && nrow(model$graph$edges) == n - count_components(model$graph))
}

count_components <- function(g) {
  adj <- adjacency_sets(g)
  seen <- character()
  k <- 0L
  for (v in g$nodes) {
    if (v %in% seen) next
    k <- k + 1L
    fr <- v
    seen <- c(seen, v)
    while (length(fr)) {
      nx <- setdiff(unique(unlist(adj[fr], use.names = FALSE)), seen)
      seen <- c(seen, nx)
      fr <- nx
    }
  }
  k
}
