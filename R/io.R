# Serialization --------------------------------------------------------------------------

require_discrete_bn <- function(model, format) {
  if (!inherits(model, "bayesian_network")) pgm_abort(sprintf("%s export supports Bayesian networks.", format))
  check_model(model)
  bad <- names(Filter(function(x) !inherits(x, "pgm_cpd"), model$cpds))
  if (length(bad)) {
    pgm_abort(sprintf("%s supports tabular CPDs only; not representable: %s. Use JSON or RDS.", format, paste(bad, collapse = ", ")),
              class = "pgm_unsupported_feature")
  }
}

bif_name_ok <- function(x) grepl("^[A-Za-z_][A-Za-z0-9_.-]*$", x)

#' Read and write models
#'
#' * BIF (`write_bif()`, `read_bif()`): the Bayesian Interchange Format used
#'   by 'bnlearn', 'pgmpy' and the Bayesian network repository. Tabular CPDs
#'   only.
#' * XMLBIF (`write_xmlbif()`, `read_xmlbif()`): XML version 0.3. Reading
#'   needs the 'xml2' package. Tabular CPDs only.
#' * UAI (`write_uai()`, `read_uai()`): the UAI competition format for
#'   Bayesian networks (`BAYES`) and Markov networks (`MARKOV`). Variable and
#'   state names are not part of the format; on reading, variables are named
#'   `V0`, `V1`, ... and states `"0"`, `"1"`, ... unless `names` and `states`
#'   are given.
#' * JSON (`write_json_model()`, `read_json_model()`): a lossless pgmcompR
#'   schema for discrete, linear Gaussian and conditional Gaussian networks
#'   and Markov networks. Needs 'jsonlite'.
#' * RDS (`save_model()`, `load_model()`): R serialization of any model.
#'
#' Features a format cannot represent raise an error of class
#' `pgm_unsupported_feature` rather than being dropped.
#'
#' @param model A model.
#' @param file Path.
#' @param names,states Optional variable names and state labels for UAI.
#' @return Writers return `file` invisibly; readers return a model.
#' @examples
#' bn <- example_network("asia")
#' f <- tempfile(fileext = ".bif")
#' write_bif(bn, f)
#' bn2 <- read_bif(f)
#' all.equal(query(bn, "lung", list(xray = "yes")), query(bn2, "lung", list(xray = "yes")))
#' @export
write_bif <- function(model, file) {
  require_discrete_bn(model, "BIF")
  vars <- nodes(model)
  bad <- c(vars[!bif_name_ok(vars)], unlist(lapply(model$registry[vars], function(r) r$states[!grepl("^[^,;{}()|\\s]+$", r$states, perl = TRUE)])))
  if (length(bad)) pgm_abort(sprintf("BIF cannot represent these names: %s.", paste(unique(bad), collapse = ", ")), class = "pgm_unsupported_feature")
  out <- c("network unknown {", "}")
  for (v in vars) {
    st <- model$registry[[v]]$states
    out <- c(out, sprintf("variable %s {", v), sprintf("  type discrete [ %d ] { %s };", length(st), paste(st, collapse = ", ")), "}")
  }
  for (v in vars) {
    x <- model$cpds[[v]]
    m <- matrix(x$values, nrow = x$card[[1L]])
    fmt <- function(p) paste(formatC(p, digits = 17, format = "g"), collapse = ", ")
    if (!length(x$parents)) {
      out <- c(out, sprintf("probability ( %s ) {", v), sprintf("  table %s;", fmt(m[, 1])), "}")
    } else {
      out <- c(out, sprintf("probability ( %s | %s ) {", v, paste(x$parents, collapse = ", ")))
      grid <- expand.grid(x$states[x$parents], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      for (j in seq_len(nrow(grid))) {
        out <- c(out, sprintf("  (%s) %s;", paste(unlist(grid[j, ]), collapse = ", "), fmt(m[, j])))
      }
      out <- c(out, "}")
    }
  }
  writeLines(out, file)
  invisible(file)
}

#' @rdname write_bif
#' @export
read_bif <- function(file) {
  txt <- paste(readLines(file, warn = FALSE), collapse = "\n")
  txt <- gsub("//[^\n]*", "", txt)
  txt <- gsub("/\\*.*?\\*/", "", txt, perl = TRUE)
  vblocks <- regmatches(txt, gregexpr("variable\\s+[^\\s{]+\\s*\\{[^}]*\\}[^}]*\\}", txt, perl = TRUE))[[1]]
  states <- list()
  for (b in vblocks) {
    nm <- sub("(?s)variable\\s+([^\\s{]+).*", "\\1", b, perl = TRUE)
    st <- sub("(?s).*type\\s+discrete\\s*\\[\\s*\\d+\\s*\\]\\s*\\{([^}]*)\\}.*", "\\1", b, perl = TRUE)
    if (identical(st, b)) pgm_abort(sprintf("Variable '%s' is not discrete; only discrete BIF is supported.", nm), class = "pgm_unsupported_feature")
    states[[nm]] <- trimws(strsplit(st, ",")[[1]])
  }
  pblocks <- regmatches(txt, gregexpr("probability\\s*\\([^)]*\\)\\s*\\{[^}]*\\}", txt, perl = TRUE))[[1]]
  cps <- list(); edges <- list()
  for (b in pblocks) {
    head <- sub("(?s)probability\\s*\\(([^)]*)\\).*", "\\1", b, perl = TRUE)
    parts <- strsplit(head, "|", fixed = TRUE)[[1]]
    child <- trimws(parts[1])
    pa <- if (length(parts) > 1L) trimws(strsplit(parts[2], ",")[[1]]) else character()
    body <- sub("^[^{]*\\{", "", b)
    body <- sub("\\}\\s*$", "", body)
    st <- states[c(child, pa)]
    q <- prod(vapply(st[pa], length, integer(1)))
    m <- matrix(NA_real_, length(st[[child]]), max(q, 1))
    if (grepl("table", body)) {
      vals <- as.numeric(trimws(strsplit(sub("(?s).*table\\s*([^;]*);.*", "\\1", body, perl = TRUE), ",")[[1]]))
      # BIF tables list the child fastest
      m[] <- vals
    }
    rows <- regmatches(body, gregexpr("\\(([^)]*)\\)\\s*([^;]*);", body, perl = TRUE))[[1]]
    for (r in rows) {
      cfg <- trimws(strsplit(sub("\\(([^)]*)\\).*", "\\1", r), ",")[[1]])
      vals <- as.numeric(trimws(strsplit(sub(".*\\)\\s*([^;]*);", "\\1", r), ",")[[1]]))
      j <- 1; stride <- 1
      for (k in seq_along(pa)) {
        j <- j + (match(cfg[k], st[[pa[k]]]) - 1) * stride
        stride <- stride * length(st[[pa[k]]])
      }
      m[, j] <- vals
    }
    if (anyNA(m)) pgm_abort(sprintf("The BIF table of '%s' is incomplete.", child))
    cps[[child]] <- cpd(child, as.numeric(m), parents = pa, states = st, normalize = FALSE)
    for (p in pa) edges[[length(edges) + 1L]] <- c(p, child)
  }
  bn <- bayesian_network(edges, nodes = names(states), states = states)
  for (c in cps) bn <- add_cpd(bn, c)
  bn
}

xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

#' @rdname write_bif
#' @export
write_xmlbif <- function(model, file) {
  require_discrete_bn(model, "XMLBIF")
  vars <- nodes(model)
  out <- c('<?xml version="1.0" encoding="UTF-8"?>', "<BIF VERSION=\"0.3\">", "<NETWORK>", "<NAME>pgmcompR</NAME>")
  for (v in vars) {
    out <- c(out, '<VARIABLE TYPE="nature">', sprintf("  <NAME>%s</NAME>", xml_escape(v)),
             sprintf("  <OUTCOME>%s</OUTCOME>", xml_escape(model$registry[[v]]$states)), "</VARIABLE>")
  }
  for (v in vars) {
    x <- model$cpds[[v]]
    # XMLBIF tables list the last GIVEN fastest and the child fastest overall.
    arr <- array(x$values, dim = unname(x$card))
    perm <- c(1L, rev(seq_along(x$parents)) + 1L)
    vals <- as.numeric(aperm(arr, perm))
    out <- c(out, "<DEFINITION>", sprintf("  <FOR>%s</FOR>", xml_escape(v)),
             if (length(x$parents)) sprintf("  <GIVEN>%s</GIVEN>", xml_escape(x$parents)),
             sprintf("  <TABLE>%s</TABLE>", paste(formatC(vals, digits = 17, format = "g"), collapse = " ")), "</DEFINITION>")
  }
  out <- c(out, "</NETWORK>", "</BIF>")
  writeLines(out, file)
  invisible(file)
}

#' @rdname write_bif
#' @export
read_xmlbif <- function(file) {
  if (!requireNamespace("xml2", quietly = TRUE)) pgm_abort("Reading XMLBIF needs the 'xml2' package.")
  doc <- xml2::read_xml(file)
  vnodes <- xml2::xml_find_all(doc, ".//VARIABLE")
  states <- list()
  for (vn in vnodes) {
    nm <- xml2::xml_text(xml2::xml_find_first(vn, "./NAME"))
    states[[nm]] <- xml2::xml_text(xml2::xml_find_all(vn, "./OUTCOME"))
  }
  defs <- xml2::xml_find_all(doc, ".//DEFINITION")
  cps <- list(); edges <- list()
  for (d in defs) {
    child <- xml2::xml_text(xml2::xml_find_first(d, "./FOR"))
    pa <- xml2::xml_text(xml2::xml_find_all(d, "./GIVEN"))
    vals <- as.numeric(strsplit(trimws(xml2::xml_text(xml2::xml_find_first(d, "./TABLE"))), "\\s+")[[1]])
    st <- states[c(child, pa)]
    card <- vapply(st, length, integer(1))
    if (length(vals) != prod(card)) pgm_abort(sprintf("The XMLBIF table of '%s' has the wrong length.", child))
    # stored with child fastest, then parents with the last parent fastest
    arr <- array(vals, dim = c(card[1], rev(card[-1])))
    perm <- c(1L, rev(seq_along(pa)) + 1L)
    v2 <- as.numeric(aperm(arr, order(perm)))
    cps[[child]] <- cpd(child, v2, parents = pa, states = st)
    for (p in pa) edges[[length(edges) + 1L]] <- c(p, child)
  }
  bn <- bayesian_network(edges, nodes = names(states), states = states)
  for (c in cps) bn <- add_cpd(bn, c)
  bn
}

#' @rdname write_bif
#' @export
write_uai <- function(model, file) {
  if (inherits(model, "bayesian_network")) {
    require_discrete_bn(model, "UAI")
    fs <- lapply(model$cpds, function(x) x)
    type <- "BAYES"
  } else if (inherits(model, c("markov_network", "factor_graph"))) {
    fs <- model_factors(model)
    type <- "MARKOV"
  } else pgm_abort("UAI export needs a discrete Bayesian network, Markov network or factor graph.")
  vars <- nodes(model)
  card <- vapply(model_states(model)[vars], length, integer(1))
  out <- c(type, length(vars), paste(card, collapse = " "), length(fs))
  scope_line <- function(f) {
    sc <- if (inherits(f, "pgm_cpd")) c(f$parents, f$child) else f$vars
    paste(c(length(sc), match(sc, vars) - 1L), collapse = " ")
  }
  out <- c(out, vapply(fs, scope_line, ""))
  for (f in fs) {
    sc <- if (inherits(f, "pgm_cpd")) c(f$parents, f$child) else f$vars
    # UAI tables list the last variable of the scope fastest.
    g <- factor_reorder(as_plain_factor(f), rev(sc))
    out <- c(out, "", length(g$values), paste(formatC(g$values, digits = 17, format = "g"), collapse = " "))
  }
  writeLines(as.character(out), file)
  invisible(file)
}

#' @rdname write_bif
#' @export
read_uai <- function(file, names = NULL, states = NULL) {
  tok <- scan(file, what = "", quiet = TRUE)
  pos <- 1L
  nxt <- function() { v <- tok[pos]; pos <<- pos + 1L; v }
  type <- nxt()
  n <- as.integer(nxt())
  card <- as.integer(vapply(seq_len(n), function(i) nxt(), ""))
  vars <- names %||% paste0("V", seq_len(n) - 1L)
  st <- states %||% lapply(card, function(k) as.character(seq_len(k) - 1L))
  names(st) <- vars
  nf <- as.integer(nxt())
  scopes <- lapply(seq_len(nf), function(i) {
    k <- as.integer(nxt())
    vars[as.integer(vapply(seq_len(k), function(j) nxt(), "")) + 1L]
  })
  fs <- lapply(scopes, function(sc) {
    m <- as.integer(nxt())
    vals <- as.numeric(vapply(seq_len(m), function(j) nxt(), ""))
    f <- new_factor(rev(sc), st[rev(sc)], vals)
    factor_reorder(f, sc)
  })
  if (toupper(type) == "BAYES") {
    edges <- list(); cps <- list()
    for (i in seq_along(scopes)) {
      sc <- scopes[[i]]
      child <- sc[length(sc)]
      pa <- sc[-length(sc)]
      f <- factor_reorder(fs[[i]], c(child, pa))
      cps[[child]] <- cpd(child, f$values, parents = pa, states = st[c(child, pa)], normalize = TRUE)
      for (p in pa) edges[[length(edges) + 1L]] <- c(p, child)
    }
    bn <- bayesian_network(edges, nodes = vars, states = st)
    for (c in cps) bn <- add_cpd(bn, c)
    return(bn)
  }
  markov_network(fs, nodes = vars, states = st)
}

#' @rdname write_bif
#' @export
write_json_model <- function(model, file) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) pgm_abort("JSON serialization needs the 'jsonlite' package.")
  enc_cpd <- function(x) {
    if (inherits(x, "pgm_cpd")) return(list(type = "tabular", variable = x$child, parents = x$parents, values = x$values))
    if (inherits(x, "pgm_lg_cpd")) return(list(type = "linear_gaussian", variable = x$variable, parents = x$parents,
                                               intercept = x$intercept, coefficients = as.list(x$coefficients), sd = x$sd))
    if (inherits(x, "pgm_cg_cpd")) return(list(type = "conditional_gaussian", variable = x$variable,
                                               discrete_parents = x$discrete_parents, continuous_parents = x$continuous_parents,
                                               intercept = x$intercept, coefficients = unname(as.list(as.data.frame(x$coefficients))),
                                               sd = x$sd))
    pgm_abort(sprintf("JSON cannot represent CPD class '%s'.", class(x)[1L]), class = "pgm_unsupported_feature")
  }
  reg <- lapply(model$registry, function(r) list(type = r$type, states = r$states))
  obj <- if (inherits(model, "bayesian_network")) {
    list(format = "pgmcompR", version = 1L, class = "bayesian_network", nodes = nodes(model),
         edges = list(from = model$graph$edges$from, to = model$graph$edges$to), registry = reg,
         cpds = unname(lapply(model$cpds, enc_cpd)))
  } else if (inherits(model, "markov_network")) {
    list(format = "pgmcompR", version = 1L, class = "markov_network", nodes = nodes(model), registry = reg,
         factors = lapply(model$factors, function(f) list(variables = f$vars, values = f$values, log = f$log)))
  } else pgm_abort("JSON export supports Bayesian and Markov networks.", class = "pgm_unsupported_feature")
  jsonlite::write_json(obj, file, auto_unbox = TRUE, digits = I(17), null = "null", pretty = TRUE)
  invisible(file)
}

#' @rdname write_bif
#' @export
read_json_model <- function(file) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) pgm_abort("JSON serialization needs the 'jsonlite' package.")
  o <- jsonlite::read_json(file, simplifyVector = TRUE)
  if (!identical(o$format, "pgmcompR")) pgm_abort("Not a pgmcompR JSON model.")
  reg <- o$registry
  st <- lapply(reg, function(r) if (identical(r$type, "discrete")) as.character(r$states) else NULL)
  st <- st[!vapply(st, is.null, logical(1))]
  if (o$class == "markov_network") {
    fs <- lapply(seq_along(o$factors$variables), function(i) {
      v <- o$factors$variables[[i]]
      new_factor(v, st[v], o$factors$values[[i]], isTRUE(o$factors$log[[i]]))
    })
    return(markov_network(fs, nodes = o$nodes, states = st))
  }
  e <- if (length(o$edges$from)) Map(c, o$edges$from, o$edges$to) else NULL
  bn <- bayesian_network(unname(e), nodes = o$nodes, states = st)
  cp <- o$cpds
  for (i in seq_len(nrow(cp))) {
    type <- cp$type[i]; v <- cp$variable[i]
    if (type == "tabular") {
      pa <- as.character(cp$parents[[i]])
      bn <- add_cpd(bn, cpd(v, cp$values[[i]], parents = pa, states = st[c(v, pa)]))
    } else if (type == "linear_gaussian") {
      co <- unlist(cp$coefficients[i, , drop = TRUE])
      co <- if (length(co)) co[!is.na(co)] else NULL
      bn <- add_cpd(bn, lg_cpd(v, cp$intercept[[i]], if (length(co)) co else NULL, cp$sd[[i]]))
    } else {
      dp <- as.character(cp$discrete_parents[[i]]); cpa <- as.character(cp$continuous_parents[[i]])
      co <- cp$coefficients[[i]]
      cm <- if (length(cpa)) do.call(cbind, lapply(co, as.numeric)) else matrix(0, length(cp$intercept[[i]]), 0)
      if (length(cpa)) colnames(cm) <- cpa
      bn <- add_cpd(bn, cg_cpd(v, dp, st[dp], as.numeric(cp$intercept[[i]]), cm, as.numeric(cp$sd[[i]]), cpa))
    }
  }
  bn
}

#' @rdname write_bif
#' @export
save_model <- function(model, file) {
  saveRDS(model, file)
  invisible(file)
}

#' @rdname write_bif
#' @export
load_model <- function(file) {
  m <- readRDS(file)
  if (!inherits(m, "pgm_model") && !inherits(m, "pgm_graph")) pgm_abort("The file does not contain a pgmcompR model.")
  m
}
