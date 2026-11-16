# Visualization --------------------------------------------------------------------------

#' Plot graphs and models
#'
#' Draws a graph, Bayesian network, Markov network, factor graph, junction
#' tree or dynamic Bayesian network with 'igraph' (a suggested package).
#' Evidence nodes are filled in one color and intervention nodes in another;
#' `edge_strength` (for example from [bootstrap_structure()]) sets edge width
#' and transparency to show uncertainty.
#'
#' @param x A graph or model.
#' @param evidence,interventions Character vectors (or named lists) of nodes
#'   to highlight.
#' @param edge_strength Optional data frame with `from`, `to`, `strength`.
#' @param layout An 'igraph' layout function.
#' @param ... Passed to `igraph::plot.igraph()`.
#' @return The 'igraph' object, invisibly.
#' @examples
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   plot(example_network("asia"), evidence = "xray", interventions = "smoke")
#' }
#' @export
plot.pgm_graph <- function(x, evidence = NULL, interventions = NULL, edge_strength = NULL, layout = NULL, ...) {
  if (!requireNamespace("igraph", quietly = TRUE)) pgm_abort("Plotting needs the 'igraph' package.")
  ig <- as_igraph(x)
  nm <- igraph::V(ig)$name
  ev <- if (is.list(evidence)) names(evidence) else evidence
  iv <- if (is.list(interventions)) names(interventions) else interventions
  col <- rep("#EAF1FB", length(nm))
  col[nm %in% ev] <- "#F2A93B"
  col[nm %in% iv] <- "#8FA2FF"
  shape <- rep("circle", length(nm))
  if (!is.null(attr(x, "factor_nodes"))) shape[nm %in% attr(x, "factor_nodes")] <- "square"
  width <- rep(1.5, igraph::ecount(ig))
  ecol <- rep("grey35", igraph::ecount(ig))
  if (!is.null(edge_strength)) {
    el <- igraph::as_edgelist(ig)
    k <- paste(el[, 1], el[, 2])
    s <- edge_strength$strength[match(k, paste(edge_strength$from, edge_strength$to))]
    s[is.na(s)] <- 0
    width <- 0.5 + 4 * s
    ecol <- grDevices::adjustcolor("grey20", alpha.f = 0.25 + 0.75 * max(0, 1) * 1)
    ecol <- vapply(s, function(a) grDevices::adjustcolor("grey20", alpha.f = 0.2 + 0.8 * a), "")
  }
  arrows <- if (igraph::is_directed(ig)) {
    tp <- igraph::E(ig)$type
    if (is.null(tp)) 0.4 else ifelse(tp == "undirected", 0, 0.4)
  } else 0
  lay <- layout %||% if (inherits(as_pgm_graph(x), "pgm_dag")) igraph::layout_with_sugiyama(ig)$layout else igraph::layout_with_fr(ig)
  igraph::plot.igraph(ig, layout = lay, vertex.color = col, vertex.shape = shape, vertex.label.color = "black",
                      vertex.frame.color = "grey30", vertex.size = 28, edge.width = width, edge.color = ecol,
                      edge.arrow.size = arrows, ...)
  invisible(ig)
}

#' @rdname plot.pgm_graph
#' @export
plot.pgm_model <- function(x, ...) {
  if (inherits(x, "factor_graph")) {
    g <- x$graph
    attr(g, "factor_nodes") <- names(x$factors)
    return(plot.pgm_graph(g, ...))
  }
  if (inherits(x, "junction_tree")) {
    labs <- vapply(x$cliques, paste, "", collapse = ",")
    labs <- make.unique(labs)
    e <- if (nrow(x$tree)) list(from = labs[x$tree[, 1]], to = labs[x$tree[, 2]]) else NULL
    g <- if (is.null(e)) undirected_graph(nodes = labs) else undirected_graph(data.frame(from = e$from, to = e$to), nodes = labs)
    return(plot.pgm_graph(g, ...))
  }
  if (inherits(x, "dynamic_bn")) return(plot.pgm_graph(transition_graph(x), ...))
  if (inherits(x, "pgm_scm")) return(plot.pgm_graph(x$bn$graph, ...))
  plot.pgm_graph(as_pgm_graph(x), ...)
}

#' Plot a discrete factor or CPD
#'
#' Bar chart of a one variable factor, or of each column of a CPD (one panel
#' group per parent configuration), using base graphics.
#'
#' @param x A `pgm_factor`.
#' @param ... Passed to `graphics::barplot()`.
#' @return `x`, invisibly.
#' @export
plot.pgm_factor <- function(x, ...) {
  f <- factor_to_linear(x)
  if (length(f$vars) == 1L) {
    graphics::barplot(setNames(f$values, f$states[[1L]]), main = f$vars, ylab = "value", col = "#8FA2FF", ...)
  } else {
    m <- matrix(f$values, nrow = f$card[[1L]], dimnames = list(f$states[[1L]], NULL))
    graphics::barplot(m, beside = TRUE, legend.text = rownames(m), main = paste(f$vars, collapse = ", "),
                      col = grDevices::hcl.colors(nrow(m), "Blue-Red 3"), ...)
  }
  invisible(x)
}
