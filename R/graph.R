# Graph classes ---------------------------------------------------------------

edge_frame <- function(from = character(), to = character(), type = character()) {
  data.frame(from = as.character(from), to = as.character(to),
             type = as.character(type), stringsAsFactors = FALSE)
}

parse_edges <- function(edges, type = "directed") {
  if (is.null(edges) || length(edges) == 0L) return(edge_frame())
  if (is.data.frame(edges)) {
    if (!all(c("from", "to") %in% names(edges))) {
      if (ncol(edges) < 2L) pgm_abort("`edges` data frame needs columns `from` and `to`.")
      names(edges)[1:2] <- c("from", "to")
    }
    ty <- if ("type" %in% names(edges)) edges$type else rep(type, nrow(edges))
    return(edge_frame(edges$from, edges$to, ty))
  }
  if (is.matrix(edges)) {
    if (ncol(edges) != 2L) pgm_abort("`edges` matrix must have two columns.")
    return(edge_frame(edges[, 1], edges[, 2], rep(type, nrow(edges))))
  }
  if (is.list(edges)) {
    ok <- vapply(edges, function(e) is.character(e) && length(e) == 2L, logical(1))
    if (!all(ok)) pgm_abort("Each element of `edges` must be a character vector of length 2.")
    m <- do.call(rbind, edges)
    return(edge_frame(m[, 1], m[, 2], rep(type, nrow(m))))
  }
  if (is.character(edges) && length(edges) == 2L) {
    return(edge_frame(edges[1], edges[2], type))
  }
  pgm_abort("`edges` must be a list of length-2 character vectors, a two column matrix, or a data frame.")
}

new_pgm_graph <- function(nodes, edges, subclass = NULL) {
  structure(list(nodes = nodes, edges = edges), class = c(subclass, "pgm_graph"))
}

# Canonical key for an edge, ignoring direction for undirected edges.
edge_key <- function(from, to, type) {
  und <- type == "undirected"
  a <- ifelse(und & from > to, to, from)
  b <- ifelse(und & from > to, from, to)
  paste(a, b, ifelse(und, "--", "->"), sep = "\r")
}

validate_graph <- function(g) {
  check_names(g$nodes, "nodes")
  e <- g$edges
  bad <- setdiff(c(e$from, e$to), g$nodes)
  if (length(bad)) pgm_abort(paste0("Edges refer to unknown nodes: ", paste(unique(bad), collapse = ", "), "."))
  if (any(e$from == e$to)) pgm_abort("Self loops are not allowed.")
  if (any(!e$type %in% c("directed", "undirected"))) pgm_abort("Edge type must be 'directed' or 'undirected'.")
  pair <- ifelse(e$from < e$to, paste(e$from, e$to, sep = "\r"), paste(e$to, e$from, sep = "\r"))
  if (anyDuplicated(pair)) {
    d <- pair[duplicated(pair)][1]
    pgm_abort(paste0("At most one edge is allowed between two nodes; duplicated pair: ",
                     sub("\r", " and ", d, fixed = TRUE), "."))
  }
  invisible(g)
}

#' Graphs for probabilistic graphical models
#'
#' Constructors for the graph classes used by every model in the package.
#' `pgm_graph()` builds a general graph that may mix directed and undirected
#' edges; `dag()` builds a directed acyclic graph; `undirected_graph()` builds
#' an undirected graph; `pdag()` builds a partially directed graph (for
#' example a CPDAG). Node order is the order of `nodes` followed by any new
#' node names found in `edges`, and is kept stable by every operation.
#'
#' @param edges Edges as a list of length-2 character vectors, a two column
#'   matrix, or a data frame with columns `from` and `to` (and, for
#'   `pgm_graph()` and `pdag()`, optionally `type`).
#' @param nodes Optional character vector of node names, including isolated
#'   nodes.
#' @param undirected For `pdag()`, edges that are undirected, in the same
#'   formats as `edges`.
#' @return An object of class `pgm_graph` (and `pgm_dag`, `pgm_ugraph` or
#'   `pgm_pdag`).
#' @examples
#' g <- dag(list(c("Smoking", "Cancer"), c("Cancer", "XRay")))
#' g
#' u <- undirected_graph(list(c("A", "B"), c("B", "C")))
#' @export
pgm_graph <- function(edges = NULL, nodes = NULL) {
  e <- parse_edges(edges)
  nodes <- unique(c(as.character(nodes), e$from, e$to))
  g <- new_pgm_graph(nodes, e)
  validate_graph(g)
  classify_graph(g)
}

classify_graph <- function(g) {
  ty <- g$edges$type
  sub <- if (length(ty) == 0L || all(ty == "directed")) {
    if (has_directed_cycle(g)) NULL else c("pgm_dag", "pgm_pdag")
  } else if (all(ty == "undirected")) {
    "pgm_ugraph"
  } else {
    if (has_directed_cycle(g)) NULL else "pgm_pdag"
  }
  class(g) <- c(sub, "pgm_graph")
  g
}

#' @rdname pgm_graph
#' @export
dag <- function(edges = NULL, nodes = NULL) {
  e <- parse_edges(edges, "directed")
  if (any(e$type != "directed")) pgm_abort("A DAG may only contain directed edges.")
  nodes <- unique(c(as.character(nodes), e$from, e$to))
  g <- new_pgm_graph(nodes, e, c("pgm_dag", "pgm_pdag"))
  validate_graph(g)
  if (has_directed_cycle(g)) {
    pgm_abort(paste0("The graph contains a directed cycle: ",
                     paste(find_cycle(g), collapse = " -> "), "."), class = "pgm_cycle_error")
  }
  g
}

#' @rdname pgm_graph
#' @export
undirected_graph <- function(edges = NULL, nodes = NULL) {
  e <- parse_edges(edges, "undirected")
  e$type <- rep("undirected", nrow(e))
  nodes <- unique(c(as.character(nodes), e$from, e$to))
  g <- new_pgm_graph(nodes, e, "pgm_ugraph")
  validate_graph(g)
}

#' @rdname pgm_graph
#' @export
pdag <- function(edges = NULL, undirected = NULL, nodes = NULL) {
  d <- parse_edges(edges, "directed")
  d$type <- rep("directed", nrow(d))
  u <- parse_edges(undirected, "undirected")
  u$type <- rep("undirected", nrow(u))
  e <- rbind(d, u)
  nodes <- unique(c(as.character(nodes), e$from, e$to))
  g <- new_pgm_graph(nodes, e, "pgm_pdag")
  validate_graph(g)
  if (has_directed_cycle(g)) pgm_abort("The directed part of a PDAG must be acyclic.", class = "pgm_cycle_error")
  g
}

#' @export
print.pgm_graph <- function(x, ...) {
  kind <- if (inherits(x, "pgm_dag")) "DAG" else if (inherits(x, "pgm_ugraph")) "Undirected graph" else if (inherits(x, "pgm_pdag")) "Partially directed acyclic graph" else "Graph"
  cat(sprintf("<%s: %d nodes, %d edges>\n", kind, length(x$nodes), nrow(x$edges)))
  if (nrow(x$edges)) {
    arrows <- ifelse(x$edges$type == "directed", " -> ", " -- ")
    txt <- paste0(x$edges$from, arrows, x$edges$to)
    cat(paste0("  ", txt, collapse = "\n"), "\n", sep = "")
  }
  iso <- setdiff(x$nodes, c(x$edges$from, x$edges$to))
  if (length(iso)) cat("  isolated: ", paste(iso, collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Extract the graph of a model
#'
#' @param x A graph or a model object.
#' @param ... Unused.
#' @return A `pgm_graph`.
#' @export
as_pgm_graph <- function(x, ...) UseMethod("as_pgm_graph")

#' @export
as_pgm_graph.pgm_graph <- function(x, ...) x

#' @export
as_pgm_graph.default <- function(x, ...) {
  if (is.list(x) && inherits(x$graph, "pgm_graph")) return(x$graph)
  pgm_abort("Cannot extract a graph from this object.")
}

#' Nodes and edges of a graph
#'
#' @param x A graph or model.
#' @return `nodes()` returns a character vector; `edges()` returns a data
#'   frame with columns `from`, `to` and `type`.
#' @examples
#' g <- dag(list(c("A", "B")), nodes = "C")
#' nodes(g)
#' edges(g)
#' @export
nodes <- function(x) UseMethod("nodes")

#' @export
nodes.default <- function(x) as_pgm_graph(x)$nodes

#' @rdname nodes
#' @export
edges <- function(x) as_pgm_graph(x)$edges

# Graph modification -----------------------------------------------------------

#' Modify a graph
#'
#' Add or remove nodes and edges. Removing a node removes its incident edges.
#' Adding an edge to a DAG that would create a directed cycle is an error.
#'
#' @param g A `pgm_graph`.
#' @param nodes Character vector of node names.
#' @param from,to Character vectors of equal length giving edge endpoints.
#' @param type Edge type, `"directed"` or `"undirected"`. Defaults to
#'   `"undirected"` for undirected graphs and `"directed"` otherwise.
#' @return The modified graph.
#' @examples
#' g <- dag(list(c("A", "B")))
#' g <- add_edges(g, "B", "C")
#' g <- remove_nodes(g, "A")
#' @export
add_nodes <- function(g, nodes) {
  stopifnot(inherits(g, "pgm_graph"))
  check_names(nodes)
  new <- setdiff(nodes, g$nodes)
  g$nodes <- c(g$nodes, new)
  g
}

#' @rdname add_nodes
#' @export
remove_nodes <- function(g, nodes) {
  stopifnot(inherits(g, "pgm_graph"))
  g$nodes <- setdiff(g$nodes, nodes)
  keep <- !(g$edges$from %in% nodes | g$edges$to %in% nodes)
  g$edges <- g$edges[keep, , drop = FALSE]
  rownames(g$edges) <- NULL
  g
}

#' @rdname add_nodes
#' @export
add_edges <- function(g, from, to, type = NULL) {
  stopifnot(inherits(g, "pgm_graph"))
  if (length(from) != length(to)) pgm_abort("`from` and `to` must have the same length.")
  if (is.null(type)) type <- if (inherits(g, "pgm_ugraph")) "undirected" else "directed"
  type <- rep_len(type, length(from))
  if (inherits(g, "pgm_dag") && any(type != "directed")) pgm_abort("A DAG may only contain directed edges.")
  if (inherits(g, "pgm_ugraph") && any(type != "undirected")) pgm_abort("An undirected graph may only contain undirected edges.")
  cls <- class(g)
  g <- add_nodes(g, unique(c(from, to)))
  g$edges <- rbind(g$edges, edge_frame(from, to, type))
  rownames(g$edges) <- NULL
  validate_graph(g)
  if (any(type == "directed") && has_directed_cycle(g)) {
    pgm_abort(paste0("Adding the edge(s) creates a directed cycle: ",
                     paste(find_cycle(g), collapse = " -> "), "."), class = "pgm_cycle_error")
  }
  class(g) <- cls
  g
}

#' @rdname add_nodes
#' @export
remove_edges <- function(g, from, to) {
  stopifnot(inherits(g, "pgm_graph"))
  if (length(from) != length(to)) pgm_abort("`from` and `to` must have the same length.")
  e <- g$edges
  drop <- rep(FALSE, nrow(e))
  for (i in seq_along(from)) {
    hit <- (e$from == from[i] & e$to == to[i]) |
      (e$type == "undirected" & e$from == to[i] & e$to == from[i])
    if (!any(hit)) pgm_abort(sprintf("Edge %s, %s is not in the graph.", from[i], to[i]))
    drop <- drop | hit
  }
  g$edges <- e[!drop, , drop = FALSE]
  rownames(g$edges) <- NULL
  g
}

# Neighborhoods ----------------------------------------------------------------

order_nodes <- function(g, x) g$nodes[g$nodes %in% x]

check_nodes_in <- function(g, x, arg = "nodes") {
  bad <- setdiff(x, g$nodes)
  if (length(bad)) pgm_abort(sprintf("Unknown node(s) in `%s`: %s.", arg, paste(bad, collapse = ", ")))
  invisible(x)
}

#' Neighborhoods in a graph
#'
#' `parents()`, `children()` and `neighbors()` (undirected neighbors) return
#' the direct relatives of a single node; `ancestors()` and `descendants()`
#' return all nodes reachable against or along directed edges (excluding the
#' nodes themselves); `markov_blanket()` returns parents, children and the
#' other parents of the children in a DAG, or the neighbors in an undirected
#' graph. Results follow the graph's node order.
#'
#' @param x A graph or model.
#' @param node A node name (for `ancestors()` and `descendants()`, one or more).
#' @return A character vector.
#' @examples
#' g <- dag(list(c("A", "C"), c("B", "C"), c("C", "D")))
#' parents(g, "C")
#' markov_blanket(g, "A")
#' ancestors(g, "D")
#' @export
parents <- function(x, node) {
  g <- as_pgm_graph(x)
  check_nodes_in(g, node, "node")
  e <- g$edges
  order_nodes(g, e$from[e$type == "directed" & e$to %in% node])
}

#' @rdname parents
#' @export
children <- function(x, node) {
  g <- as_pgm_graph(x)
  check_nodes_in(g, node, "node")
  e <- g$edges
  order_nodes(g, e$to[e$type == "directed" & e$from %in% node])
}

#' @rdname parents
#' @export
neighbors <- function(x, node) {
  g <- as_pgm_graph(x)
  check_nodes_in(g, node, "node")
  e <- g$edges[g$edges$type == "undirected", , drop = FALSE]
  order_nodes(g, setdiff(c(e$to[e$from %in% node], e$from[e$to %in% node]), node))
}

adjacent <- function(g, node) {
  e <- g$edges
  order_nodes(g, setdiff(c(e$to[e$from %in% node], e$from[e$to %in% node]), node))
}

#' @rdname parents
#' @export
ancestors <- function(x, node) {
  g <- as_pgm_graph(x)
  check_nodes_in(g, node, "node")
  reach_directed(g, node, upward = TRUE)
}

#' @rdname parents
#' @export
descendants <- function(x, node) {
  g <- as_pgm_graph(x)
  check_nodes_in(g, node, "node")
  reach_directed(g, node, upward = FALSE)
}

reach_directed <- function(g, start, upward) {
  e <- g$edges[g$edges$type == "directed", , drop = FALSE]
  src <- if (upward) e$to else e$from
  dst <- if (upward) e$from else e$to
  seen <- character()
  frontier <- start
  while (length(frontier)) {
    nxt <- setdiff(unique(dst[src %in% frontier]), c(seen, start))
    seen <- c(seen, nxt)
    frontier <- nxt
  }
  order_nodes(g, seen)
}

#' @rdname parents
#' @export
markov_blanket <- function(x, node) {
  g <- as_pgm_graph(x)
  check_string(node)
  check_nodes_in(g, node, "node")
  if (inherits(g, "pgm_ugraph")) return(neighbors(g, node))
  ch <- children(g, node)
  sp <- if (length(ch)) parents(g, ch) else character()
  order_nodes(g, setdiff(unique(c(parents(g, node), ch, sp, neighbors(g, node))), node))
}

# Order and cycles -------------------------------------------------------------

#' Topological order and acyclicity
#'
#' `topological_sort()` returns the nodes of a DAG so that every parent
#' precedes its children; ties are broken by the graph's node order, so the
#' result is deterministic. `has_directed_cycle()` reports whether the
#' directed edges contain a cycle and `is_dag()` whether the graph is a DAG
#' (only directed edges, no cycle).
#'
#' @param x A graph or model.
#' @return `topological_sort()` returns a character vector; the predicates
#'   return a logical scalar.
#' @examples
#' g <- dag(list(c("B", "C"), c("A", "B")))
#' topological_sort(g)
#' is_dag(g)
#' @export
topological_sort <- function(x) {
  g <- as_pgm_graph(x)
  e <- g$edges[g$edges$type == "directed", , drop = FALSE]
  indeg <- setNames(integer(length(g$nodes)), g$nodes)
  tab <- table(factor(e$to, levels = g$nodes))
  indeg[] <- as.integer(tab)
  out <- character()
  remaining <- g$nodes
  while (length(remaining)) {
    ready <- remaining[indeg[remaining] == 0L]
    if (!length(ready)) pgm_abort("The graph contains a directed cycle.", class = "pgm_cycle_error")
    v <- ready[1L]
    out <- c(out, v)
    remaining <- remaining[remaining != v]
    ch <- e$to[e$from == v]
    for (c in ch) indeg[c] <- indeg[c] - 1L
  }
  out
}

#' @rdname topological_sort
#' @export
has_directed_cycle <- function(x) {
  g <- as_pgm_graph(x)
  length(find_cycle(g)) > 0L
}

#' @rdname topological_sort
#' @export
is_dag <- function(x) {
  g <- as_pgm_graph(x)
  all(g$edges$type == "directed") && !has_directed_cycle(g)
}

# Returns one directed cycle as a node path, or character(0).
find_cycle <- function(g) {
  e <- g$edges[g$edges$type == "directed", , drop = FALSE]
  if (!nrow(e)) return(character())
  adj <- split(e$to, factor(e$from, levels = g$nodes))
  state <- setNames(integer(length(g$nodes)), g$nodes) # 0 new, 1 active, 2 done
  parent <- setNames(rep(NA_character_, length(g$nodes)), g$nodes)
  for (s in g$nodes) {
    if (state[s] != 0L) next
    stack <- s
    iter <- list()
    iter[[s]] <- 1L
    state[s] <- 1L
    while (length(stack)) {
      v <- stack[length(stack)]
      nb <- adj[[v]]
      i <- iter[[v]]
      if (i <= length(nb)) {
        iter[[v]] <- i + 1L
        w <- nb[i]
        if (state[w] == 1L) {
          path <- v
          while (path[1L] != w) path <- c(parent[path[1L]], path)
          return(unname(c(path, w)))
        }
        if (state[w] == 0L) {
          state[w] <- 1L
          parent[w] <- v
          iter[[w]] <- 1L
          stack <- c(stack, w)
        }
      } else {
        state[v] <- 2L
        stack <- stack[-length(stack)]
      }
    }
  }
  character()
}

# Moralization and triangulation ----------------------------------------------

#' Moral graph
#'
#' Marries the parents of every node and drops edge directions.
#'
#' @param x A DAG or Bayesian network.
#' @return A `pgm_ugraph`.
#' @examples
#' moralize(dag(list(c("A", "C"), c("B", "C"))))
#' @export
moralize <- function(x) {
  g <- as_pgm_graph(x)
  if (!inherits(g, "pgm_dag")) pgm_abort("`moralize()` requires a DAG.")
  from <- g$edges$from
  to <- g$edges$to
  for (v in g$nodes) {
    pa <- parents(g, v)
    if (length(pa) > 1L) {
      cmb <- combn(pa, 2L)
      from <- c(from, cmb[1, ])
      to <- c(to, cmb[2, ])
    }
  }
  undirected_from_pairs(g$nodes, from, to)
}

undirected_from_pairs <- function(nodes, from, to) {
  if (length(from)) {
    a <- pmin(from, to)
    b <- pmax(from, to)
    keep <- !duplicated(paste(a, b, sep = "\r"))
    from <- a[keep]
    to <- b[keep]
  }
  new_pgm_graph(nodes, edge_frame(from, to, rep("undirected", length(from))), "pgm_ugraph")
}

skeleton_of <- function(g) {
  e <- g$edges
  undirected_from_pairs(g$nodes, e$from, e$to)
}

#' Skeleton of a graph
#'
#' @param x A graph or model.
#' @return The undirected graph obtained by dropping edge directions.
#' @export
skeleton <- function(x) skeleton_of(as_pgm_graph(x))

adjacency_sets <- function(g) {
  e <- g$edges
  a <- split(c(e$to, e$from), factor(c(e$from, e$to), levels = g$nodes))
  lapply(a, unique)
}

# Greedy elimination ordering on an undirected graph. `weights` gives node
# cardinalities for the weighted criteria.
greedy_order <- function(adj, heuristic = c("min_fill", "min_degree", "weighted_min_fill", "min_weight"),
                         weights = NULL, eliminate = names(adj)) {
  heuristic <- match.arg(heuristic)
  nodes <- names(adj)
  n <- length(nodes)
  w <- if (is.null(weights)) rep(2, n) else as.numeric(weights[nodes])
  A <- matrix(FALSE, n, n)
  for (i in seq_len(n)) A[i, match(adj[[i]], nodes)] <- TRUE
  alive <- rep(TRUE, n)
  remaining <- match(eliminate, nodes)
  order <- integer()
  cost_of <- function(v) {
    nb <- which(A[v, ] & alive)
    k <- length(nb)
    switch(heuristic,
      min_degree = k,
      min_weight = w[v] * prod(w[nb]),
      min_fill = if (k < 2L) 0 else (k * (k - 1) - sum(A[nb, nb])) / 2,
      weighted_min_fill = if (k < 2L) 0 else {
        wn <- w[nb]
        (sum(wn)^2 - sum(wn^2) - sum(outer(wn, wn)[A[nb, nb]])) / 2
      }
    )
  }
  while (length(remaining)) {
    cost <- vapply(remaining, cost_of, numeric(1))
    v <- remaining[which.min(cost)]
    nb <- which(A[v, ] & alive)
    A[nb, nb] <- TRUE
    diag(A) <- FALSE
    alive[v] <- FALSE
    order <- c(order, v)
    remaining <- remaining[remaining != v]
  }
  nodes[order]
}

#' Triangulate an undirected graph
#'
#' Eliminates nodes in the given or heuristic order, adding fill-in edges, and
#' returns the chordal graph together with the elimination order and the
#' maximal cliques it induces.
#'
#' @param x An undirected graph (a DAG is moralized first).
#' @param order Optional elimination order (all nodes).
#' @param heuristic Ordering heuristic used when `order` is `NULL`:
#'   `"min_fill"`, `"min_degree"`, `"weighted_min_fill"` or `"min_weight"`.
#' @param cardinality Optional named numeric vector of node cardinalities used
#'   by the weighted heuristics (default 2 for every node).
#' @return A list with `graph` (a chordal `pgm_ugraph`), `order` and `cliques`
#'   (a list of character vectors).
#' @examples
#' cycle4 <- undirected_graph(list(c("A","B"), c("B","C"), c("C","D"), c("D","A")))
#' tri <- triangulate(cycle4)
#' tri$cliques
#' @export
triangulate <- function(x, order = NULL, heuristic = "min_fill", cardinality = NULL) {
  g <- as_pgm_graph(x)
  if (inherits(g, "pgm_dag")) g <- moralize(g)
  if (any(g$edges$type != "undirected")) pgm_abort("`triangulate()` needs an undirected graph or a DAG.")
  adj <- adjacency_sets(g)
  if (is.null(order)) {
    order <- greedy_order(adj, heuristic, cardinality)
  } else {
    if (!setequal(order, g$nodes) || anyDuplicated(order)) pgm_abort("`order` must list every node exactly once.")
  }
  from <- g$edges$from
  to <- g$edges$to
  cliques <- list()
  work <- adj
  for (v in order) {
    nb <- work[[v]]
    cl <- c(v, nb)
    if (length(nb) > 1L) {
      cmb <- combn(nb, 2L)
      from <- c(from, cmb[1, ])
      to <- c(to, cmb[2, ])
    }
    for (a in nb) work[[a]] <- union(setdiff(work[[a]], v), setdiff(nb, a))
    work[[v]] <- NULL
    cliques[[length(cliques) + 1L]] <- order_nodes(g, cl)
  }
  # Keep maximal cliques only.
  keep <- vapply(seq_along(cliques), function(i) {
    !any(vapply(seq_along(cliques), function(j) {
      j != i && all(cliques[[i]] %in% cliques[[j]]) &&
        (length(cliques[[j]]) > length(cliques[[i]]) || j < i)
    }, logical(1)))
  }, logical(1))
  list(graph = undirected_from_pairs(g$nodes, from, to), order = order, cliques = cliques[keep])
}

#' Chordality test
#'
#' Uses maximum cardinality search.
#'
#' @param x An undirected graph.
#' @return Logical scalar.
#' @export
is_chordal <- function(x) {
  g <- as_pgm_graph(x)
  if (any(g$edges$type != "undirected")) pgm_abort("`is_chordal()` needs an undirected graph.")
  adj <- adjacency_sets(g)
  n <- length(g$nodes)
  if (n <= 3L) return(TRUE)
  numbered <- character()
  weight <- setNames(integer(n), g$nodes)
  for (i in seq_len(n)) {
    cand <- setdiff(g$nodes, numbered)
    v <- cand[which.max(weight[cand])]
    prev <- intersect(adj[[v]], numbered)
    # prev must form a clique
    if (length(prev) > 1L) {
      for (a in prev) if (!all(setdiff(prev, a) %in% adj[[a]])) return(FALSE)
    }
    numbered <- c(numbered, v)
    nb <- setdiff(adj[[v]], numbered)
    weight[nb] <- weight[nb] + 1L
  }
  TRUE
}

# Separation -------------------------------------------------------------------

#' Separation and d-separation
#'
#' `d_separated()` tests whether every node in `x` is d-separated from every
#' node in `y` given `z` in a DAG (Koller and Friedman, Algorithm 3.1).
#' `active_trail_nodes()` returns the nodes reachable from `x` by an active
#' trail given `z`. `separated()` tests ordinary vertex separation in an
#' undirected graph: every path from `x` to `y` passes through `z`.
#'
#' @param g A DAG, undirected graph or model.
#' @param x,y Character vectors of nodes.
#' @param z Character vector of conditioning nodes (may be empty).
#' @return `d_separated()` and `separated()` return a logical scalar;
#'   `active_trail_nodes()` returns a character vector.
#' @examples
#' g <- dag(list(c("A", "C"), c("B", "C"), c("C", "D")))
#' d_separated(g, "A", "B")
#' d_separated(g, "A", "B", "D")
#' active_trail_nodes(g, "A", "C")
#' @export
d_separated <- function(g, x, y, z = character()) {
  g <- as_pgm_graph(g)
  if (!inherits(g, "pgm_dag")) pgm_abort("`d_separated()` requires a DAG.")
  check_nodes_in(g, c(x, y, z))
  if (length(intersect(x, y))) return(FALSE)
  !any(y %in% active_trail_nodes(g, x, z))
}

#' @rdname d_separated
#' @export
active_trail_nodes <- function(g, x, z = character()) {
  g <- as_pgm_graph(g)
  if (!inherits(g, "pgm_dag")) pgm_abort("`active_trail_nodes()` requires a DAG.")
  check_nodes_in(g, c(x, z))
  e <- g$edges
  pa <- split(e$from, factor(e$to, levels = g$nodes))
  ch <- split(e$to, factor(e$from, levels = g$nodes))
  # Ancestors of z (including z)
  anc_z <- c(z, if (length(z)) reach_directed(g, z, upward = TRUE) else character())
  visited <- new.env(hash = TRUE)
  reachable <- character()
  # direction: "up" means arrived from a child, "down" from a parent
  queue <- lapply(x, function(v) c(v, "up"))
  while (length(queue)) {
    item <- queue[[1L]]
    queue <- queue[-1L]
    v <- item[1L]
    d <- item[2L]
    key <- paste(v, d)
    if (!is.null(visited[[key]])) next
    assign(key, TRUE, envir = visited)
    if (!(v %in% z)) reachable <- c(reachable, v)
    if (d == "up" && !(v %in% z)) {
      for (p in pa[[v]]) queue[[length(queue) + 1L]] <- c(p, "up")
      for (cc in ch[[v]]) queue[[length(queue) + 1L]] <- c(cc, "down")
    } else if (d == "down") {
      if (!(v %in% z)) for (cc in ch[[v]]) queue[[length(queue) + 1L]] <- c(cc, "down")
      if (v %in% anc_z) for (p in pa[[v]]) queue[[length(queue) + 1L]] <- c(p, "up")
    }
  }
  order_nodes(g, unique(reachable))
}

#' @rdname d_separated
#' @export
separated <- function(g, x, y, z = character()) {
  g <- as_pgm_graph(g)
  check_nodes_in(g, c(x, y, z))
  if (length(intersect(x, y))) return(FALSE)
  adj <- adjacency_sets(g)
  seen <- setdiff(x, z)
  frontier <- seen
  while (length(frontier)) {
    nxt <- setdiff(unique(unlist(adj[frontier], use.names = FALSE)), c(seen, z))
    if (any(nxt %in% y)) return(FALSE)
    seen <- c(seen, nxt)
    frontier <- nxt
  }
  TRUE
}

# Equivalence ------------------------------------------------------------------

#' V-structures, CPDAGs and Markov equivalence
#'
#' `v_structures()` lists unshielded colliders `a -> c <- b` with `a` and `b`
#' nonadjacent. `cpdag()` returns the completed partially directed acyclic
#' graph (essential graph) of a DAG: v-structures are oriented and Meek's
#' rules propagate compelled orientations. `markov_equivalent()` tests whether
#' two DAGs have the same skeleton and v-structures (Verma and Pearl, 1990).
#'
#' @param x,y DAGs or Bayesian networks.
#' @return `v_structures()` returns a data frame with columns `a`, `c`, `b`;
#'   `cpdag()` returns a `pgm_pdag`; `markov_equivalent()` a logical scalar.
#' @examples
#' g1 <- dag(list(c("A", "B"), c("B", "C")))
#' g2 <- dag(list(c("C", "B"), c("B", "A")))
#' markov_equivalent(g1, g2)
#' cpdag(dag(list(c("A", "C"), c("B", "C"))))
#' @export
v_structures <- function(x) {
  g <- as_pgm_graph(x)
  adj <- adjacency_sets(g)
  out <- list()
  for (cnode in g$nodes) {
    pa <- parents(g, cnode)
    if (length(pa) < 2L) next
    cmb <- combn(pa, 2L)
    for (k in seq_len(ncol(cmb))) {
      a <- cmb[1, k]
      b <- cmb[2, k]
      if (!(b %in% adj[[a]])) out[[length(out) + 1L]] <- c(a, cnode, b)
    }
  }
  if (!length(out)) return(data.frame(a = character(), c = character(), b = character(), stringsAsFactors = FALSE))
  m <- do.call(rbind, out)
  data.frame(a = m[, 1], c = m[, 2], b = m[, 3], stringsAsFactors = FALSE)
}

#' @rdname v_structures
#' @export
cpdag <- function(x) {
  g <- as_pgm_graph(x)
  if (!inherits(g, "pgm_dag")) pgm_abort("`cpdag()` requires a DAG.")
  vs <- v_structures(g)
  e <- g$edges
  compelled <- paste(e$from, e$to) %in% c(paste(vs$a, vs$c), paste(vs$b, vs$c))
  e$type <- ifelse(compelled, "directed", "undirected")
  p <- new_pgm_graph(g$nodes, e, "pgm_pdag")
  meek_rules(p)
}

# Apply Meek's orientation rules R1 to R4 until no change.
meek_rules <- function(p) {
  repeat {
    changed <- FALSE
    e <- p$edges
    und <- which(e$type == "undirected")
    if (!length(und)) break
    adj <- adjacency_sets(p)
    is_dir <- function(a, b) any(e$from == a & e$to == b & e$type == "directed")
    is_und <- function(a, b) any(e$type == "undirected" & ((e$from == a & e$to == b) | (e$from == b & e$to == a)))
    is_adj <- function(a, b) b %in% adj[[a]]
    for (i in und) {
      for (orient in list(c(e$from[i], e$to[i]), c(e$to[i], e$from[i]))) {
        a <- orient[1L]
        b <- orient[2L]
        rule <- FALSE
        # R1: c -> a, c not adjacent to b  =>  a -> b
        cs <- e$from[e$to == a & e$type == "directed"]
        if (any(!vapply(cs, function(cc) is_adj(cc, b) || cc == b, logical(1)))) rule <- TRUE
        # R2: a -> c -> b  =>  a -> b
        if (!rule) {
          mids <- e$to[e$from == a & e$type == "directed"]
          if (any(vapply(mids, function(cc) is_dir(cc, b), logical(1)))) rule <- TRUE
        }
        # R3: a - c, a - d, c -> b, d -> b, c and d nonadjacent  =>  a -> b
        if (!rule) {
          und_nb <- adj[[a]][vapply(adj[[a]], function(cc) is_und(a, cc), logical(1))]
          into_b <- und_nb[vapply(und_nb, function(cc) is_dir(cc, b), logical(1))]
          if (length(into_b) >= 2L) {
            cmb <- combn(into_b, 2L)
            if (any(apply(cmb, 2L, function(cd) !is_adj(cd[1], cd[2])))) rule <- TRUE
          }
        }
        # R4: a - c, c -> d, d -> b, a adjacent to d (d and c ...), c not adjacent to b  =>  a -> b
        if (!rule) {
          und_nb <- adj[[a]][vapply(adj[[a]], function(cc) is_und(a, cc), logical(1))]
          for (cc in und_nb) {
            ds <- e$to[e$from == cc & e$type == "directed"]
            for (d in ds) {
              if (is_dir(d, b) && is_adj(a, d) && !is_adj(cc, b)) {
                rule <- TRUE
                break
              }
            }
            if (rule) break
          }
        }
        if (rule) {
          e$from[i] <- a
          e$to[i] <- b
          e$type[i] <- "directed"
          p$edges <- e
          changed <- TRUE
          break
        }
      }
      if (changed) break
    }
    if (!changed) break
  }
  p$edges <- e
  if (all(p$edges$type == "directed") && !has_directed_cycle(p)) class(p) <- c("pgm_dag", "pgm_pdag", "pgm_graph")
  else class(p) <- c("pgm_pdag", "pgm_graph")
  p
}

#' @rdname v_structures
#' @export
markov_equivalent <- function(x, y) {
  gx <- as_pgm_graph(x)
  gy <- as_pgm_graph(y)
  if (!setequal(gx$nodes, gy$nodes)) return(FALSE)
  sk <- function(g) sort(paste(pmin(g$edges$from, g$edges$to), pmax(g$edges$from, g$edges$to)))
  if (!identical(sk(gx), sk(gy))) return(FALSE)
  vk <- function(g) {
    v <- v_structures(g)
    sort(paste(pmin(v$a, v$b), v$c, pmax(v$a, v$b)))
  }
  identical(vk(gx), vk(gy))
}

# Comparison -------------------------------------------------------------------

#' Compare two graphs
#'
#' Counts edges that agree, are reversed, missing or extra in `estimated`
#' relative to `true`, and returns the structural Hamming distance (SHD:
#' missing plus extra plus wrongly oriented edges, where an undirected edge
#' matching a directed one counts as one orientation error).
#'
#' @param true,estimated Graphs or models over the same nodes.
#' @return A list of class `pgm_graph_comparison` with elements `tp`,
#'   `reversed`, `fn`, `fp`, `shd`, `precision`, `recall`.
#' @examples
#' compare_graphs(dag(list(c("A", "B"), c("B", "C"))), dag(list(c("B", "A"), c("B", "C"))))
#' @export
compare_graphs <- function(true, estimated) {
  gt <- as_pgm_graph(true)
  ge <- as_pgm_graph(estimated)
  if (!setequal(gt$nodes, ge$nodes)) pgm_abort("The graphs must have the same nodes.")
  pk <- function(e) paste(pmin(e$from, e$to), pmax(e$from, e$to), sep = "\r")
  kt <- pk(gt$edges)
  ke <- pk(ge$edges)
  common <- intersect(kt, ke)
  same <- 0L
  wrong <- 0L
  for (k in common) {
    it <- which(kt == k)
    ie <- which(ke == k)
    a <- gt$edges[it, ]
    b <- ge$edges[ie, ]
    agree <- a$type == b$type && (a$type == "undirected" || (a$from == b$from && a$to == b$to))
    if (agree) same <- same + 1L else wrong <- wrong + 1L
  }
  fn <- length(setdiff(kt, ke))
  fp <- length(setdiff(ke, kt))
  structure(list(
    tp = same, reversed = wrong, fn = fn, fp = fp, shd = fn + fp + wrong,
    precision = if (length(ke)) same / length(ke) else NA_real_,
    recall = if (length(kt)) same / length(kt) else NA_real_
  ), class = "pgm_graph_comparison")
}

#' @export
print.pgm_graph_comparison <- function(x, ...) {
  cat(sprintf("<Graph comparison> SHD = %d; correct = %d, misoriented = %d, missing = %d, extra = %d\n",
              x$shd, x$tp, x$reversed, x$fn, x$fp))
  invisible(x)
}

#' Graph isomorphism
#'
#' Tests whether two graphs are isomorphic, respecting edge directions. Uses
#' a degree refined backtracking search, suitable for graphs with up to a few
#' dozen nodes.
#'
#' @param x,y Graphs or models.
#' @return Logical scalar with attribute `mapping` (a named character vector
#'   from nodes of `x` to nodes of `y`) when `TRUE`.
#' @examples
#' is_isomorphic(dag(list(c("A", "B"))), dag(list(c("X", "Y"))))
#' @export
is_isomorphic <- function(x, y) {
  gx <- as_pgm_graph(x)
  gy <- as_pgm_graph(y)
  nx <- gx$nodes
  ny <- gy$nodes
  if (length(nx) != length(ny) || nrow(gx$edges) != nrow(gy$edges)) return(FALSE)
  if (!identical(sort(table(gx$edges$type)), sort(table(gy$edges$type)))) return(FALSE)
  am <- function(g) {
    m <- matrix(0L, length(g$nodes), length(g$nodes), dimnames = list(g$nodes, g$nodes))
    e <- g$edges
    for (i in seq_len(nrow(e))) {
      if (e$type[i] == "directed") m[e$from[i], e$to[i]] <- 1L
      else {
        m[e$from[i], e$to[i]] <- 2L
        m[e$to[i], e$from[i]] <- 2L
      }
    }
    m
  }
  ax <- am(gx)
  ay <- am(gy)
  sig <- function(m) paste(rowSums(m == 1L), colSums(m == 1L), rowSums(m == 2L))
  sx <- setNames(sig(ax), nx)
  sy <- setNames(sig(ay), ny)
  if (!identical(sort(unname(sx)), sort(unname(sy)))) return(FALSE)
  map <- setNames(character(length(nx)), nx)
  used <- setNames(logical(length(ny)), ny)
  ord <- nx[order(-rowSums(ax > 0) - colSums(ax > 0))]
  search <- function(k) {
    if (k > length(ord)) return(TRUE)
    v <- ord[k]
    for (w in ny[!used & sy == sx[v]]) {
      ok <- TRUE
      if (k > 1L) {
        prev <- ord[seq_len(k - 1L)]
        if (any(ax[v, prev] != ay[w, map[prev]]) || any(ax[prev, v] != ay[map[prev], w])) ok <- FALSE
      }
      if (ok) {
        map[v] <<- w
        used[w] <<- TRUE
        if (search(k + 1L)) return(TRUE)
        used[w] <<- FALSE
        map[v] <<- ""
      }
    }
    FALSE
  }
  res <- search(1L)
  if (res) structure(TRUE, mapping = map) else FALSE
}

#' Adjacency matrix of a graph
#'
#' @param x A graph or model.
#' @return An integer matrix with `m[i, j] = 1` for a directed edge `i -> j`
#'   and `m[i, j] = m[j, i] = 1` for an undirected edge.
#' @export
adjacency_matrix <- function(x) {
  g <- as_pgm_graph(x)
  m <- matrix(0L, length(g$nodes), length(g$nodes), dimnames = list(g$nodes, g$nodes))
  e <- g$edges
  if (nrow(e)) {
    m[cbind(e$from, e$to)] <- 1L
    u <- e$type == "undirected"
    if (any(u)) m[cbind(e$to[u], e$from[u])] <- 1L
  }
  m
}

#' Build a graph from an adjacency matrix
#'
#' Symmetric entries become undirected edges; asymmetric entries become
#' directed edges.
#'
#' @param m A square matrix with row and column names.
#' @return A `pgm_graph` (classified as DAG, undirected or PDAG when possible).
#' @export
graph_from_adjacency <- function(m) {
  if (!is.matrix(m) || nrow(m) != ncol(m)) pgm_abort("`m` must be a square matrix.")
  nm <- rownames(m) %||% colnames(m) %||% paste0("V", seq_len(nrow(m)))
  m <- m != 0
  from <- to <- type <- character()
  for (i in seq_len(nrow(m))) for (j in seq_len(ncol(m))) {
    if (!m[i, j] || i == j) next
    if (m[j, i]) {
      if (i < j) {
        from <- c(from, nm[i]); to <- c(to, nm[j]); type <- c(type, "undirected")
      }
    } else {
      from <- c(from, nm[i]); to <- c(to, nm[j]); type <- c(type, "directed")
    }
  }
  pgm_graph(edge_frame(from, to, type), nodes = nm)
}
