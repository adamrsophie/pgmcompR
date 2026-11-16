# pgmcompR <img src="man/figures/logo.png" align="right" height="139" alt="pgmcompR logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/adamrsophie/pgmcompR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/adamrsophie/pgmcompR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

pgmcompR is a framework for probabilistic graphical models written in R.
It covers representation, exact and approximate inference, parameter and
structure learning, dynamic models and causal queries for discrete,
Gaussian and conditional Gaussian models, with one consistent interface.

Every algorithm has a pure R implementation that can be read and
modified; an optional C++ backend for factor algebra is checked against it
to 1e-12. Exact algorithms are tested against an independent brute force
enumeration, and approximate results carry diagnostics and their own class
so they are never mistaken for exact ones.

## Installation

```r
# install.packages("remotes")
remotes::install_github("adamrsophie/pgmcompR")
```

## Example

```r
library(pgmcompR)

bn <- bayesian_network(list(c("Smoking", "Cancer"), c("Cancer", "XRay")))
bn <- add_cpd(bn, "Smoking", c(no = 0.7, yes = 0.3))
bn <- add_cpd(bn, "Cancer", matrix(c(0.99, 0.01, 0.90, 0.10), 2,
                                   dimnames = list(c("no", "yes"), NULL)))
bn <- add_cpd(bn, "XRay", matrix(c(0.95, 0.05, 0.20, 0.80), 2,
                                 dimnames = list(c("negative", "positive"), NULL)))

query(bn, "Cancer", evidence = list(XRay = "positive"))
#> <Factor: Cancer>
#>  Cancer     value
#>      no 0.7795455
#>     yes 0.2204545

do_query(bn, "Smoking", interventions = list(XRay = "positive"))
#> <Factor: Smoking>
#>  Smoking value
#>       no   0.7
#>      yes   0.3
```

## What is included

| Area | Functions |
|---|---|
| Graphs | `dag()`, `undirected_graph()`, `pdag()`, `d_separated()`, `moralize()`, `triangulate()`, `cpdag()`, `compare_graphs()` |
| Models | `bayesian_network()`, `markov_network()`, `factor_graph()`, `junction_tree()`, `dynamic_bn()`, `hmm()`, `scm()` |
| Distributions | `cpd()`, `noisy_or_cpd()`, `deterministic_cpd()`, `lg_cpd()`, `cg_cpd()`, custom CPD classes |
| Exact inference | `query()` with variable elimination, junction tree or enumeration; `map_query()`, `mpe_query()`, `belief_propagation()` |
| Approximate inference | rejection sampling, likelihood weighting, importance sampling, Gibbs sampling, `loopy_bp()` |
| Gaussian models | `gaussian_joint()`, `gaussian_query()`, `gaussian_bp()`, `kalman_filter()`, `kalman_smoother()` |
| Dynamic models | `filter_dbn()`, `smooth_dbn()`, `viterbi_dbn()`, `forecast_dbn()`, `particle_filter()` |
| Learning | `fit_parameters()` (MLE, Bayesian, EM), `pc_algorithm()`, `hill_climb()`, `tabu_search()`, `mmhc()`, `bootstrap_structure()` |
| Causal | `do_query()`, `is_backdoor_set()`, `adjustment_sets()`, `identify_effect()`, `counterfactual()` |
| Exchange | BIF, XMLBIF, UAI, JSON; `as_bnlearn()`, `from_bnlearn()`, `as_grain()`, `from_grain()`, `as_igraph()` |

The vignettes walk through each area, from a first Bayesian network to
causal inference and extension with custom methods.

## Related packages

[bnlearn](https://www.bnlearn.com/) is the reference implementation of
structure learning in R and is much faster for it;
[gRain](https://CRAN.R-project.org/package=gRain) provides fast junction
tree inference. pgmcompR converts models to and from both, and is inspired
by the scope of the Python library [pgmpy](https://pgmpy.org/).

## License

GPL-3
