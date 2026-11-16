# pgmcompR 0.1.0

First release.

## Representation

* Directed, undirected and partially directed graphs with d-separation,
  moralization, triangulation, CPDAGs, Meek rules and graph comparison.
* Discrete factors with explicit state labels, linear and log space
  algebra, and an R and a C++ backend that agree to 1e-12.
* Tabular, noisy-OR, deterministic, custom, linear Gaussian and
  conditional Gaussian CPDs, plus an S3 contract for user defined CPDs.
* Bayesian networks with a variable registry, Markov networks, factor
  graphs, junction trees, dynamic Bayesian networks and hidden Markov
  models.

## Inference

* Enumeration, variable elimination with four greedy orderings, junction
  tree calibration, sum product and max product belief propagation, loopy
  belief propagation with damping.
* MAP and MPE queries, soft (virtual) evidence, interventions, and an error
  class for evidence of probability zero.
* Rejection sampling, likelihood weighting, importance sampling and Gibbs
  sampling with split R-hat, bulk ESS and MCSE.
* Exact Gaussian and conditional Gaussian inference, Gaussian belief
  propagation, filtering, smoothing, Viterbi decoding, forecasting,
  particle filtering and the Kalman filter and smoother.

## Learning

* Maximum likelihood, Bayesian (BDeu, K2, custom pseudo counts) and EM
  parameter estimation, with explicit policies for unobserved parent
  configurations and structural zeros.
* PC (stable), hill climbing, tabu search and MMHC with blacklists,
  whitelists and tiers; BIC, AIC, BDeu, K2 and log likelihood scores;
  bootstrap edge strengths.

## Causal inference

* `do_query()` by truncated factorization, back door and front door
  criteria, adjustment set search, effect identification and
  counterfactuals in discrete structural causal models.

## Interoperability

* BIF, XMLBIF, UAI and JSON input and output; conversion to and from
  bnlearn, gRain and igraph.
