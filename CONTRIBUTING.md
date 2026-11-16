# Contributing to pgmcompR

Thank you for considering a contribution.

## Reporting problems

Open an issue at <https://github.com/adamrsophie/pgmcompR/issues> with a
minimal reproducible example: the model (or code that builds it), the
query and the result you expected. For numerical disagreements, please
include the answer from an independent implementation (for example
gRain, bnlearn or a brute force enumeration) if you have one.

## Development

* Code is documented with roxygen2; run `devtools::document()` after
  changing documentation.
* Tests use testthat (edition 3). Exact algorithms are tested against an
  independent brute force oracle in `tests/testthat/helper-oracle.R`;
  please add a comparison of the same kind for new exact methods.
* Every algorithm must have a pure R implementation. Compiled code is an
  optional backend and must agree with the R implementation to 1e-12.
* Do not change results silently: invalid inputs, impossible evidence and
  unsupported features raise classed errors.
* Run `devtools::check()` before opening a pull request.

## Extensions

New inference methods, scores, estimators and samplers can often live in
your own package through the registries (`?register_inference_method`)
and the CPD contract (`?cpd_contract`), without changes here.
