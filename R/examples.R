# Example networks --------------------------------------------------------------

#' Example networks
#'
#' Small, well known networks used in examples, tests and vignettes.
#'
#' * `"asia"`: the chest clinic network of Lauritzen and Spiegelhalter (1988)
#'   with 8 binary nodes (`asia`, `tub`, `smoke`, `lung`, `bronc`, `either`,
#'   `xray`, `dysp`), states `"no"`/`"yes"`.
#' * `"sprinkler"`: the four node rain and sprinkler network (Murphy's
#'   parameterization) with nodes `Cloudy`, `Sprinkler`, `Rain`, `WetGrass`,
#'   states `"F"`/`"T"`.
#' * `"cancer"`: the five node cancer network of Korb and Nicholson (2010):
#'   `Pollution`, `Smoker`, `Cancer`, `Xray`, `Dyspnoea`.
#' * `"student"`: the student network of Koller and Friedman (2009) with
#'   `Difficulty`, `Intelligence`, `Grade` (three states), `SAT`, `Letter`.
#' * `"gaussian"`: a linear Gaussian chain `X1 -> X2 -> X3` with `X1 -> X3`.
#'
#' @param name Network name.
#' @return A `bayesian_network`.
#' @references
#' Lauritzen, S. L. and Spiegelhalter, D. J. (1988). Local computations with
#' probabilities on graphical structures and their application to expert
#' systems. *Journal of the Royal Statistical Society B*, 50(2), 157 to 224.
#'
#' Koller, D. and Friedman, N. (2009). *Probabilistic Graphical Models:
#' Principles and Techniques*. MIT Press.
#' @examples
#' example_network("asia")
#' @export
example_network <- function(name = c("asia", "sprinkler", "cancer", "student", "gaussian")) {
  name <- match.arg(name)
  yn <- c("no", "yes")
  switch(name,
    asia = {
      st <- list(asia = yn, tub = yn, smoke = yn, lung = yn, bronc = yn, either = yn, xray = yn, dysp = yn)
      bn <- bayesian_network(list(c("asia", "tub"), c("smoke", "lung"), c("smoke", "bronc"),
                                  c("tub", "either"), c("lung", "either"), c("either", "xray"),
                                  c("either", "dysp"), c("bronc", "dysp")),
                             nodes = c("asia", "tub", "smoke", "lung", "bronc", "either", "xray", "dysp"))
      bn <- add_cpd(bn, cpd("asia", c(no = 0.99, yes = 0.01)))
      bn <- add_cpd(bn, cpd("smoke", c(no = 0.5, yes = 0.5)))
      bn <- add_cpd(bn, cpd("tub", c(0.99, 0.01, 0.95, 0.05), "asia", st))
      bn <- add_cpd(bn, cpd("lung", c(0.99, 0.01, 0.9, 0.1), "smoke", st))
      bn <- add_cpd(bn, cpd("bronc", c(0.7, 0.3, 0.4, 0.6), "smoke", st))
      # either is the deterministic OR of tub and lung
      bn <- add_cpd(bn, cpd("either", c(1, 0, 0, 1, 0, 1, 0, 1), c("tub", "lung"), st))
      bn <- add_cpd(bn, cpd("xray", c(0.95, 0.05, 0.02, 0.98), "either", st))
      bn <- add_cpd(bn, cpd("dysp", c(0.9, 0.1, 0.3, 0.7, 0.2, 0.8, 0.1, 0.9), c("either", "bronc"), st))
      bn
    },
    sprinkler = {
      tf <- c("F", "T")
      st <- list(Cloudy = tf, Sprinkler = tf, Rain = tf, WetGrass = tf)
      bn <- bayesian_network(list(c("Cloudy", "Sprinkler"), c("Cloudy", "Rain"),
                                  c("Sprinkler", "WetGrass"), c("Rain", "WetGrass")))
      bn <- add_cpd(bn, cpd("Cloudy", c(F = 0.5, T = 0.5)))
      bn <- add_cpd(bn, cpd("Sprinkler", c(0.5, 0.5, 0.9, 0.1), "Cloudy", st))
      bn <- add_cpd(bn, cpd("Rain", c(0.8, 0.2, 0.2, 0.8), "Cloudy", st))
      bn <- add_cpd(bn, cpd("WetGrass", c(1, 0, 0.1, 0.9, 0.1, 0.9, 0.01, 0.99), c("Sprinkler", "Rain"), st))
      bn
    },
    cancer = {
      st <- list(Pollution = c("low", "high"), Smoker = c("False", "True"), Cancer = c("False", "True"),
                 Xray = c("negative", "positive"), Dyspnoea = c("False", "True"))
      bn <- bayesian_network(list(c("Pollution", "Cancer"), c("Smoker", "Cancer"),
                                  c("Cancer", "Xray"), c("Cancer", "Dyspnoea")))
      bn <- add_cpd(bn, cpd("Pollution", c(low = 0.9, high = 0.1)))
      bn <- add_cpd(bn, cpd("Smoker", c(False = 0.7, True = 0.3)))
      bn <- add_cpd(bn, cpd("Cancer", c(0.999, 0.001, 0.98, 0.02, 0.97, 0.03, 0.95, 0.05),
                            c("Pollution", "Smoker"), st))
      bn <- add_cpd(bn, cpd("Xray", c(0.8, 0.2, 0.1, 0.9), "Cancer", st))
      bn <- add_cpd(bn, cpd("Dyspnoea", c(0.7, 0.3, 0.35, 0.65), "Cancer", st))
      bn
    },
    student = {
      st <- list(Difficulty = c("easy", "hard"), Intelligence = c("low", "high"),
                 Grade = c("A", "B", "C"), SAT = c("low", "high"), Letter = c("weak", "strong"))
      bn <- bayesian_network(list(c("Difficulty", "Grade"), c("Intelligence", "Grade"),
                                  c("Intelligence", "SAT"), c("Grade", "Letter")))
      bn <- add_cpd(bn, cpd("Difficulty", c(easy = 0.6, hard = 0.4)))
      bn <- add_cpd(bn, cpd("Intelligence", c(low = 0.7, high = 0.3)))
      bn <- add_cpd(bn, cpd("Grade", c(0.3, 0.4, 0.3, 0.05, 0.25, 0.7, 0.9, 0.08, 0.02, 0.5, 0.3, 0.2),
                            c("Difficulty", "Intelligence"), st))
      bn <- add_cpd(bn, cpd("SAT", c(0.95, 0.05, 0.2, 0.8), "Intelligence", st))
      bn <- add_cpd(bn, cpd("Letter", c(0.1, 0.9, 0.4, 0.6, 0.99, 0.01), "Grade", st))
      bn
    },
    gaussian = {
      bn <- bayesian_network(list(c("X1", "X2"), c("X2", "X3"), c("X1", "X3")))
      bn <- add_cpd(bn, lg_cpd("X1", intercept = 1, sd = 1))
      bn <- add_cpd(bn, lg_cpd("X2", intercept = -0.5, coefficients = c(X1 = 0.8), sd = 0.6))
      bn <- add_cpd(bn, lg_cpd("X3", intercept = 0.2, coefficients = c(X2 = 1.5, X1 = -0.4), sd = 0.5))
      bn
    }
  )
}
