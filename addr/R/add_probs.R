#' Add Probabilities to Data Frames
#'
#' This is a generic function to append event probabilities to a data
#' frame. An event probability, such as Pr(Response < 10), is
#' generated for the fitted value of each observation in
#' \code{tb}. These probabilities are then appended to \code{tb} and
#' returned to the user as a tibble.
#'
#' @import lme4
#' @import tidyverse
#' @import merTools 
#' @import arm
#'
#' @param tb A tibble or Data Frame on which to append probabilities
#' @param fit An object of class lm, glm, or lmerMod. Predictions are
#'     made with this object.
#' @param quant A double. A quantile of the response distribution.
#' @param probName NULL or character vector of length one. If
#'     \code{NULL}, probabilities will automatically be named by
#'     \code{add_probs()}, otherwise, the probabilities will be named
#'     \code{probName} in the returned tibble
#' @param comparison Default is "<". Must be "<" or ">" for
#'     linear, log-linear and linear mixed models. If \code{fit} is a
#'     glm, then comparison may also be "<=", ">=", or "=".
#' @param ... Additional arguments
#' @return A tibble, \code{tb}, with predicted values and
#'     probabilities attached.
#' 
#' @examples
#' # linear regression
#' fit1 <- lm(dist ~ speed, data = cars)
#' # Calculate Pr(dist < 40) for each observation in cars
#' add_probs(cars, fit1, quant = 40)
#' # Poisson regression
#' fit2 <- glm(dist ~ speed, data = cars, family = "poisson")
#' # Calculate Pr(dist >= 40) for each observation in cars
#' add_probs(cars, fit2, quant = 40, comparison = ">=")
#' 
#' @export
add_probs <- function(tb, fit, quant, probName = NULL, comparison, ...){
    UseMethod("add_probs", fit)
}


