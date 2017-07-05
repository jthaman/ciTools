#' Add Confidence Intervals for Predictions to Data Frames.
#'
#' This is a generic function to append confidence intervals for
#' predictions of a model fit to a data frame. A confidence interval
#' is generated for the fitted value of each observation in
#' \code{tb}. These confidence intervals are then appended to
#' \code{tb} and returned to the user as a tibble.
#'
#'
#' @import lme4
#' @import tidyverse
#' @import merTools 
#' @import arm
#'
#' @param tb A tibble or Data Frame on which to make predictions.
#' @param fit An object of class lm, glm, or lmerMod. Predictions are
#'     made with this object.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param ciNames NULL or character vector of length two. If
#'     \code{NULL}, confidence bounds will automatically be named by
#'     \code{add_ci}, otherwise, the lower confidence bound will be
#'     named \code{ciNames[1]} and the upper confidence bound will be
#'     named \code{ciNames[2]}.
#' @param ... Additional arguments
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     confidence bounds attached.
#'
#' @examples
#' # linear regression
#' # Append a 50% confidence interval for the expected response to
#' #cars
#' fit1 <- lm(dist ~ speed, data = cars)
#' add_ci(cars, fit1, alpha = 0.5)
#' 
#' # Poisson regression
#' # Append a 95% confidence interval for the expected response to
#' #cars
#' fit2 <- glm(dist ~ speed, data = cars, family = "poisson")
#' add_ci(cars, fit2)
#' 
#' @export
add_ci <- function(tb, fit, alpha = 0.05, ciNames = NULL, ...){
  UseMethod("add_ci", fit)
}
