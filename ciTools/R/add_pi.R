#' Add Prediction Intervals to Data Frames
#'
#' This is a generic function to append prediction intervals to a data
#' frame. A prediction interval is made for each observation in
#' \code{tb} with respect to the model \code{fit}. These intervals are
#' then appended to \code{tb} and returned to the user as a
#' tibble. \code{fit} can be a linear model, log-linear mode, linear
#' mixed model, or Poisson model.
#'
#' @param tb A tibble or Data Frame on which to make predictions.
#' @param fit An object of class lm, glm, or lmerMod. Predictions are
#'     made with this object.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names NULL or character vector of length two. If
#'     \code{NULL}, prediction bounds will automatically be named by
#'     \code{add_pi}, otherwise, the lower prediction bound will be
#'     named \code{piNames[1]} and the upper prediction bound will be
#'     named \code{piNames[2]}.
#' @param ... Additional arguments
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#'
#' @examples
#' ## linear regression
#' fit1 <- lm(dist ~ speed, data = cars)
#' add_pi(cars, fit1, alpha = 0.5)
#' ## Poisson regression
#' fit2 <- glm(dist ~ speed, data = cars, family = "poisson")
#' add_pi(cars, fit2)
#' 
#' @export

add_pi <- function(tb, fit, alpha = 0.05, names = NULL, ...){
  UseMethod("add_pi", fit)
}


