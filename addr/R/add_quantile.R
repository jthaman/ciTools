#' Add Regression Quantiles to Data Frames
#'
#' This is a generic function to append regression quantiles to a data
#' frame. A regression quantile \emph{q} is a point such that
#' Pr(Response | Covariates < q) = p. These quantiles are generated
#' for the fitted value of each observation in \code{tb}. Quantiles
#' are then appended to \code{tb} and returned to the user as a
#' tibble.
#'
#' @import lme4
#' @import tidyverse
#' @import merTools 
#' @import arm
#'
#' @param tb A tibble or Data Frame on which to append probabilities
#' @param fit An object of class lm, glm, or lmerMod. Predictions are
#'     made with this object.
#' @param p A double. A probability that determines the
#'     quantile. Must be between 0 and 1.
#' @param name NULL or character vector of length one. If \code{NULL},
#'     quantiles will automatically be named by \code{add_quantile()},
#'     otherwise, the quantiles will be named \code{quantileName} in
#'     the returned tibble
#' @param ... Additional arguments
#' @return A tibble, \code{tb}, with predicted values and quantiles
#'     attached.
#' 
#' @examples
#' # linear regression
#' fit1 <- lm(dist ~ speed, data = cars)
#' # Calculate the 0.2 quantile of dist
#' add_quantile(cars, fit1, prob = 0.2)
#' # Poisson regression
#' fit2 <- glm(dist ~ speed, data = cars, family = "poisson")
#' # Calculate the median of the response variable, dist.
#' add_quantile(cars, fit2, prob = 0.5)
#' 
#' @export

add_quantile <- function(tb, fit, p, name = NULL, ...){
  UseMethod("add_quantile", fit)
}

