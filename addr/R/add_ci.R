#' Add CI generic function
#'
#' This is a generic function to append confidence intervals for the predicted
#' values of a model fit to a data frame.
#'
#' @param tb A tibble or Data Frame
#' @param fit An object or class lm, glm, or lmerMod
#' @param alpha A real number between 0 and 1
#' @param ciNames NULL or character vector of length two
#' @return A data frame, tb, with predicted values, upper and lower confidence
#'   bounds

add_ci <- function(tb, fit, alpha = 0.05, ciNames = NULL, ...){
  UseMethod("add_ci", fit)
}
