#' Add Confidence Intervals for Predictions to Data Frames.
#'
#' This is a generic function to append confidence intervals for
#' predictions of a model fit to a data frame. A confidence interval
#' is generated for the fitted value of each observation in
#' \code{tb}. These confidence intervals are then appended to
#' \code{tb} and returned to the user as a tibble.
#'
#' For more specific information about the arguments that are useful
#' in each method, consult
#'
#' \itemize{
#'   \item \code{\link{add_ci.lm}} for linear regression confidence intervals
#'   \item \code{\link{add_ci.glm}} for generalized linear regression confidence intervals
#'   \item \code{\link{add_ci.lmerMod}} for linear mixed models confidence intervals
#' }
#'
#' @import stats
#' @import lme4
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @importFrom tibble as_data_frame
#' @importFrom merTools predictInterval
#' @importFrom arm se.ranef
#' @importFrom arm sim
#'
#' @param tb A tibble or data frame on which to make predictions.
#' @param fit An object of class lm, glm, or lmerMod. Predictions are
#'     made with this object.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names NULL or character vector of length two. If
#'     \code{NULL}, confidence bounds will automatically be named by
#'     \code{add_ci}, otherwise, the lower confidence bound will be
#'     named \code{ciNames[1]} and the upper confidence bound will be
#'     named \code{ciNames[2]}.
#' @param ... Additional arguments. 
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     confidence bounds attached.
#'
#' @examples
#' fit <- lm(dist ~ speed, data = cars)
#' add_ci(cars, fit)
#' fit2 <- glm(dist ~ speed, family = "poisson", data = cars) 
#' add_ci(cars, fit2)
#' fit3 <- lme4::lmer(Reaction ~ Days + (1|Subject), data = lme4::sleepstudy)
#' add_ci(lme4::sleepstudy, fit3)
#' fit4 <- glm(I(dist > 20) ~ speed, family = "binomial", data = cars)
#' add_ci(cbind(cars, I(cars$dist > 20)), fit4)
#'
#'
#' @seealso \code{\link{add_pi}} for prediction intervals,
#'     \code{\link{add_probs}} for response level probabilities, and
#'     \code{\link{add_quantile}} for quantiles of the conditional
#'     response distribution.
#'
#' @export

add_ci <- function(tb, fit, alpha = 0.05, names = NULL, ...){
  UseMethod("add_ci", fit)
}
