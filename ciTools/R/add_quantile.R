#' Add Regression Quantiles to Data Frames
#'
#' This is a generic function to append regression quantiles to a data
#' frame. A regression quantile \emph{q} is a point such that
#' Pr(Response | Covariates < q) = p. These quantiles are generated
#' for the fitted value of each observation in \code{tb}. Quantiles
#' are then appended to \code{tb} and returned to the user as a
#' tibble.
#'
#' For more specific information about the arguments that are useful
#' in each method, consult
#'
#' \itemize{
#'   \item \code{\link{add_quantile.lm}} for linear regression response quantiles
#'   \item \code{\link{add_quantile.glm}} for generalized linear regression response quantiles
#'   \item \code{\link{add_quantile.lmerMod}} for linear mixed models response quantiles
#' }
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
#' fit <- lm(dist ~ speed, data = cars)
#' add_quantile(cars, fit, p = 0.4)
#' fit2 <- glm(dist ~ speed, family = "poisson", data = cars) 
#' add_quantile(cars, fit2, p = 0.4)
#' fit3 <- lme4::lmer(Reaction ~ Days + (1|Subject), data = lme4::sleepstudy)
#' add_quantile(lme4::sleepstudy, fit3, p = 0.4)
#'
#'
#' @seealso \code{\link{add_ci}} for confidence intervals,
#'     \code{\link{add_probs}} for response level probabilities, and
#'     \code{\link{add_pi}} for prediction intervals
#' 
#' @export

add_quantile <- function(tb, fit, p, name = NULL, ...){
  UseMethod("add_quantile", fit)
}

