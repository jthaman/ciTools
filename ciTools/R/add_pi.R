# Copyright (C) 2017 Institute for Defense Analyses
#
# This file is part of ciTools.
#
# ciTools is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ciTools is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ciTools. If not, see <http://www.gnu.org/licenses/>.

#' Add Prediction Intervals to Data Frames
#'
#' This is a generic function to append prediction intervals to a data
#' frame. A prediction interval is made for each observation in
#' \code{tb} with respect to the model \code{fit}. These intervals are
#' then appended to \code{tb} and returned to the user as a
#' tibble. \code{fit} can be a linear, log-linear, linear
#' mixed, or generalized linear models.
#' 
#' For more specific information about the arguments that are useful
#' in each method, consult
#'
#' \itemize{
#'   \item \code{\link{add_pi.lm}} for linear regression prediction intervals
#'   \item \code{\link{add_pi.glm}} for generalized linear regression prediction intervals
#'   \item \code{\link{add_pi.lmerMod}} for linear mixed models prediction intervals
#' }
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
#' fit <- lm(dist ~ speed, data = cars)
#' add_pi(cars, fit)
#' fit2 <- glm(dist ~ speed, family = "poisson", data = cars) 
#' add_pi(cars, fit2)
#' fit3 <- lme4::lmer(Reaction ~ Days + (1|Subject), data = lme4::sleepstudy)
#' add_pi(lme4::sleepstudy, fit3)
#'
#'
#' @seealso \code{\link{add_ci}} for confidence intervals,
#'     \code{\link{add_probs}} for response level probabilities, and
#'     \code{\link{add_quantile}} for quantiles of the conditional
#'     response distribution.
#' 
#' @export

add_pi <- function(tb, fit, alpha = 0.05, names = NULL, ...){
  UseMethod("add_pi", fit)
}


