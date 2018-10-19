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
#' tibble. \code{fit} can be a linear, log-linear, linear mixed,
#' generalized linear, generalized linear mixed, or accelerated
#' failure time model.
#'
#' For more specific information about the arguments that are applicable
#' in each method, consult:
#'
#' \itemize{
#'   \item \code{\link{add_pi.lm}} for linear regression prediction intervals
#'   \item \code{\link{add_pi.glm}} for generalized linear regression prediction intervals
#'   \item \code{\link{add_pi.lmerMod}} for linear mixed models prediction intervals
#'   \item \code{\link{add_pi.glmerMod}} for generalized linear mixed model prediction intervals
#'   \item \code{\link{add_pi.lmerMod}} for accelerated failure time model prediction intervals
#' }
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{lm}, \code{glm}, or
#'     \code{lmerMod}. Predictions are made with this object.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, prediction bounds automatically will be named by
#'     \code{add_pi}, otherwise, the lower prediction bound will be
#'     named \code{piNames[1]} and the upper prediction bound will be
#'     named \code{piNames[2]}.
#' @param yhatName A string. Name of the predictions vector.
#' @param ... Additional arguments
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#'
#' @examples
#' # Fit a linear model
#' fit <- lm(dist ~ speed, data = cars)
#' # Define some new data
#' new_data <- cars[sample(NROW(cars), 10), ]
#' # Add fitted values and prediction intervals to new_data
#' add_pi(new_data, fit)
#'
#' # Fit a Poisson model
#' fit2 <- glm(dist ~ speed, family = "poisson", data = cars)
#' # Add approximate prediction intervals to the fitted values of
#' # new_data
#' add_pi(new_data, fit2)
#'
#' # Fit a linear mixed model
#' fit3 <- lme4::lmer(Reaction ~ Days + (1|Subject), data = lme4::sleepstudy)
#' # Add parametric prediction intervals for the fitted values to the
#' # original data
#' add_pi(lme4::sleepstudy, fit3, type = "parametric")
#'
#'
#' @seealso \code{\link{add_ci}} for confidence intervals,
#'     \code{\link{add_probs}} for response level probabilities, and
#'     \code{\link{add_quantile}} for quantiles of the conditional
#'     response distribution.
#'
#' @export

add_pi <- function(tb, fit, alpha = 0.05, names = NULL, yhatName = "pred", ...){
  UseMethod("add_pi", fit)
}
