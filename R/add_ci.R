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

#' Add Confidence Intervals for Fitted Values to Data Frames
#'
#' This is a generic function to append confidence intervals for
#' predictions of a model fit to a data frame. A confidence interval
#' is generated for the fitted value of each observation in
#' \code{tb}. These confidence intervals are then appended to
#' \code{tb} and returned to the user as a tibble. The \code{fit} may
#' be a linear, log-linear, linear mixed, or generalized linear model.
#'
#' For more specific information about the arguments that are applicable
#' in each method, consult:
#'
#' \itemize{
#'   \item \code{\link{add_ci.lm}} for linear regression confidence intervals
#'   \item \code{\link{add_ci.glm}} for generalized linear regression confidence intervals
#'   \item \code{\link{add_ci.lmerMod}} for linear mixed models confidence intervals
#' }
#'
#' Note that \code{add_ci} calculates confidence intervals for
#' \emph{fitted values}, not model coefficients.
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
#' @param tb A tibble or data frame of new data. \code{tb} can be
#'     the original data or new data.
#' @param fit An object of class \code{lm}, \code{glm}, or
#'     \code{lmerMod}. Predictions are made with this object.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names NULL or character vector of length two. If
#'     \code{NULL}, confidence bounds automatically will be named by
#'     \code{add_ci}, otherwise, the lower confidence bound will be
#'     named \code{ciNames[1]} and the upper confidence bound will be
#'     named \code{ciNames[2]}.
#' @param yhatName A string. Name of the vector of the predictions
#'     made for each observation in tb
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     confidence bounds attached.
#'
#' @examples
#' # Fit a linear model
#' fit <- lm(dist ~ speed, data = cars)
#' # Make a confidence interval for each observation in cars, and
#' # append to the data frame
#' add_ci(cars, fit)
#'
#' # Make new data
#' new_data <- cars[sample(NROW(cars), 10), ]
#' add_ci(new_data, fit)
#'
#' # Fit a Poisson model
#' fit2 <- glm(dist ~ speed, family = "poisson", data = cars)
#' # Append CIs
#' add_ci(cars, fit2)
#'
#' # Fit a linear mixed model using lme4
#' fit3 <- lme4::lmer(Reaction ~ Days + (1|Subject), data = lme4::sleepstudy)
#' # Append CIs
#' # Generally, you should use more than 100 bootstrap replicates
#' add_ci(lme4::sleepstudy, fit3, nSims = 100)
#'
#' # Fit a logistic model
#' fit4 <- glm(I(dist > 20) ~ speed, family = "binomial", data = cars)
#' # Append CIs
#' add_ci(cbind(cars, I(cars$dist > 20)), fit4)
#'
#'
#' @seealso \code{\link{add_pi}} for prediction intervals,
#'     \code{\link{add_probs}} for response level probabilities, and
#'     \code{\link{add_quantile}} for quantiles of the conditional
#'     response distribution.
#'
#' @export

add_ci <- function(tb, fit, alpha = 0.05, names = NULL, yhatName = "pred", ...){
  UseMethod("add_ci", fit)
}
