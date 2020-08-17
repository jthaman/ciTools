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
#' \code{df}. These confidence intervals are then appended to
#' \code{df} and returned to the user as a data frame. The \code{fit} may
#' be a linear, log-linear, linear mixed, generalized linear model,
#' generalized linear mixed, or accelerated failure time model.
#'
#' For more specific information about the arguments that are applicable
#' in each method, consult:
#'
#' \itemize{
#'   \item \code{\link{add_ci.lm}} for linear model confidence intervals
#'   \item \code{\link{add_ci.glm}} for generalized linear model confidence intervals
#'   \item \code{\link{add_ci.lmerMod}} for linear mixed model confidence intervals
#'   \item \code{\link{add_ci.glmerMod}} for generalized linear mixed model confidence intervals
#'   \item \code{\link{add_ci.survreg}} for accelerated failure time confidence intervals
#' }
#'
#' Note that \code{add_ci} calculates confidence intervals for
#' \emph{fitted values}, not model coefficients. For confidence
#' intervals of model coefficients, see \code{confint}.
#'
#' @import stats
#' @import lme4
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @importFrom arm se.ranef
#' @importFrom arm sim
#' @importFrom MASS rnegbin
#' @importFrom MASS mvrnorm
#' @importFrom MASS glm.nb
#' @importFrom methods new
#' @importFrom boot boot
#' @importFrom boot boot.ci
#' @importFrom utils packageVersion
#'
#' @param df A data frame of new data. \code{df} can be the
#'     original data or new data.
#' @param fit An object of class \code{lm}, \code{glm},
#'     \code{lmerMod}, \code{glmerMod}, or \code{survreg}. Predictions
#'     are made with this object.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, confidence bounds automatically will be named by
#'     \code{add_ci}, otherwise, the lower confidence bound will be
#'     named \code{names[1]} and the upper confidence bound will be
#'     named \code{names[2]}.
#' @param yhatName A character vector of length one. Name of the
#'     vector of the predictions made for each observation in \code{df}
#' @param ... Additional arguments.
#' @return A dataframe, \code{df}, with predicted values, upper and lower
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

add_ci <- function(df, fit, alpha = 0.05, names = NULL, yhatName = "pred", ...){
    UseMethod("add_ci", fit)
}
