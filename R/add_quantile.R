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

#' Add Regression Quantiles to Data Frames
#'
#' This is a generic function to append regression quantiles to a data
#' frame. A regression quantile \emph{q} is a point such that
#' \eqn{Pr(Response | Covariates < q) = p}. These quantiles are
#' generated for the fitted value of each observation in
#' \code{tb}. Quantiles are then appended to \code{tb} and returned to
#' the user as a tibble.
#'
#' For more specific information about the arguments that are applicable
#' for each type of model, consult:
#'
#' \itemize{
#'   \item \code{\link{add_quantile.lm}} for linear regression response quantiles
#'   \item \code{\link{add_quantile.glm}} for generalized linear regression response quantiles
#'   \item \code{\link{add_quantile.lmerMod}} for linear mixed models response quantiles
#' }
#'
#' Note: the quantiles that \code{add_quantile} calculates are on the
#' distribution of \eqn{Y|x}, not \eqn{E[Y|x]}. That is, they use the
#' same distribution that determines a prediction interval, not the
#' distribution that determines a confidence interval.

#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{lm}, \code{glm}, or
#'     \code{lmerMod}. Predictions are made with this object.
#' @param p A double. A probability that determines the quantile. Must
#'     be between 0 and 1.
#' @param name \code{NULL} or a string. If \code{NULL},
#'     quantiles automatically will be named by \code{add_quantile()},
#'     otherwise, the quantiles will be named \code{name} in the
#'     returned tibble.
#' @param yhatName A string. Name of the vector of predictions.
#' @param ... Additional arguments
#' @return A tibble, \code{tb}, with predicted values and
#'     level-\emph{p} quantiles attached.
#' 
#' @examples
#'
#' # Fit a linear model
#' fit <- lm(dist ~ speed, data = cars)
#'
#' # Find the 0.4 quantile (or 40th percentile) of new distances for
#' # each observations in cars, conditioned on the linear model.
#' add_quantile(cars, fit, p = 0.4)
#'
#' # Fit a Poisson model
#' fit2 <- glm(dist ~ speed, family = "poisson", data = cars)
#' # Find the 0.4 quantile (or 40th percentile) of new distances for
#' # each observation in cars, conditioned on the Poisson model.
#' add_quantile(cars, fit2, p = 0.4)
#'
#' # Fit a random intercept linear mixed model
#' fit3 <- lme4::lmer(Reaction ~ Days + (1|Subject), data = lme4::sleepstudy)
#' # Find the 0.4 quantile (or 40 percentile) of reaction times for
#' # each observation in the sleepstudy data. Condition on the model and random effects.
#' add_quantile(lme4::sleepstudy, fit3, p = 0.4, type = "parametric")
#'
#'
#' @seealso \code{\link{add_ci}} for confidence intervals,
#'     \code{\link{add_probs}} for response level probabilities, and
#'     \code{\link{add_pi}} for prediction intervals
#' 
#' @export

add_quantile <- function(tb, fit, p, name = NULL, yhatName = "pred", ...){
  UseMethod("add_quantile", fit)
}

