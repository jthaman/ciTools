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

#' Add Regression Probabilities to Data Frames
#'
#' This is a generic function to append response level probabilities
#' to a data frame. A response level probability (conditioned on the
#' model and covariates), such as \eqn{Pr(Response|Covariates < 10)},
#' is generated for the fitted value of each observation in
#' \code{tb}. These probabilities are then appended to \code{tb} and
#' returned to the user as a tibble.
#' 
#' For more specific information about the arguments that are useful
#' in each method, consult:
#'
#' \itemize{
#'   \item \code{\link{add_probs.lm}} for linear regression response probabilities
#'   \item \code{\link{add_probs.glm}} for generalized linear regression response probabilities
#'   \item \code{\link{add_probs.lmerMod}} for linear mixed models response probabilities
#' }
#'
#' Note that the probabilities calculated by \code{add_probs} are on
#' the distribution of \eqn{Y|x}, not \eqn{E[Y|x]}. That is, they use
#' the same distribution from which a prediction interval is
#' determined, not the distribution that determines a confidence
#' interval.
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{lm}, \code{glm}, or
#'     \code{lmerMod}. Predictions are made with this object.
#' @param q A real number. A quantile of the conditional response
#'     distribution.
#' @param name \code{NULL} or character vector of length one. If
#'     \code{NULL}, probabilities automatically will be named by
#'     \code{add_probs}, otherwise, the probabilities will be named
#'     \code{name} in the returned tibble.
#' @param yhatName A character vector of length one. Names of the
#' @param comparison A string. If \code{comparison = "<"}, then
#'     \eqn{Pr(Y|x < q)} is calculated for each observation in
#'     \code{tb}. Default is "<". Must be "<" or ">" for objects of
#'     class \code{lm} or \code{lmerMod}. If \code{fit} is a
#'     \code{glm}, then \code{comparison} also may be \code{"<="} ,
#'     \code{">="} , or \code{"="}.
#' @param ... Additional arguments
#' @return A tibble, \code{tb}, with predicted values and
#'     probabilities attached.
#' 
#' @examples
#' # Define a model 
#' fit <- lm(dist ~ speed, data = cars)
#'
#' # Calculate the probability that the probability that a new
#' # dist is less than 20 (Given the model).
#' add_probs(cars, fit, q = 20)
#'
#' # Calculate the probability that a new
#' # dist is greater than 20 (Given the model).
#' add_probs(cars, fit, q = 20, comparison = ">")
#' 
#' # Try a different model fit.
#' fit2 <- glm(dist ~ speed, family = "poisson", data = cars) 
#' add_probs(cars, fit2, q = 20)
#'
#' # Try a different model fit.
#' fit3 <- lme4::lmer(Reaction ~ Days + (1|Subject), data = lme4::sleepstudy)
#' add_probs(lme4::sleepstudy, fit3, q = 300, type = "parametric")
#'
#' # As above, but do not condition on the random effects.
#' add_probs(lme4::sleepstudy, fit3, q = 300, type = "parametric", includeRanef = FALSE)
#'
#' @seealso \code{\link{add_ci}} for confidence intervals,
#'     \code{\link{add_quantile}} for response level quantiles, and
#'     \code{\link{add_pi}} for prediction intervals.
#'
#' 
#' @export

add_probs <- function(tb, fit, q, name = NULL, yhatName = "pred", comparison, ...){
    UseMethod("add_probs", fit)
}


