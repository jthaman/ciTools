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

#' Response Probabilities for Generalized Linear Mixed Model Predictions
#'
#' This function is one of the methods for \code{add_probs}, and is
#' called automatically when \code{add_probs} is used on a \code{fit}
#' of class \code{glmerMod}. Probabilities are approximate and
#' determined via a simulation. This function is experimental.
#'
#' @param df A data frame of new data.
#' @param fit An object of class \code{glmerMod}.
#' @param q A double. A quantile of the response distribution.
#' @param name \code{NULL} or character vector of length one. If
#'     \code{NULL}, response probabilities automatically will be named
#'     by \code{add_probs},
#' @param yhatName \code{NULL} or a string. Name of the predictions
#'     vector.
#' @param type A string. Must be \code{"boot"}, If \code{type =
#'     "boot"}, then \code{add_ci} calls \code{lme4::simulate} to
#'     calculate the probabilities.
#' @param comparison A string. If \code{comparison = "<"}, then
#'     \eqn{Pr(Y|x < q)} is calculated for each observation in
#'     \code{df}. Default is "<". Must be "<" or ">" for objects of
#'     class \code{lm} or \code{lmerMod}. If \code{fit} is a
#'     \code{glm} or \code{glmerMod}, then \code{comparison} also may
#'     be \code{"<="} , \code{">="} , or \code{"="}.
#' @param includeRanef A logical. Default is \code{TRUE}. Set whether
#'     the predictions and intervals should be made conditional on the
#'     random effects. If \code{FALSE}, random effects will not be
#'     included.
#' @param nSims A positive integer.  Controls the number of bootstrap
#'     replicates if \code{type = "boot"}.
#' @param ... Additional arguments.
#' @return A dataframe, \code{df}, with predicted values and estimated
#'     probabilities attached.
#'
#' @seealso \code{\link{add_pi.glmerMod}} for prediction intervals
#'     of \code{glmerMod} objects, \code{\link{add_ci.glmerMod}} for
#'     confidence intervals of \code{glmerMod} objects, and
#'     \code{\link{add_quantile.glmerMod}} for response quantiles of
#'     \code{glmerMod} objects.
#'
#' @details If \code{IncludeRanef} is False, random slopes and intercepts are set to 0. Unlike in
#'   `lmer` fits, settings random effects to 0 does not mean they are marginalized out. Consider
#'   generalized estimating equations if this is desired.
#'
#' @examples
#' n <- 300
#' x <- runif(n)
#' f <- factor(sample(1:5, size = n, replace = TRUE))
#' y <- rpois(n, lambda = exp(1 - 0.05 * x * as.numeric(f) + 2 * as.numeric(f)))
#' df <- data.frame(x = x, f = f, y = y)
#' fit <- lme4::glmer(y ~ (1+x|f), data=df, family = "poisson")
#'
#' add_probs(df, fit, name = "p0.6", q = 0.6, nSims = 500)
#'
#' @export

add_probs.glmerMod <- function(df, fit,
                               q, name = NULL, yhatName = "pred", comparison = "<",
                               type = "boot", includeRanef = TRUE,
                               nSims = 10000, ...){

    if (!is.null(fit@optinfo$conv$lme4$code))
        warning ("Coverage probabilities may be inaccurate if the model failed to converge")

    if(fit@resp$family$family == "binomial")
        stop("Prediction Intervals are not useful if the response is Bernoulli")

    if (is.null(name) & (comparison == "<"))
        name <- paste("prob_less_than", q, sep="")
    if (is.null(name) & (comparison == ">"))
        name <- paste("prob_greater_than", q, sep="")
    if (is.null(name) & (comparison == "<="))
        name <- paste("prob_less_than_or_equal_to", q, sep="")
    if (is.null(name) & (comparison == ">="))
        name <- paste("prob_greater_than_or_equal_to", q, sep="")
    if (is.null(name) & (comparison == "="))
        name <- paste("prob_equal_to", q, sep="")

    if ((name %in% colnames(df))) {
        warning ("These probabilitiess may have already been appended to your dataframe. Overwriting.")
    }

    if (type == "boot")
        bootstrap_probs_glmermod(df, fit, q, name, includeRanef, nSims, yhatName, comparison)
    else
        stop("Incorrect type specified!")
}

bootstrap_probs_glmermod <- function(df, fit, q, name, includeRanef, nSims, yhatName, comparison) {

    if (includeRanef) {
        rform = NULL
    } else {
        rform = NA
    }

  gg <- simulate(fit, newdata = df, re.form = rform, nsim = nSims)
  gg <- as.matrix(gg)
  probs <- apply(gg, 1, FUN = calc_prob, quant = q, comparison = comparison)
  out <- predict(fit, df, re.form = rform, type = "response")

  if(is.null(df[[yhatName]]))
    df[[yhatName]] <- out
  df[[name]] <- probs
  data.frame(df)
}
