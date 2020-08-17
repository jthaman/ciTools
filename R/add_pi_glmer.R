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

#' Prediction Intervals for Generalized Linear Mixed Model Predictions
#'
#' This function is one of the methods for \code{add_pi}, and is
#' called automatically when \code{add_pi} is used on a \code{fit} of
#' class \code{glmerMod}. This function is experimental.
#'
#' Prediction intervals are approximate and determined by simulation
#' through the \code{simulate} function distributed with \code{lme4}.
#'
#' @param df A data frame of new data.
#' @param fit An object of class \code{glmerMod}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, prediction bounds automatically will be named by
#'     \code{add_pi}, otherwise, the lower prediction bound will be
#'     named \code{names[1]} and the upper prediction bound will be
#'     named \code{names[2]}.
#' @param yhatName \code{NULL} or a string. The name of the
#'     predictions vector.
#' @param type A string. Must be \code{"boot"}, If \code{type =
#'     "boot"}, then \code{add_ci} calls \code{lme4::bootMer} to
#'     calculate the confidence intervals.
#' @param includeRanef A logical. Default is \code{TRUE}. Set whether
#'     the predictions and intervals should be made conditional on the
#'     random effects. If \code{FALSE}, random effects will not be
#'     included.
#' @param nSims A positive integer.  Controls the number of bootstrap
#'     replicates.
#' @param ... Additional arguments.
#' @return A dataframe, \code{df}, with predicted values, upper and lower
#'     prediction bounds attached.
#'
#' @seealso \code{\link{add_ci.glmerMod}} for confidence intervals
#'     of \code{glmerMod} objects, \code{\link{add_probs.glmerMod}} for
#'     conditional probabilities of \code{glmerMod} objects, and
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
#' add_pi(df, fit, names = c("LPB", "UPB"), nSims = 500)
#'
#' @export

add_pi.glmerMod <- function(df, fit,
                            alpha = 0.05, names = NULL, yhatName = "pred",
                            type = "boot", includeRanef = TRUE,
                            nSims = 10000, ...){

    if (!is.null(fit@optinfo$conv$lme4$code))
        warning ("Coverage probabilities may be inaccurate if the model failed to converge")

    if(fit@resp$family$family == "binomial")
        warning("Prediction Intervals are not useful if the response is Bernoulli")

    if (is.null(names)){
        names[1] <- paste("LPB", alpha/2, sep = "")
        names[2] <- paste("UPB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(df))) {
        warning ("These PIs may have already been appended to your dataframe. Overwriting.")
    }

    if(fit@resp$family$family %in% c("poisson", "quasipoisson", "binomial"))
        warning("The response is not continuous, so Prediction Intervals are approximate")

    if (type == "boot")
        bootstrap_pi_glmermod(df, fit, alpha, names, includeRanef, nSims, yhatName)
    else
        stop("Incorrect type specified!")
}

bootstrap_pi_glmermod <- function(df, fit, alpha, names, includeRanef, nSims, yhatName) {

    if (includeRanef) {
        rform = NULL
    } else {
        rform = NA
    }

    gg <- simulate(fit, newdata = df, re.form = rform, nsim = nSims)

    gg <- as.matrix(gg)
    lwr <- apply(gg, 1, FUN = quantile, probs = alpha/2)
    upr <- apply(gg, 1, FUN = quantile, probs = 1 - alpha / 2)

  out <- predict(fit, df, re.form = rform, type = "response")

  if(is.null(df[[yhatName]]))
    df[[yhatName]] <- out
  df[[names[1]]] <- lwr
  df[[names[2]]] <- upr

  data.frame(df)
}
