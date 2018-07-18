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

#' Response Quantiles for Generalized Linear Mixed Model Predictions
#'
#' This function is one of the methods for \code{add_pi}, and is
#' called automatically when \code{add_pi} is used on a \code{fit} of
#' class \code{glmerMod}.
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{glmerMod}.
#' @param p A real number between 0 and 1. Sets the probability level
#'     of the quantiles.
#' @param name \code{NULL} or a string. If \code{NULL}, quantile
#'     automatically will be named by \code{add_quantile}
#' @param yhatName \code{NULL} or a string. Name of the predictions vector.
#' @param type A string. Must be \code{"boot"}, If \code{type =
#'     "boot"}, then \code{add_ci} calls \code{lme4::simulate} to
#'     calculate the confidence intervals. This method may be time
#'     consuming, but is applicable with random slope and random
#'     intercept models.
#' @param includeRanef A logical. Default is \code{TRUE}. Set whether
#'     the predictions and intervals should be made conditional on the
#'     random effects. If \code{FALSE}, random effects will not be
#'     included.
#' @param nSims A positive integer.  Controls the number of bootstrap
#'     replicates.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values and quantiles attached.
#'
#' @seealso \code{\link{add_pi.glmerMod}} for prediction intervals
#'     of \code{glmerMod} objects, \code{\link{add_probs.glmerMod}} for
#'     conditional probabilities of \code{glmerMod} objects, and
#'     \code{\link{add_ci.glmerMod}} for confidence intervals of
#'     \code{glmerMod} objects.
#'
#' @examples
#' n <- 300
#' x <- runif(n)
#' f <- factor(sample(1:5, size = n, replace = TRUE))
#' y <- rpois(n, lambda = exp(1 - 0.05 * x * as.numeric(f) + 2 * as.numeric(f)))
#' tb <- tibble::tibble(x = x, f = f, y = y)
#' fit <- lme4::glmer(y ~ (1+x|f), data=tb, family = "poisson")
#'
#' add_quantile(tb, fit, name = "quant0.6", p = 0.6, nSims = 500)
#'
#' @export

add_quantile.glmerMod <- function(tb, fit,
                                  p, name = NULL, yhatName = "pred",
                                  type = "boot", includeRanef = TRUE,
                                  nSims = 10000, ...){

    if (!is.null(fit@optinfo$conv$lme4$code))
        warning ("Coverage probabilities may be inaccurate if the model failed to converge")
    if(fit@resp$family$family == "binomial")
        warning("Quantiles are not useful if the response is Bernoulli")
    if(fit@resp$family$family %in% c("poisson", "quasipoisson", "binomial"))
        warning("The response is not continuous, so regression quantiles are approximate")
    if (p <= 0 || p >= 1)
        stop ("p should be in (0,1)")
    if (is.null(name))
        name <- paste("quantile", p, sep="")
    if ((name %in% colnames(tb)))
        warning ("These quantiles may have already been appended to your dataframe. Overwriting.")
    if (type == "boot")
        bootstrap_quant_glmermod(tb, fit, p, name, includeRanef, nSims, yhatName)
    else
        stop("Incorrect type specified!")
}

bootstrap_quant_glmermod <- function(tb, fit, p, name, includeRanef, nSims, yhatName) {

    if (includeRanef) {
        rform = NULL
    } else {
        rform = NA
    }

    gg <- simulate(fit, newdata = tb, re.form = rform, nsim = nSims)
    gg <- as.matrix(gg)
    quant_out <- apply(gg, 1, FUN = quantile, probs = p)
    out <- predict(fit, tb, re.form = rform, type = "response")

    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[name]] <- quant_out

    tibble::as_data_frame(tb)
}
