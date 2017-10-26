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
#' called automatically when \code{add_probs} is used on a \code{fit} of
#' class \code{glmerMod}. 
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{lmerMod}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param name \code{NULL} or character vector of length one. If
#'     \code{NULL}, response probabilities automatically will be named
#'     by \code{add_probs},
#' @param yhatName A string. Name of the predictions vector.
#' @param type A string. Must be \code{"boot"}, If \code{type =
#'     "boot"}, then \code{add_ci} calls \code{lme4::simulate} to
#'     calculate the probabilities.
#' @param includeRanef A logical. Default is \code{TRUE}. Set whether
#'     the predictions and intervals should be made conditional on the
#'     random effects. If \code{FALSE}, random effects will not be
#'     included.
#' @param nSims A positive integer.  Controls the number of bootstrap
#'     replicates if \code{type = "boot"}.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     confidence bounds attached.
#'
#' @seealso \code{\link{add_pi.glmerMod}} for prediction intervals
#'     of \code{glmerMod} objects, \code{\link{add_ci.glmerMod}} for
#'     confidence intervals of \code{glmerMod} objects, and
#'     \code{\link{add_quantile.glmerMod}} for response quantiles of
#'     \code{glmerMod} objects.
#'
#' @references
#' TODO
#'
#' @examples
#' TODO
#'
#' @export

add_probs.glmerMod <- function(tb, fit, 
                               q, names = NULL, yhatName = "pred",
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

    if ((name %in% colnames(tb))) {
        warning ("These PIs may have already been appended to your dataframe. Overwriting.")
    }

    if (type == "boot")
        bootstrap_probs_glmermod(tb, fit, alpha, names, includeRanef, nSims, yhatName)
    else
        stop("Incorrect type specified!")
}

bootstrap_probs_glmermod <- function(tb, fit, alpha, names, includeRanef, nSims, yhatName) {

    if (includeRanef) { 
        rform = NULL
        my_pred <- my_pred_full
    } else {
        rform = NA
        my_pred <- my_pred_fixed
    }
        
    gg <- simulate(fit, newdata = tb, re.form = rform, nsim = nSims)
    gg <- as.matrix(gg)
    probs <- apply(gg, 1, FUN = calc_prob, quant = q, comparison = comparison)
    out <- predict(fit, tb, re.form = rform, type = "response")

    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[name]] <- probs
    tibble::as_data_frame(tb)
}

