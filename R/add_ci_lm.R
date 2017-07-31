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

#' Confidence Intervals for Linear Model Predictions.
#'
#' This function is one of the methods in \code{add_ci} and is
#' automatically called when an object is class \code{lm} is passed to
#' to \code{add_ci}.
#'
#' Confidence intervals for \code{lm} objects are calculated
#' parametrically, this function is essentially a wrapper for
#' \code{predict(fit, tb, interval = "confidence")} if \code{fit} is a
#' linear model. If \code{log_response = TRUE}, confidence intervals
#' for the response are calculated using Wald's Method. See Meeker and
#' Escobar (1998) for details.
#'
#' @param tb A tibble or data frame.
#' @param fit An object of class \code{lm}. Predictions are made with this
#'     object.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, confidence bounds will automatically be named by
#'     \code{add_ci}, otherwise, the lower confidence bound will be
#'     named \code{names[1]} and the upper confidence bound will be
#'     named \code{names[2]}.
#' @param log_response logical. Default is \code{FALSE}. If
#'     \code{TRUE}, confidence intervals will be generated for the
#'     \emph{response level} of a log-linear model:  \eqn{\log(Y) =
#'     X\beta + \epsilon}.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     confidence bounds attached.
#'
#' @seealso \code{\link{add_pi.lm}} for prediction intervals for
#'     \code{lm} objects. \code{\link{add_probs.lm}} for conditional
#'     probabilities of \code{lm} objects, and
#'     \code{\link{add_quantile.lm}} for response quantiles of
#'     \code{lm} objects.
#'
#' @examples
#' # Fit a Linear model
#' fit <- lm(dist ~ speed, data = cars)
#' # Get fitted values for each observation in cars, and append
#' # confidence intervals
#' add_ci(cars, fit)
#' # Try a different confidence level
#' add_ci(cars, fit, alpha = 0.5)
#' # Try custom names for the confidence bounds 
#' add_ci(cars, fit, alpha = 0.5, names = c("lwr", "upr"))
#' 
#' @export

add_ci.lm <- function(tb, fit, alpha = 0.05, names = NULL, log_response = FALSE, ...){
    if (log_response)
        add_ci_lm_log(tb, fit, alpha, names)
    else {

        if (is.null(names)){
            names[1] <- paste("LCB", alpha/2, sep = "")
            names[2] <- paste("UCB", 1 - alpha/2, sep = "")
        }
        if ((names[1] %in% colnames(tb))) {
            warning ("These CIs may have already been appended to your dataframe. Overwriting.")
        }
        out <- predict(fit, tb, interval = "confidence", level = 1 - alpha)
        if(is.null(tb[["pred"]]))
            tb[["pred"]] <- out[, 1]
        if (is.null(tb[[names[1]]]))
            tb[[names[1]]] <- out[, 2]
        if (is.null(tb[[names[2]]]))
            tb[[names[2]]] <- out[, 3]
        tibble::as_data_frame(tb)
    } 
}

