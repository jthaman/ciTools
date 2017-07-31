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

#' Prediction Intervals for Linear Model Predictions
#'
#' This function is one of the methods for \code{add_pi} and is
#' automatically called when an object of class \code{lm} is passed to
#' to \code{add_pi}.
#'
#' Prediction intervals for \code{lm} objects are calculated
#' parametrically. This function is essentially just a wrapper for
#' \code{predict(fit, tb, interval = "prediction")} if \code{fit} is a
#' linear model. If \code{log_response = TRUE}, prediction intervals
#' for the response are calculated parametrically, then the
#' exponential function is applied to transform them to the original
#' scale.
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class lm. Predictions are made with this
#'     object.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, prediction bounds automatically will be named by
#'     \code{add_pi}, otherwise, the lower prediction bound will be
#'     named \code{names[1]} and the upper prediction bound will be
#'     named \code{names[2]}.
#' @param yhatName A string. Name of the predictions vector.
#' @param log_response A logical. If TRUE, prediction intervals will be
#'     generated at the \emph{response level} of a log-linear model:
#'     \eqn{\log(Y) = X\beta + \epsilon}. Again, these intervals will
#'     be on the scale of the original response, Y.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#' 
#' @seealso \code{\link{add_ci.lm}} for confidence intervals for
#'     \code{lm} objects. \code{\link{add_probs.lm}} for conditional
#'     probabilities of \code{lm} objects, and
#'     \code{\link{add_quantile.lm}} for response quantiles of
#'     \code{lm} objects.
#'
#' @examples
#' # Fit a linear model
#' fit <- lm(dist ~ speed, data = cars)
#' # Add prediction intervals and fitted values to the original data
#' add_pi(cars, fit)
#' 
#' # Try to add predictions to a data frame of new data
#' new_data <- cars[sample(NROW(cars), 10), ]
#' add_pi(new_data, fit)
#' 
#' # Try a different confidence level
#' add_pi(cars, fit, alpha = 0.5)
#' 
#' # Add custom names to the prediction bounds.
#' add_pi(cars, fit, alpha = 0.5, names = c("lwr", "upr"))
#' @export

add_pi.lm <- function(tb, fit, alpha = 0.05, names = NULL,
                      yhatName = "pred", log_response = FALSE, ...){
    
    if (log_response)
        add_pi_lm_log(tb, fit, alpha, names, yhatName)

    else {
        if (is.null(names)){
            names[1] <- paste("LPB", alpha/2, sep = "")
            names[2] <- paste("UPB", 1 - alpha/2, sep = "")
        }
        if ((names[1] %in% colnames(tb))) {
            warning ("These PIs may have already been appended to your dataframe. Overwriting.")
        }
        out <- predict(fit, tb, interval = "prediction", level = 1 - alpha)
        if(is.null(tb[[yhatName]]))
            tb[[yhatName]] <- out[, 1]
        if (is.null(tb[[names[1]]]))
            tb[[names[1]]] <- out[, 2]
        if (is.null(tb[[names[2]]]))
            tb[[names[2]]] <- out[, 3]
        tibble::as_data_frame(tb)
    }
}

