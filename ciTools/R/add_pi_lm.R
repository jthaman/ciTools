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

#' Prediction Intervals for Linear Model Predictions.
#'
#' This function is one of the methods for \code{add_pi}. This
#' function is essentially a wrapper for the \code{predict} function
#' in \code{R}.
#'
#' @param tb A tibble or Data Frame.
#' @param fit An object of class lm. Predictions are made with this
#'     object.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, prediction bounds will automatically be named by
#'     \code{add_pi}, otherwise, the lower prediction bound will be
#'     named \code{names[1]} and the upper prediction bound will be
#'     named \code{names[2]}.
#' @param log_response logical. If TRUE, prediction intervals will be
#'     generated for the prediction made with a log-linear model:
#'     \eqn{\log(Y) = X\beta + \epsilon}. These intervals will be on
#'     the scale of the original response, Y.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#' 
#' @export

add_pi.lm <- function(tb, fit, alpha = 0.05, names = NULL, log_response = FALSE, ...){
    
    if (log_response)
        add_pi_lm_log(tb, fit, alpha, names)

    else {
        if (is.null(names)){
            names[1] <- paste("LPB", alpha/2, sep = "")
            names[2] <- paste("UPB", 1 - alpha/2, sep = "")
        }
        if ((names[1] %in% colnames(tb))) {
            warning ("These PIs may have already been appended to your dataframe. Overwriting.")
        }
        out <- predict(fit, tb, interval = "prediction", level = 1 - alpha)
        if(is.null(tb[["pred"]]))
            tb[["pred"]] <- out[, 1]
        if (is.null(tb[[names[1]]]))
            tb[[names[1]]] <- out[, 2]
        if (is.null(tb[[names[2]]]))
            tb[[names[2]]] <- out[, 3]
        tibble::as_data_frame(tb)
    }
}

