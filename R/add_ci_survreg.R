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

#' Confidence Intervals for the Expected Response of Parametric
#' Survival Models.
#'
#' This function is one of the methods of \code{add_ci}.
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{survreg}. Predictions are made
#'     with this object.
#' @param names \code{NULL} or a string. If \code{NULL}, quantiles
#'     automatically will be named by \code{add_quantile}, otherwise,
#'     they will be named \code{name}.
#' @param yhatName A string. Name of the vector of predictions.
#' @param alpha A number between 0 and 1. 1 - \code{alpha} is the
#'     confidence level of the intervals.
#' @param method A string. One of either \code{"parametric"} or
#'     \code{"boot"}.
#' @param nSims A positive integer. Set the number of simulated draws
#'     to use.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted expected values and level
#'     \emph{alpha} confidence levels attached.
#'
#' @seealso \code{\link{add_quantile.survreg}} for quantiles of
#'     \code{survreg} objects, \code{\link{add_pi.survreg}} for
#'     prediction intervals of \code{survreg} objects, and
#'     \code{\link{add_probs.survreg}} for response probabilities of
#'     \code{survreg} objects.
#'
#' @examples TODO
#'
#' @references TODO
#'
#' @export

add_ci.survreg <- function(tb, fit, p = 0.5,
                           names = NULL,
                           yhatName = "mean_pred",
                           alpha = 0.1,
                           method = "parametric", nSims = 2000,
                           ...){

    if (is.null(names)){
        names[1] <- paste("LCB", alpha/2, sep = "")
        names[2] <- paste("UCB", 1 - alpha/2, sep = "")
    }

    if ((names[1] %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe. Overwriting.")
    }


    if(method == "boot")
        boot_ci_survreg_expectation(tb, fit, yhatName,
                                    confint, alpha, names, nSims)
    else if(method == "parametric")
        parametric_ci_survreg_expectation(tb, fit, yhatName,
                                          alpha, names)
    else
        stop("method must be either 'boot' or 'parametric'")
}

## TODO checks for missing values and convergence
parametric_ci_survreg_expectation <- function(tb, fit, yhatName,
                                              alpha, names){
    distr <- fit$dist
    form <- formula(fit)
    m <- model.frame(form, tb)
    mat <- model.matrix(form, m)

    if(any(is.na(mat)))
        stop("Check tb for missingness")

    nPred <- dim(tb)[1]
    beta <- coef(fit)
    scale <- fit$scale

    if (distr == "loglogistic" && (scale >= 1))
        stop("Expected value is undefined for loglogistic distribution with scale >= 1")

    if (distr == "weibull")
        pred <- exp(mat %*% beta) * gamma(1 + scale)
    else if (distr == "exponential")
        pred <- exp(mat %*% beta)
    else if (distr == "loglogistic")
        pred <- exp(mat %*% beta) * gamma(1 + scale) * gamma(1 - scale)
    else if (distr == "lognormal")
        pred <- exp(c(mat %*% beta) + (scale^2) / 2)

    crit_val <- qnorm(p = 1 - alpha/2, mean = 0, sd = 1)
    cov_mat <- vcov(fit)

    if (fit$dist == "exponential")
        cov_mat <- cbind(rbind(cov_mat, 0), 0)

    d_g <- rep(NA, nPred)
    seYhat <- rep(NA, nPred)

    for(i in 1:nPred){
        if (distr == "weibull"){
            d_g_beta <- c(exp(mat[i,] %*% beta)) * mat[i,] * gamma(1 + scale)
            d_g_delta <- exp(mat[i,] %*% beta) *
                digamma(1 + scale) * gamma (1 + scale) * scale
        }
        if (distr == "exponential"){
            d_g_beta <- c(exp(mat[i,] %*% beta)) * mat[i,]
            d_g_delta <- 0
        }
        if (distr == "loglogistic"){
            d_g_beta <- c(exp(mat[i,] %*% beta)) * mat[i,] *
                gamma(1 + scale) * gamma(1 - scale)
            d_g_delta <- c(exp(mat[i,] %*% beta)) *
                (gamma(1 + scale) * digamma(1 - scale) * gamma(1 - scale) * (-scale) +
                 gamma(1 - scale) * digamma(1 + scale) * gamma(1 + scale) * scale)
        }
        if (distr == "lognormal"){
            d_g_beta <- exp(c(mat[i,] %*% beta) + (scale^2) / 2) * mat[i,]
            d_g_delta <- exp(c(mat[i,] %*% beta) + (scale^2) / 2) * (scale^2)
        }
        d_g_vec <- c(d_g_beta, d_g_delta)
        seYhat[i] <- sqrt(t(d_g_vec) %*% cov_mat %*% d_g_vec)
    }

    w <- exp(crit_val * seYhat / pred)
    lwr <- pred / w
    upr <- pred * w

    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- c(pred)

    tb[[names[1]]] <- as.numeric(lwr)
    tb[[names[2]]] <- as.numeric(upr)

    tibble::as_data_frame(tb)
}
