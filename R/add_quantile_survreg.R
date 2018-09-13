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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ciTools. If not, see <http://www.gnu.org/licenses/>.

#' Confidence Intervals for a Quantile of a Parametric Survival Model
#'
#' This function is one of the methods of \code{add_quantile}.
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{survreg}. Predictions are made
#'     with this object.
#' @param p A real number between 0 and 1. Sets the probability level
#'     of the quantiles.
#' @param name \code{NULL} or a string. If \code{NULL}, quantiles
#'     automatically will be named by \code{add_quantile}, otherwise,
#'     they will be named \code{name}.
#' @param yhatName A string. Name of the vector of predictions.
#' @param confint A logical. If \code{TRUE}, confidence intervals for
#'     the quantiles are also appended to \code{tb}.
#' @param method A string. One of either \code{"parametric"} or
#'     \code{"boot"}.
#' @param nSims A positive integer. Set the number of simulated draws
#'     to use.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values and level
#'     \emph{p} quantiles attached.
#'
#' @seealso \code{\link{add_ci.survreg}} for confidence intervals for
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

add_quantile.survreg <- function(tb, fit, p = 0.5,
                                 yhatName = "median_pred", confint = TRUE,
                                 alpha = 0.1, name = NULL,
                                 method = "parametric", nSims = 2000,
                                 ...){
    if (p <= 0 || p >= 1)
        stop ("p should be in (0,1)")
    if (is.null(name)){
        name[1] <- paste("quantile", p, sep="")
        name[2] <- paste("lcb")
        name[3] <- paste("ucb")
    }
    if ((name[1] %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe. Overwriting.")
    }

    if(method == "boot")
        boot_ci_survreg_quantile(tb, fit, p, yhatName,
                                 confint, alpha, name, nSims)
    else if(method == "parametric")
        parametric_ci_survreg_quantile(tb, fit, p, yhatName,
                                       confint, alpha, name)
    else
        stop("method must be either 'boot' or 'parametric'")
}

boot_ci_survreg_quantile <- function(tb, fit, p, yhatName, confint,
                                     alpha, name, nSims){
    nPred <- dim(tb)[1]
    out <- survival:::predict.survreg(fit, tb, se.fit = TRUE,
                                      type = "quantile", p = p)
    pred <- out$fit

    if (confint){
        boot_mat <- matrix(NA, nrow = nSims, ncol = nPred)
        for (i in 1:nSims){
            temp <- tb[sample(1:nPred, size = nPred, replace = TRUE),]
            boot_fit <- survival::survreg(formula(fit$terms), data = temp,
                                          dist = fit$dist)
            boot_pred <- survival:::predict.survreg(boot_fit, tb,
                                                    type = "quantile", p = p)
            boot_mat[i,] <- boot_pred
        }
        lwr = apply(boot_mat, 2, quantile, probs = alpha / 2)
        upr = apply(boot_mat, 2, quantile, probs = 1 - alpha / 2)
    }

    if (is.null(tb[[name[1]]]))
        tb[[name[1]]] <- pred

    if (confint){
        tb[[name[2]]] <- lwr
        tb[[name[3]]] <- upr
    }
    tibble::as_data_frame(tb)
}

## this assumes log link function
parametric_ci_survreg_quantile <- function(tb, fit, p, yhatName, confint,
                                           alpha, name){
    out <- predict(fit, tb, se.fit = TRUE, type = "quantile", p = p)

    if (confint){
        crit_val <- qnorm(p = 1 - alpha/2, mean = 0, sd = 1)
        pred <- out$fit
        se <- out$se.fit
        w <- exp(crit_val * se / pred)
        upr <- pred * w
        lwr <- pred / w
    }

    if (is.null(tb[[name[1]]]))
        tb[[name[1]]] <- pred

    if (confint){
        tb[[name[2]]] <- lwr
        tb[[name[3]]] <- upr
    }
    tibble::as_data_frame(tb)
}
