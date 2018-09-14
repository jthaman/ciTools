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

#' Confidence Intervals for the Mean Survival Time of Parametric
#' Survival Models.
#'
#' This function is one of the methods for \code{add_ci}, and is
#' called automatically when \code{add_ci} is used on a \code{fit} of
#' class \code{survreg}.
#'
#' \code{add_ci.survreg} calculates confidence intervals for the mean
#' survival time of several accelerated failure time (AFT) models
#' including exponential, lognormal, weibull, and loglogistic
#' models. AFT models must be fit with the \code{survreg} function in
#' the \code{survival} package. Confidence intervals are formed either
#' parametrically via the Delta method, or through a non-parametric
#' bootstrap resampling procedure. Generally, both of the these
#' methods perform well under mild to moderate right censoring.
#'
#' TODO: finish coding bootstrap method.
#'
#' Traditionally, survival time predictions are made with the median
#' survival time. For forming confidence intervals for the median
#' survival time (or any quantile of the survival time distribution),
#' see \code{\link{add_quantile.survreg}}.
#'
#' Note: The expected survival time of a loglogistic model with scale
#' >= 1 does not exist. Otherwise, expected survival times exist for
#' each of the four AFT models considered in \code{add.ci_survreg}.
#'
#' Note: Due to a limitation, the \code{Surv} object must be specified in
#' \code{survreg} function call. See the examples section for one way
#' to do this.
#'
#' Note: \code{add_ci.survreg} cannot inspect the convergence of
#' \code{fit}. Poor maximum likelihood estimates will result in poor
#' confidence intervals. Inspect any warning messages given from
#' \code{survreg}.
#'
#' @param tb A tibble or data frame of new data on which to form
#'     predictions and confidence interval.
#' @param fit An object of class \code{survreg}. Predictions are made
#'     with this object.
#' @param names \code{NULL} or a string. If \code{NULL}, quantiles
#'     automatically will be named by \code{add_quantile}, otherwise,
#'     they will be named \code{name}.
#' @param yhatName A string. Name of the vector of predictions. The
#'     default name is \code{mean_pred}.
#' @param alpha A number between 0 and 1. 1 - \code{alpha} is the
#'     confidence level of the intervals.
#' @param method A string. One of either \code{"parametric"} or
#'     \code{"boot"}. If \code{method = "parametric"}, Delta method
#'     intervals are calculated. If \code{method = "boot"}
#'     nonparametric bootstrap intervals are calculated.
#' @param nSims A positive integer. Set the number of simulated draws
#'     to use. A value greater than or equal to 2000 is recommended.
#' @param ... Additional arguments.
#'
#' @return A tibble, \code{tb}, with predicted expected values and
#'     level \emph{alpha} confidence levels attached.
#'
#' @seealso \code{\link{add_quantile.survreg}} for quantiles of the
#'     survival time distribution of \code{survreg} objects,
#'     \code{\link{add_pi.survreg}} for prediction intervals of
#'     \code{survreg} objects, and \code{\link{add_probs.survreg}} for
#'     survival probabilities of \code{survreg} objects.
#'
#' @examples
#' ## Define a data set.
#' tb <- survival::stanford2
#' ## remove a covariate with missing values.
#' tb <- tb[, 1:4]
#' ## next, create the Surv object inside the survreg call:
#' fit <- survival::survreg(survival::Surv(time, status) ~ age + I(age^2),
#'                          data = tb, dist = "lognormal")
#' add_ci(tb, fit, alpha = 0.1, names = c("lwr", "upr"))
#'
#' ## Try a different model:
#' fit2 <- survival::survreg(survival::Surv(time, status) ~ age + I(age^2),
#'                           data = tb, dist = "weibull")
#' add_ci(tb, fit2, alpha = 0.1, names = c("lwr", "upr"))
#'
#' @references
#' For descriptions of the log-location scale models supported:
#' Meeker, William Q., and Luis A. Escobar. Statistical methods for reliability data. John Wiley & Sons, 2014. (Chapter 4)
#'
#' For a description of the multivariate Delta method:
#' Meeker, William Q., and Luis A. Escobar. Statistical methods for reliability data. John Wiley & Sons, 2014. (Appendix B.2)
#'
#' @export

add_ci.survreg <- function(tb, fit,
                           alpha = 0.1,
                           names = NULL,
                           yhatName = "mean_pred",
                           method = "parametric", nSims = 2000,
                           ...){

    if (is.null(names)){
        names[1] <- paste("LCB", alpha/2, sep = "")
        names[2] <- paste("UCB", 1 - alpha/2, sep = "")
    }

    if ((names[1] %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe. Overwriting.")
    }

    if (!(fit$dist %in%
          c("loglogistic", "lognormal", "loggaussian", "exponential", "weibull")))
        stop("Unsupported distribution")

    if (!is.null(fit$weights))
        if (var(fit$weights) != 0)
            stop("weighted regression is unsupported.")

    if(method == "boot")
        stop("not yet implemented")
    ## boot_ci_survreg_expectation(tb, fit, yhatName,
    ##                             confint, alpha, names, nSims)
    else if(method == "parametric")
        parametric_ci_survreg_expectation(tb, fit,
                                          alpha,
                                          names,
                                          yhatName)
    else
        stop("method must be either 'boot' or 'parametric'")
}

parametric_ci_survreg_expectation <- function(tb, fit,
                                              alpha,
                                              names,
                                              yhatName){
    distr <- fit$dist

    if (distr == "loglogistic" && (scale >= 1))
        stop("Expected value is undefined for loglogistic distribution with scale >= 1")

    form <- formula(fit)
    m <- model.frame(form, tb)
    mat <- model.matrix(form, m)

    if(any(is.na(mat)))
        stop("Check tb for missingness")

    nPred <- dim(tb)[1]
    beta <- coef(fit)
    scale <- fit$scale

    if (distr == "weibull")
        pred <- exp(mat %*% beta) * gamma(1 + scale)
    else if (distr == "exponential")
        pred <- exp(mat %*% beta)
    else if (distr == "loglogistic")
        pred <- exp(mat %*% beta) * gamma(1 + scale) * gamma(1 - scale)
    else if ((distr == "lognormal") || (distr == "loggaussian"))
        pred <- exp(mat %*% beta + (scale^2) / 2)

    crit_val <- qnorm(p = 1 - alpha/2, mean = 0, sd = 1)
    cov_mat <- vcov(fit)

    if (distr == "exponential")
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
