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

#' Confidence Intervals for Predicted Survival Time Quantiles of
#' Accelerated Failure Time Models
#'
#' This function is one of the methods of \code{add_quantile} and is
#' automatically called when an object of class \code{survreg} is
#' passed to \code{add_quantile}.
#'
#' \code{add_quantile.survreg} produces quantiles for the estimated
#' distribution of survival times from a survreg object. Estimated
#' quantiles (such as the median survival time) may be calculated for
#' a range of distributions including lognormal, exponential, weibull,
#' and loglogistic models. \code{add_quantile.survreg} can compute
#' quantiles through a parametric method based on the Delta
#' Method. Generally, this method performs well under a mild to
#' moderate amount of censoring. Parametric intervals are calculated
#' using a transformation of the confidence intervals produced by
#' \code{predict.survreg} and are mathematically identical to intervals
#' calculated by a manual Delta Method.
#'
#' Unlike other \code{add_quantile} methods,
#' \code{add_quantile.survreg} additionally produces confidence
#' intervals for predicted quantiles by default. This may optionally
#' be disabled by switching the \code{confint} argument.
#'
#' Note: Due to a limitation, the \code{Surv} object must be specified in
#' \code{survreg} function call. See the examples section for one way
#' to do this.
#'
#' Note: \code{add_quantile.survreg} cannot inspect the convergence of
#' \code{fit}. Poor maximum likelihood estimates will result in poor
#' confidence intervals. Inspect any warning messages given from
#' \code{survreg}.
#'
#' @param df A data frame of new data.
#' @param fit An object of class \code{survreg}. Predictions are made
#'     with this object.
#' @param p A real number between 0 and 1. Sets the probability level
#'     of the quantiles.
#' @param name \code{NULL} or a character vector of length 3. If
#'     \code{NULL}, quantiles automatically will be named by
#'     \code{add_quantile}, otherwise, they will be named \code{name}.
#' @param yhatName A string. Name of the vector of predictions.
#' @param confint A logical. If \code{TRUE}, confidence intervals for
#'     the quantiles are also appended to \code{df}.
#' @param alpha A number. Controls the confidence level of the
#'     confidence intervals if \code{confint = TRUE}.
## #' @param method A string. One of either \code{"parametric"} or
## #'     \code{"boot"}.
## #' @param nSims A positive integer. Set the
## #'     number of simulated draws to use.
#' @param ... Additional arguments.
#' @return A dataframe, \code{df}, with predicted medians, level \eqn{p}
#'     quantiles, and confidence intervals attached.
#'
#' @seealso \code{\link{add_ci.survreg}} for confidence intervals
#'     \code{survreg} objects, \code{\link{add_pi.survreg}} for
#'     prediction intervals of \code{survreg} objects, and
#'     \code{\link{add_probs.survreg}} for survival probabilities of
#'     \code{survreg} objects.
#'
#' @references
#' For descriptions of the log-location scale models supported:
#' Meeker, William Q., and Luis A. Escobar. Statistical methods for reliability data. John Wiley & Sons, 2014. (Chapter 4)
#'
#' For a description of the multivariate Delta method:
#' Meeker, William Q., and Luis A. Escobar. Statistical methods for reliability data. John Wiley & Sons, 2014. (Appendix B.2)
#'
#' For a description of Delta Method Confidence Intervals:
#' Meeker, William Q., and Luis A. Escobar. Statistical methods for reliability data. John Wiley & Sons, 2014. (Chapter 8)
#'
#' @examples
#' ## Define a data set:
#' df <- survival::stanford2
#' ## remove a covariate with missing values:
#' df <- df[, 1:4]
#' ## next, create the Surv object inside the survreg call:
#' fit <- survival::survreg(survival::Surv(time, status) ~ age + I(age^2),
#'                          data = df, dist = "lognormal")
#' ## Calculate the level 0.75 quantile wit CIs for that quantile
#' add_quantile(df, fit, p = 0.75, name = c("quant", "lwr", "upr"))
#'
#' ## Try a weibull model for the same data:
#' fit2 <- survival::survreg(survival::Surv(time, status) ~ age + I(age^2),
#'                           data = df, dist = "weibull")
#' ## Calculate the level 0.75 quantile with CIs for the quantile
#' add_quantile(df, fit2, p = 0.75, name = c("quant", "lwr", "upr"))
#'
#' @export

add_quantile.survreg <- function(df, fit, p = 0.5,
                                 name = NULL,
                                 yhatName = "median_pred",
                                 confint = TRUE,
                                 alpha = 0.1,
                                 ...){
    if (p <= 0 || p >= 1)
        stop ("p should be in (0,1)")
    if (is.null(name)){
        name[1] <- paste("quantile", p, sep="")
        name[2] <- paste("lcb")
        name[3] <- paste("ucb")
    }
    if ((name[1] %in% colnames(df))) {
        warning ("These quantiles may have already been appended to your dataframe. Overwriting.")
    }

    if (!(fit$dist %in%
          c("loglogistic", "lognormal", "loggaussian", "exponential", "weibull")))
        stop("Unsupported distribution")

    if (!is.null(fit$weights))
        if (var(fit$weights) != 0)
            stop("weighted regression is unsupported.")

    if(any(is.na(df)))
        stop("Check df for missingness")

    ## if(method == "boot")
    ##     boot_ci_survreg_quantile(df, fit, p, name, yhatName,
    ##                              confint, alpha, nSims)
    ## else if(method == "parametric")
        parametric_ci_survreg_quantile(df, fit, p, name, yhatName,
                                       confint, alpha)
    ## else
    ##     stop("method must be either 'boot' or 'parametric'")
}

boot_ci_survreg_quantile <- function(df, fit, p, name, yhatName,
                                     confint, alpha, nSims){
    nPred <- dim(df)[1]
    out <- predict(fit, df, se.fit = TRUE,
                   type = "quantile", p = p)
    med <- predict(fit, df, se.fit = TRUE,
                   type = "quantile", p = 0.5)
    pred <- out$fit

    if (confint){
        boot_mat <- matrix(NA, nrow = nSims, ncol = nPred)
        for (i in 1:nSims){
            temp <- df[sample(1:nPred, size = nPred, replace = TRUE),]
            boot_fit <- survival::survreg(formula(fit$terms), data = temp,
                                          dist = fit$dist)
            boot_pred <- predict(boot_fit, df,
                                 type = "quantile", p = p)
            boot_mat[i,] <- boot_pred
        }
        lwr = apply(boot_mat, 2, quantile, probs = alpha / 2)
        upr = apply(boot_mat, 2, quantile, probs = 1 - alpha / 2)
    }
    if (is.null(df[[yhatName]]))
        df[[yhatName]] <- med$fit

    df[[name[1]]] <- pred

    if (confint){
        df[[name[2]]] <- lwr
        df[[name[3]]] <- upr
    }
    data.frame(df)
}

parametric_ci_survreg_quantile <- function(df, fit, p, name, yhatName,
                                           confint, alpha){
    out <- predict(fit, df, se.fit = TRUE, type = "quantile", p = p)
    med <- predict(fit, df, se.fit = TRUE, type = "quantile", p = 0.5)
    pred <- out$fit

    if (confint){
        crit_val <- qnorm(p = 1 - alpha/2, mean = 0, sd = 1)
        se <- out$se.fit
        w <- exp(crit_val * se / pred)
        upr <- pred * w
        lwr <- pred / w
    }

    if (is.null(df[[yhatName]]))
        df[[yhatName]] <- med$fit

    df[[name[1]]] <- pred

    if (confint){
        df[[name[2]]] <- lwr
        df[[name[3]]] <- upr
    }
    data.frame(df)
}
