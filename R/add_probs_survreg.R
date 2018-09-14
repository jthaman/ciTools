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

#' Confidence Intervals for the Survivor Function of Parametric Survival Models
#'
#' TODO add a description
#'
#' TODO add details
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{survreg}. Predictions are made
#'     with this object.
#' @param q A double. A quantile of the response distribution. In
#'     survival applications this is the time of event.
#' @param confint A logical. If \code{TRUE}, confidence intervals for
#'     the estimated probabilities will be calculated and appended to
#'     \code{tb}.
#' @param alpha A number. Control the confidence level of the
#'     confidence intervals if \code{confint = TRUE}.
#' @param name \code{NULL} or a string. If \code{NULL}, probabilities
#'     automatically will be named by \code{add_probs()}, otherwise,
#'     the probabilities will be named \code{name} in the returned
#'     tibble
#' @param method A string. One of either \code{"parametric"} or
#'     \code{"bootstrap"}.
#' @param yhatName A string. Name of the vector of predictions.
#' @param comparison A character vector of length one. If
#'     \code{comparison = "<"}, then \eqn{Pr(Y|X < q)} is
#'     calculated. If \code{comparison = ">"}, the survivor function
#'     at time \code{q} is calculated.
#' @param nSims Number of simulations used for bootstrap confidence
#'     intervals.
#' @param ... Additional arguments.
#'
#' @return A tibble, \code{tb}, with predicted values and
#'     probabilities attached.
#'
#' @seealso \code{\link{add_ci.survreg}} for confidence intervals for
#'     \code{survreg} objects, \code{\link{add_pi.survreg}} for
#'     prediction intervals of \code{survreg} objects, and
#'     \code{\link{add_quantile.survreg}} for response quantiles of
#'     \code{survreg} objects.
#'
#' @examples
#'
#' @export


add_probs.survreg <- function(tb, fit, q,
                              name = NULL, yhatName = "median_pred",
                              comparison = "<",
                              method = "parametric",
                              confint = TRUE,
                              alpha = 0.05,
                              nSims = 2000, ...){

    if (is.null(name) & (comparison == "<" || comparison == ">=")){

        name[1] <- paste("prob_less_than", q, sep="")
        name[2] <- paste("lcb")
        name[3] <- paste("ucb")
    }
    if (is.null(name) & (comparison == ">" || comparison == ">=")){

        name[1] <- paste("prob_greater_than", q, sep="")
        name[2] <- paste("lcb")
        name[3] <- paste("ucb")
    }
    if ((name[1] %in% colnames(tb)))
        warning ("These probabilities may have already been appended to your dataframe. Overwriting.")

    if (method == "parametric") {
        parametric_ci_survreg_prob(tb, fit, q, name, yhatName, comparison,
                                   confint, alpha)
    }
    else if (method == "boot") {
        boot_ci_survreg_prob(tb, fit, q, name, yhatName, comparison,
                             confint, alpha, nSims)
    }
    else
        stop("method must be either boot or parametric")
}

survreg_calc_probs <- function(tb, fit, q, comparison){
    form <- formula(fit)
    m <- model.frame(form, tb)
    mat <- model.matrix(form, m)

    dist <- survival::survreg.distributions[[fit$dist]][["dist"]]
    fn_list <- survival::survreg.distributions[[dist]][["density"]]
    cdf <- function(x) fn_list(x)[,1]

    pred <- predict(fit, tb, type = "linear")
    scale <- fit$scale
    zeta <- (log(q) - pred) / scale

    if (comparison == "<" || comparison == "<=")
        F <- cdf(zeta)
    else if (comparison == ">" || comparison == ">=")
        F <- 1 - cdf(zeta)
    else
        stop("invalid comparison")

    return(list(
        F = F,
        mat = mat,
        zeta = zeta
    ))
}

parametric_ci_survreg_prob <- function(tb, fit, q, confint,
                                       alpha, name, yhatName,
                                       comparison){

    collect <- survreg_calc_probs(tb = tb, fit = fit, q = q,
                                  comparison = comparison)

    F <- collect[["F"]]
    zeta <- collect[["zeta"]]
    mat <- collect[["mat"]]

    if (confint){
        dist <- survival::survreg.distributions[[fit$dist]][["dist"]]
        fn_list <- survival::survreg.distributions[[dist]][["density"]]
        dens <- function(x) fn_list(x)[,3]
        scale <- fit$scale
        crit_val <- qnorm(p = 1 - alpha/2, mean = 0, sd = 1)
        cov_mat <- vcov(fit)

        if (fit$dist == "exponential")
            cov_mat <- cbind(rbind(cov_mat, 0), 0)

        d_g <- rep(NA, dim(mat)[1])
        seF <- rep(NA, dim(mat)[1])

        for(i in 1:dim(mat)[1]){
            d_g_beta <- dens(zeta[i]) * (-mat[i,] / scale)
            d_g_delta <- dens(zeta[i]) * (-zeta[i])
            d_g_vec <- c(d_g_beta, d_g_delta)
            seF[i] <- sqrt(t(d_g_vec) %*% cov_mat %*% d_g_vec)
        }

        w <- exp(crit_val * seF / (F * (1 - F)))
        lwr <- F / (F + (1 - F) * w)
        upr <- F / (F + (1 - F) / w)
    }

    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- predict(fit, tb, type = "quantile" , p = 0.5)

    tb[[name[1]]] <- as.numeric(F)

    if (confint){
        tb[[name[2]]] <- as.numeric(lwr)
        tb[[name[3]]] <- as.numeric(upr)
    }

    tibble::as_data_frame(tb)
}

boot_ci_survreg_prob <- function(tb, fit, q, confint,
                                 alpha, name, yhatName,
                                 comparison, nSims){
    pred <- predict(fit, tb, type = "quantile", p = 0.5)

    collect <- survreg_calc_probs(tb = tb, fit = fit, q = q,
                                  comparison = comparison)

    F <- collect[["F"]]

    if (confint){
        nPred <- NROW(tb)
        boot_mat <- matrix(NA, nrow = nSims, ncol = nPred)
        for (i in 1:nSims){
            temp <- tb[sample(1:nPred, size = nPred, replace = TRUE),]
            boot_fit <- survival::survreg(formula(fit$terms), data = temp,
                                          dist = fit$dist)
            boot_mat[i,] <- survreg_calc_probs(tb = tb, fit = boot_fit, q = q,
                                               comparison = comparison)[["F"]]
        }
        lwr = apply(boot_mat, 2, quantile, probs = alpha / 2)
        upr = apply(boot_mat, 2, quantile, probs = 1 - alpha / 2)
    }

    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- pred

    if (confint){
        tb[[name[2]]] <- lwr
        tb[[name[3]]] <- upr
    }
    tb[[name[1]]] <- F

    tibble::as_data_frame(tb)
}
