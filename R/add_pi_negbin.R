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

#' Prediction Intervals for Negative Binomial Linear Models
#'
#' This function is one of the methods for \code{add_pi}, and is
#' called automatically when \code{add_pi} is used on a \code{fit} of
#' class \code{negbin}. 
#'
#' Prediction intervals for negative binomial fits are formed through
#' a two part simulation scheme:
#'
#' 1. Model coefficients are generated through a parametric bootstrap
#' procedure that simulates the uncertainty in the regression
#' coefficients.
#'
#' 2. Random draws from the negative binomial
#' distribution are taken with a mean that varies based on the model
#' coefficients determined in step (1) and over-dispersion parameter
#' that is taken from the original fitted model.
#'
#' Quantiles of the simulated responses are taken at the end to
#' produce intervals of the desired level.
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{negbin}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, prediction bounds automatically will be named by
#'     \code{add_pi}, otherwise, the lower prediction bound will be
#'     named \code{names[1]} and the upper prediction bound will be
#'     named \code{names[2]}.
#' @param yhatName A string. Name of the predictions vector.
#' @param nSims A positive integer. Determines the number of
#'     simulations to run.
#' @param ... Additional arguments.
#' 
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#'
#' @seealso \code{\link{add_ci.negbin}} for confidence intervals for
#'     \code{negbin} objects, \code{\link{add_probs.negbin}} for conditional
#'     probabilities of \code{negbin} objects, and
#'     \code{\link{add_quantile.negbin}} for response quantiles of
#'     \code{negbin} objects.
#'
#' @examples
#' x1 <- rnorm(100, mean = 1)
#' y <- MASS::rnegbin(n = 100, mu = exp(1 + x1), theta = 5)
#' df <- data.frame(x1 = x1, y = y)
#' fit <- MASS::glm.nb(y ~ x1, data = df)
#' add_pi(df, fit, names = c("lpb", "upb"))
#' 
#' @export


add_pi.negbin <- function(tb, fit, alpha = 0.05, names = NULL, yhatName = "pred", 
                          nSims = 2000, ...){

    if (is.null(names)) {
        names[1] <- paste("LPB", alpha/2, sep = "")
        names[2] <- paste("UPB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(tb))) 
        warning ("These PIs may have already been appended to your dataframe. Overwriting.")
    
    sim_pi_negbin(tb, fit, alpha, names, yhatName, nSims)
}


sim_pi_negbin <- function(tb, fit, alpha, names, yhatName, nSims){
    out <- predict(fit, newdata = tb, type = "response")
    disp <- fit$dispersion
    sim_response <- get_sim_response_nb(tb, fit, nSims)

    lwr <- apply(sim_response, 1, FUN = quantile, probs = alpha/2, type = 1)
    upr <- apply(sim_response, 1, FUN = quantile, probs = 1 - alpha / 2, type = 1)
    
    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[names[1]]] <- lwr
    tb[[names[2]]] <- upr
    tibble::as_data_frame(tb)

}

## modification of arm::sim that will accept an negative binomial fit.
get_negbin_sims <- function(fit, nSims) {
    summ <- summary (fit, correlation=TRUE, dispersion = fit$dispersion)
    coef <- summ$coef[,1:2,drop=FALSE]
    dimnames(coef)[[2]] <- c("coef.est","coef.sd")
    beta.hat <- coef[,1,drop=FALSE]
    sd.beta <- coef[,2,drop=FALSE]
    corr.beta <- summ$corr
    n <- summ$df[1] + summ$df[2]
    k <- summ$df[1]
    V.beta <- corr.beta * array(sd.beta,c(k,k)) * t(array(sd.beta,c(k,k)))
    beta <- array (NA, c(nSims,k))
    dimnames(beta) <- list (NULL, dimnames(beta.hat)[[1]])

    for (s in 1:nSims){
        beta[s,] <- MASS::mvrnorm (1, beta.hat, V.beta)
    }

    beta2 <- array (0, c(nSims,length(coefficients(fit))))
    dimnames(beta2) <- list (NULL, names(coefficients(fit)))
    beta2[,dimnames(beta2)[[2]]%in%dimnames(beta)[[2]]] <- beta
    sigma <- rep (sqrt(summ$dispersion), nSims)

    ans <- new("sim",
               coef = beta2,
               sigma = sigma)
    return(ans)
}

## link get_sim_response, but only for negative binomial model
get_sim_response_nb <- function(tb, fit, nSims){
    nPreds <- NROW(tb)
    modmat <- model.matrix(fit, data = tb)
    response_distr <- fit$family$family
    inverselink <- fit$family$linkinv
    sims <- get_negbin_sims(fit, nSims)
    sim_response <- matrix(0, ncol = nSims, nrow = nPreds)

    for (i in 1:nSims){
        yhat <- inverselink(modmat %*% sims@coef[i,])
        sim_response[,i] <- rnegbin(n = nPreds,
                                    mu = yhat,
                                    theta = fit$theta)
    }
    sim_response
}
