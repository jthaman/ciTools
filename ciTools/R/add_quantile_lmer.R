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

#' Quantiles for the Response of a Linear Mixed Model
#'
#' This function is one of the methods for \code{add_quantile} and it
#' automatically called when \code{add_quantile} is called on an
#' object of class \code{lmerMod}.
#'
#' \code{add_qauntile} may use three different for determining
#' quantiles: a parametric method, a simulation method (via
#' \code{merTools::predictInterval}), or a parametric bootstrap method
#' (via \code{lme4::simulate}). The default and recommened method is
#' parametric bootstrap, which corresponds to setting \code{type =
#' "boot"}.
#' 
#' @param tb A tibble or Data Frame.
#' @param fit An object of class lm. Predictions are made with this
#'     object.
#' @param p A real number between 0 and 1. Sets the probability for
#'     which to calculate the quantiles.
#' @param name NULL or character vector of length one. If \code{NULL},
#'     quantiles will automatically be named by \code{add_quantile},
#'     otherwise, they will be named \code{name}
#' @param includeRanef Should the quantiles be calculated condition on
#'     the random effects?
#' @param type A character vector of length one. Options are
#'     \code{"parametric"} or \code{"sim"}
#' @param nSims A positive integer. Set the number of simulations to
#'     perform.
#' @param log_response A logical. Set to \code{TRUE} if the model is a
#'     log-linear mixed model.
#' @param yhatName A string. Determines the name of column of
#'     predictions.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#' 
#' @seealso \code{{\link{add_ci.lmerMod}}} for confidence intervals
#'     for \code{lmerMod} objects. \code{\link{add_pi.lmerMod}} for
#'     prediction intervals of \code{lmerMod} objects, and
#'     \code{\link{add_probs.lmerMod}} for response probabilities of
#'     \code{lmerMod} objects.
#'
#' @examples
#' dat <- lme4::sleepstudy
#' fit <- lme4::lmer(Reaction ~ Days + (1|Subject), data = lme4::sleepstudy)
#' add_quantile(dat, fit, p = 0.5)
#' add_quantile(dat, fit, p = 0.5, type = "parametric", includeRanef = FALSE)
#' add_quantile(dat, fit, p = 0.5, type = "sim", name = "my_quantile", nSims = 1000)
#' 
#' @export

add_quantile.lmerMod <- function(tb, fit, 
                                 p, name = NULL, includeRanef = TRUE,
                                 type = "boot",
                                 nSims = 200, log_response = FALSE,
                                 yhatName = "pred", ...) {

    if (p <= 0 || p >= 1)
        stop ("p should be in (0,1)")
    if (is.null(name)){
        name <- paste("quantile", p, sep = "")
    }
    if ((name %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe. Overwriting.")
    }
    if (type == "parametric")
        parametric_quantile_mermod(tb, fit, p, name, includeRanef, log_response, yhatName)
    else if (type == "boot")
        boot_quantile_mermod(tb, fit, p, name, includeRanef, nSims, log_response, yhatName)
    else if (type == "sim")
        sim_quantile_mermod(tb, fit, p, name, includeRanef, nSims, log_response, yhatName)
    else
        stop ("Incorrect type specified")
 }

parametric_quantile_mermod <- function(tb, fit, p, name, includeRanef, log_response, yhatName){
    
    rdf <- get_resid_df_mermod(fit)
    seGlobal <- get_pi_mermod_var(tb, fit, includeRanef)
    
    if(includeRanef)
        re.form <- NULL
    else
        re.form <- NA

    out <- predict(fit, tb, re.form = re.form)
    quant <- out + qt(p ,df = rdf) * seGlobal
    if (log_response)
        out <- exp(out)
    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[name]] <- quant
    if(log_response)
        tb[[name]]<- exp(quant)
    tibble::as_data_frame(tb)
}


boot_quantile_mermod <- function(tb, fit, p, name, includeRanef, nSims, log_response, yhatName) {

    if (includeRanef) 
        reform = NULL
    else 
        reform = NA

    gg <- simulate(fit, re.form = reform, nsim = nSims)
    gg <- as.matrix(gg)
    quant <- apply(gg, 1, FUN = quantile, probs = p)

    out <- predict(fit, tb, re.form = reform)
    if (log_response)
        out <- exp(out)
    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[name]] <- quant
    if (log_response)
        tb[[name]]<- exp(quant)
    tibble::as_data_frame(tb)
}


sim_quantile_mermod <- function(tb, fit, p, name, includeRanef, nSims, log_response, yhatName) {

    if (includeRanef) {
        which = "full"
        reform = NULL
    } else {
        which = "fixed"
        reform = NA
    }

    pi <- suppressWarnings(predictInterval(fit, tb, which = which, level = 0.95,
                              n.sims = nSims,
                              stat = "median",
                              include.resid.var = TRUE,
                              returnSims = TRUE))

    store_sim <- attributes(pi)$sim.results
    quant <- apply(store_sim, 1, FUN = quantile, p = p)
    out <- predict(fit, tb, re.form = reform)

    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[name]] <- quant

    if (log_response){
        tb[[yhatName]] <- exp(out)
        tb[[name]] <- exp(quant)
    }

    tibble::as_data_frame(tb)
    
}
