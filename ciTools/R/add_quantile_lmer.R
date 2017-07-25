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
#' This function is one of the methods for
#' \code{add_quantile}. 
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
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#' 
#' @export

add_quantile.lmerMod <- function(tb, fit, 
                                 p, name = NULL, includeRanef = TRUE,
                                 type = "parametric",
                                 nSims = 200, ...) {

    if (p <= 0 || p >= 1)
        stop ("p should be in (0,1)")
    if (is.null(name)){
        name <- paste("quantile", p, sep = "")
    }
    if ((name %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe. Overwriting.")
    }
    if (type == "parametric")
        parametric_quantile_mermod(tb, fit, p, name, includeRanef)
    else if (type == "sim_lme4")
        simulate_quantile_mermod(tb, fit, p, name, includeRanef, nSims)
    else
        stop ("Incorrect type specified")
 }

parametric_quantile_mermod <- function(tb, fit, p, name, includeRanef){
    
    rdf <- get_resid_df_mermod(fit)
    seGlobal <- get_pi_mermod_var(tb, fit, includeRanef)
    
    if(includeRanef)
        re.form <- NULL
    else
        re.form <- NA

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- predict(fit, tb, re.form = re.form)
    tb[[name]] <- tb[["pred"]] + qt(p ,df = rdf) * seGlobal
    tibble::as_data_frame(tb)
}


simulate_quantile_mermod <- function(tb, fit, p, name, includeRanef, nSims = 200) {

    if (includeRanef) 
        reform = NULL
    else 
        reform = NA

    gg <- simulate(fit, re.form = reform, nsim = nSims)
    gg <- as.matrix(gg)
    quant <- apply(gg, 1, FUN = quantile, probs = p)

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- predict(fit, tb, re.form = reform)
    tb[[name]] <- quant
    tibble::as_data_frame(tb)
}
