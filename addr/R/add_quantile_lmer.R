#' Quantiles for the Response of a Linear Mixed Model
#'
#' This function is one of the methods for
#' \code{add_quantile}. Currently, only a parametric method is
#' supported.
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
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#' 
#' @export

add_quantile.lmerMod <- function(tb, fit, 
                          p, includeRanef = TRUE,
                          name = NULL, ...) {

    if (p <= 0 || p >= 1)
        stop ("p should be in (0,1)")
    if (is.null(name)){
        name <- paste("quantile", p, sep = "")
    }
    if ((name %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe")
        return(tb)
    }
    parametric_quantile_mermod(tb, fit, p, name, includeRanef, ...)
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
    as_data_frame(tb)
}


