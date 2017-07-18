#' Quantiles for the Response of a Linear Model
#'
#' This function is one of the methods for \code{add_quantile}.
#' 
#' @param tb A tibble or Data Frame.
#' @param fit An object of class lm. Predictions are made with this
#'     object.
#' @param p A real number between 0 and 1. Sets the probability for
#'     which to calculate the quantiles.
#' @param name NULL or character vector of length one. If \code{NULL},
#'     quantiles will automatically be named by \code{add_quantile},
#'     otherwise, they will be named \code{name}
#' @param log_response logical. If TRUE, quantiles will be generated
#'     for the prediction made with a log-linear model: \eqn{\log(Y) =
#'     X\beta + \epsilon}. These quantiles will be on the scale of the
#'     original response, Y.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#'
#' @examples
#' # linear regression
#' fit1 <- lm(dist ~ speed, data = cars)
#' add_quantile.lm(cars, fit1, alpha = 0.5)
#' 
#' @export

add_quantile.lm <- function(tb, fit, p, name = NULL, log_response = FALSE){
    if (log_response)
        add_quantile_lm_log(tb, fit, p, name)
    else {
        if (p == 0.5)
            warning ("The 0.5 quantile is equal to the fitted values")
        if (p <= 0 || p >= 1)
            stop ("p should be in (0,1)")
        if (is.null(name))
            name <- paste("quantile", p, sep="")
        if (name %in% colnames(tb)) {
            warning ("These quantiles may have already been appended to your dataframe")
            return(tb)
        }
        out <- predict(fit, tb, interval = "prediction", se.fit = TRUE)
        fitted <- out$fit[,1]
        residual_df <- out$df
        se_fitted <- out$se.fit
        resid_var <- out$residual.scale^2
        se_pred <- sqrt(resid_var + se_fitted^2)
        t_quantile <- qt(p = p, df = residual_df)
        out_quantiles <- fitted + se_pred * t_quantile
        if (is.null(tb[["pred"]]))
            tb[["pred"]] <- fitted
        tb[[name]] <- out_quantiles
        tibble::as_data_frame(tb)
    }
}

