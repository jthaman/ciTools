#' Prediction Intervals for Linear Model Predictions.
#'
#' This function is one of the methods for \code{add_pi}. This
#' function is essentially a wrapper for the \code{predict} function
#' in \code{R}.
#'
#' @param tb A tibble or Data Frame.
#' @param fit An object of class lm. Predictions are made with this
#'     object.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, prediction bounds will automatically be named by
#'     \code{add_pi}, otherwise, the lower prediction bound will be
#'     named \code{names[1]} and the upper prediction bound will be
#'     named \code{names[2]}.
#' @param log_response logical. If TRUE, prediction intervals will be
#'     generated for the prediction made with a log-linear model:
#'     \eqn{\log(Y) = X\beta + \epsilon}. These intervals will be on
#'     the scale of the original response, Y.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#'
#' @examples
#' # linear regression
#' fit1 <- lm(dist ~ speed, data = cars)
#' add_pi.lm(cars, fit1, alpha = 0.5)
#' 
#' @export

add_pi.lm <- function(tb, fit, alpha = 0.05, names = NULL, log_response = FALSE){
    
    if (log_response)
        add_pi_lm_log(tb, fit, alpha, names)

    else {
        if (is.null(names)){
            names[1] <- paste("LPB", alpha/2, sep = "")
            names[2] <- paste("UPB", 1 - alpha/2, sep = "")
        }
        if ((names[1] %in% colnames(tb))) {
            warning ("These PIs may have already been appended to your dataframe")
            return(tb)
        }
        out <- predict(fit, tb, interval = "prediction", level = 1 - alpha)
        if(is.null(tb[["pred"]]))
            tb[["pred"]] <- out[, 1]
        if (is.null(tb[[names[1]]]))
            tb[[names[1]]] <- out[, 2]
        if (is.null(tb[[names[2]]]))
            tb[[names[2]]] <- out[, 3]
        tibble::as_data_frame(tb)
    }
}

