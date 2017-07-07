#' Confidence Intervals for Linear Model Predictions.
#'
#' This function is one of the methods for \code{add_ci}. 
#'
#' @param tb A tibble or Data Frame.
#' @param fit An object of class lm, glm, or lmerMod. Predictions are
#'     made with this object.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names NULL or character vector of length two. If
#'     \code{NULL}, confidence bounds will automatically be named by
#'     \code{add_ci}, otherwise, the lower confidence bound will be
#'     named \code{names[1]} and the upper confidence bound will be
#'     named \code{names[2]}.
#' @param log_response logical. If TRUE, confidence intervals will be
#'     generated for the prediction made with a log-linear model:
#'     \eqn{\log(Y) = X\beta + \epsilon}
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     confidence bounds attached.
#'
#' @examples
#' # linear regression
#' fit1 <- lm(dist ~ speed, data = cars)
#' add_ci(cars, fit1, alpha = 0.5)
#' 
#' @export

add_ci.lm <- function(tb, fit, alpha = 0.05, names = NULL, log_response = FALSE){
    if (log_response)
        add_ci_lm_log(tb, fit, alpha, names)
    else {

        if (is.null(names)){
            names[1] <- paste("LCB", alpha/2, sep = "")
            names[2] <- paste("UCB", 1 - alpha/2, sep = "")
        }
        if ((names[1] %in% colnames(tb))) {
            warning ("These CIs may have already been appended to your dataframe")
            return(tb)
        }
        out <- predict(fit, tb, interval = "confidence", level = 1 - alpha)
        if(is.null(tb[["pred"]]))
            tb[["pred"]] <- out[, 1]
        if (is.null(tb[[names[1]]]))
            tb[[names[1]]] <- out[, 2]
        if (is.null(tb[[names[2]]]))
            tb[[names[2]]] <- out[, 3]
        as_data_frame(tb)
    } 
}

