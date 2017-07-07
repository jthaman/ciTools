#' Quantiles for the Response of a Generalized Linear Model
#'
#' This function is one of the methods for
#' \code{add_quantile}. Currently, you can only use this function to
#' compute the quantiles of the response of a Poisson regression with
#' the log link function.
#' 
#' @param tb A tibble or Data Frame.
#' @param fit An object of class lm. Predictions are made with this
#'     object.
#' @param p A real number between 0 and 1. Sets the probability for
#'     which to calculate the quantiles.
#' @param name NULL or character vector of length one. If \code{NULL},
#'     quantiles will automatically be named by \code{add_quantile},
#'     otherwise, they will be named \code{name}
#' @param nSims A positive integer. Set the number of simulated draws
#'     to use.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#'
#' @examples
#' # linear regression
#' fit1 <- lm(dist ~ speed, data = cars)
#' # Calculate the 0.30 quantile of dist | speed.
#' add_quantile.lm(cars, fit1, p = 0.3)
#' 
#' @export

add_quantile.glm <- function(tb, fit, p, name = NULL, nSims = 200){
    if (p <= 0 || p >= 1)
        stop ("p should be in (0,1)")
    if (is.null(name))
        name <- paste("quantile", p, sep="")
    if (name %in% colnames(tb)) {
        warning ("These quantiles may have already been appended to your dataframe")
        return(tb)
    }
    if ((name %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe")
        return(tb)
    }
    if (fit$family$family == "binomial"){
       stop ("Quantiles for Logistic Regression don't make sense") 
    }
    if (fit$family$family == "poisson"){
        warning ("The response is not continuous, so estimated quantiles are only approximate")
        sim_quantile_pois(tb, fit, p, name, nSims)
    }
}

sim_quantile_pois <- function(tb, fit, p, name, nSims){
    nPreds <- NROW(tb)
    modmat <- model.matrix(fit)
    response_distr <- fit$family$family
    inverselink <- fit$family$linkinv
    out <- predict(fit, newdata = tb, type = "response")
    sims <- arm::sim(fit, n.sims = nSims)
    sim_response <- matrix(0, ncol = nSims, nrow = nPreds)

    for (i in 1:nPreds){
        if(response_distr == "poisson"){
            sim_response[i,] <- rpois(n = nSims, lambda = inverselink(rnorm(nPreds,sims@coef[i,] %*% modmat[i,], sd = sims@sigma[i])))
            }
    }

    quants <- apply(sim_response, 1, FUN = quantile, probs = p, type = 1)

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- out
    tb[[name]] <- quants
    as_data_frame(tb)


}
