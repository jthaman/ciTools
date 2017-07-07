#' Event Probabilities for Generalized Linear Models
#'
#' This is the method \code{add_probs} uses if the model fit is an
#' object of class \code{glm}. Probabilities are determined through
#' simulation, using the same method as \code{add_pi.glm}. Currently,
#' only logistic and poisson regression are supported.
#' 
#' @param tb A tibble or Data Frame on which to append probabilities
#' @param fit An object of class \code{glm}. Predictions are
#'     made with this object.
#' @param q A double. A quantile of the response distribution.
#' @param name NULL or character vector of length one. If \code{NULL},
#'     probabilities will automatically be named by
#'     \code{add_probs()}, otherwise, the probabilities will be named
#'     \code{name} in the returned tibble
#' @param comparison A character vector of length one. If
#'     \code{comparison = "<"}, then Pr(Y|X < q) is
#'     calculated. Must be "<" or ">" for linear, log-linear and
#'     linear mixed models. If \code{fit} is a glm, then
#'     \code{comparison} may also be "<=", ">=", or "=".
#' @param nSims A positive integer. 
#' @return A tibble, \code{tb}, with predicted values and
#'     probabilities attached.
#' 
#' @examples
#' fit1 <- glm(dist ~ speed, data = cars, family = "poisson")
#' add_probs.glm(cars, fit1, q = 40, comparison = "<")
#' add_probs.glm(cars, fit1, q = 40, comparison = ">=")
#' 
#' @export

add_probs.glm <- function(tb, fit, q, name = NULL, comparison = "<",
                          nSims = 200, ...){

    if (is.null(name) && comparison == "<")
        name <- paste("prob_less_than", q, sep="")
    else if (is.null(name) && comparison == ">")
        name <- paste("prob_greater_than", q, sep="")
    else if (is.null(name) && comparison == "<=")
        name <- paste("prob_less_than_or_equal", q, sep="")
    else if (is.null(name) && comparison == ">=")
        name <- paste("prob_greater_than_or_equal", q, sep="")
    else if (is.null(name) && comparison == "=")
        name <- paste("prob_equal_to", q, sep="")
    else
        stop ("Cannot understand this probability statement")

    if ((name %in% colnames(tb))) {
        warning ("These probabilities may have already been appended to your dataframe")
        return(tb)
    }

    if (fit$family$family == "binomial"){
        warning ("Be careful. You should only be asking probabilities that are equivalent to Pr(Y = 0) or Pr(Y = 1).")
        probs_logistic(tb, fit, q, name, comparison)
    }

    else if (fit$family$family == "poisson"){
        warning("The response is not continuous, so estimated probabilities are only approximate")
        sim_probs_pois(tb, fit, q, name, nSims, comparison)
    }

    else
        stop("This family is not supported")
}

probs_logistic <- function(tb, fit, q, name, comparison, ...){
    inverselink <- fit$family$linkinv
    out <- predict(fit, tb, se.fit = TRUE)
    out <- inverselink(out$fit)
    if (((comparison == "=") && (q == 0)) || ((comparison == "<") && (q < 1) && (q > 0)))
        probs <- 1 - out
    else
        probs <- out
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- out
    tb[[name]] <- probs
    as_data_frame(tb)
}

sim_probs_pois <- function(tb, fit, q, name, nSims, comparison){
    nPreds <- NROW(tb)
    modmat <- model.matrix(fit)
    response_distr <- fit$family$family
    inverselink <- fit$family$linkinv
    out <- inverselink(predict(fit, tb))
    fitted_values <- fit$family$linkinv
    sims <- arm::sim(fit, n.sims = nSims)
    sim_response <- matrix(0, ncol = nSims, nrow = nPreds)

    for (i in 1:nPreds){
        if(response_distr == "poisson"){
            sim_response[i,] <- rpois(n = nSims, lambda = inverselink(rnorm(nPreds,sims@coef[i,] %*% modmat[i,], sd = sims@sigma[i])))
            }
    }

    probs <- apply(sim_response, 1, FUN = calc_prob, quant = q, comparison = comparison)
    
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- out
    tb[[name]] <- probs
    as_data_frame(tb)

}


