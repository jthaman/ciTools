#' @export

add_probs.glm <- function(tb, fit, quant, name = NULL, comparison = "<",
                          nSims = 200, ...){

    if (is.null(name) && comparison == "<")
        name <- paste("Pr(Y < ", quant, ")", sep="")
    else if (is.null(name) && comparison == ">")
        name <- paste("Pr(Y > ", quant, ")", sep="")
    else if (is.null(name) && comparison == "<=")
        name <- paste("Pr(Y <= ", quant, ")", sep="")
    else if (is.null(name) && comparison == ">=")
        name <- paste("Pr(Y >= ", quant, ")", sep="")
    else if (is.null(name) && comparison == "=")
        name <- paste("Pr(Y = ", quant, ")", sep="")
    else
        stop ("Cannot understand this probability statement")

    if ((name %in% colnames(tb))) {
        warning ("These probabilities may have already been appended to your dataframe")
        return(tb)
    }

    if (fit$family$family == "binomial"){
        warning ("Be careful. You should only be asking probabilities that are equivalent to Pr(Y = 0) or Pr(Y = 1).")
        probs_logistic(tb, fit, quant, name, comparison)
    }

    else if (fit$family$family == "poisson"){
        warning("The response is not continuous, so estimated probabilities are only approximate")
        sim_probs_pois(tb, fit, quant, name, nSims, comparison)
    }

    else
        stop("This family is not supported")
}

probs_logistic <- function(tb, fit, quant, name, comparison, ...){
    inverselink <- fit$family$linkinv
    out <- predict(fit, tb, se.fit = TRUE)
    out <- inverselink(out$fit)
    if (((comparison == "=") && (quant == 0)) || ((comparison == "<") && (quant < 1) && (quant > 0)))
        probs <- 1 - out
    else
        probs <- out
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- out
    tb[[name]] <- probs
    as_data_frame(tb)
}

sim_probs_pois <- function(tb, fit, quant, name, nSims, comparison){
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

    probs <- apply(sim_response, 1, FUN = calc_prob, quant = quant, comparison = comparison)
    
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- out
    tb[[name]] <- probs
    as_data_frame(tb)

}


