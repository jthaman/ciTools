## add_probs method for glm objects

add_probs.glm <- function(tb, fit, quant, probName = NULL, comparison = "<",
                          nSims = 200, ...){

    if (is.null(probName) && comparison == "<")
        probName <- paste("Pr(Y < ", quant, ")", sep="")
    else if (is.null(probName) && comparison == ">")
        probName <- paste("Pr(Y > ", quant, ")", sep="")
    else if (is.null(probName) && comparison == "<=")
        probName <- paste("Pr(Y <= ", quant, ")", sep="")
    else if (is.null(probName) && comparison == ">=")
        probName <- paste("Pr(Y >= ", quant, ")", sep="")
    else if (is.null(probName) && comparison == "=")
        probName <- paste("Pr(Y = ", quant, ")", sep="")
    else
        stop ("Cannot understand this probability statement")

    if ((probName %in% colnames(tb))) {
        warning ("These probabilities may have already been appended to your dataframe")
        return(tb)
    }

    if (fit$family$family == "binomial"){
        warning ("Be careful. You should only be asking for Pr(Y = 0) or Pr(Y = 1).")
        probs_logistic(tb, fit, quant, probName, comparison)
    }

    else if (fit$family$family == "poisson"){
        warning("The response is not continuous, so estimated probabilities are only approximate")
        sim_probs_pois(tb, fit, quant, probName, nSims, comparison)
    }

    else
        stop("This family is not supported")
}

probs_logistic <- function(tb, fit, quant, probName, comparison, ...){
    inverselink <- fit$family$linkinv
    out <- predict(fit, tb, se.fit = TRUE)
    out <- inverselink(out$fit)
    if (((comparison == "=") && (quant == 0)) || ((comparison == "<") && (quant < 1) && (quant > 0)))
        probs <- 1 - out
    else
        probs <- out
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- out
    tb[[probName]] <- probs
    return(tb)
}

calc_prob <- function(x, quant, comparison){
    if (comparison == "<")
        mean(x < quant)
    else if (comparison == ">")
        mean(x > quant)
    else if (comparison == "<=")
        mean(x <= quant)
    else
        mean(x >= quant)
}

sim_probs_pois <- function(tb, fit, quant, probName, nSims, comparison){
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
    tb[[probName]] <- probs
    tb

}


