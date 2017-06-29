## add_pi method for glm objects
add_pi.glm <- function(tb, fit, alpha = 0.05, piNames = NULL,
                       type = "response", piType = "sim", nSims = 1000){

    if (is.null(piNames)){
        piNames[1] <- paste("LPB", alpha/2, sep = "")
        piNames[2] <- paste("UPB", 1 - alpha/2, sep = "")
    }
    if ((piNames[1] %in% colnames(tb))) {
        warning ("These PIs may have already been appended to your dataframe")
        return(tb)
    }

    if(fit$family$family == "binomial"){
        stop("Prediction interval for Bernoulli response doesn't make sense")
    }

    if(fit$family$family == "poisson"){
        warning("The response is not continuous, so Prediction Intervals are only approximate")
    }

    if(piType == "sim"){
        sim_pi_glm(tb, fit, alpha, piNames, type, nSims)
    }
    else if(!(method %in% c("sim")))
        stop("Only Simulated prediction intervals are implemented for glm objects")
}

## TODO : hardcode more response distributions

sim_pi_glm <- function(tb, fit, alpha, piNames, type, nSims){
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

    lwr <- apply(sim_response, 1, FUN = quantile, probs = alpha / 2)
    upr <- apply(sim_response, 1, FUN = quantile, probs = 1 - alpha / 2)
    
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- out
    tb[[piNames[1]]] <- lwr
    tb[[piNames[2]]] <- upr
    tb

}

