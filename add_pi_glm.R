## add_pi method for glm objects
add_pi.glm <- function(tb, fit, alpha = 0.05, piNames = c("LPB", "UPB"),
                       type = "response", method = "sim", nSims = 1000){
    if(fit$family$family == "binomial"){
        stop("Prediction interval for Binomial response doesn't make sense")
    }
    if(ciType == "sim"){
        sim_pi_glm(tb, fit, alpha, piNames, type, nSims, ...)
    }
    else if(!(method %in% c("sim")))
        stop("Only Simulated prediction intervals are implemented for glm objects")
}

## add_pi simulation method for glm
## Generate fake data from the model to calculate prediction intervals
## TODO : hardcode more response distributions

sim_pi_glm <- function(tb, fit, alpha, piNames, type, nSims, ...){
    nPreds <- NROW(tb)
    modmat <- model.matrix(fit)
    response_distr <- fit$family$family
    inverselink <- fit$family$linkinv
    out <- inverselink(predict(fit, tb))
    fitted_values <- fit$family$linkinv
    sims <- arm::sim(fit, n.sims = nSims)
    
    sim_response <- matrix(0, ncol = nSims, nrow = nPreds)

    ## TODO Fix this up, eliminate for loop

    for (i in 1:nPreds){
        if(response_distr == "poisson"){
            sim_response[i,] <- rpois(n = nSims, lambda = inverselink(rnorm(nPreds,sims@coef[i,] %*% modmat[i,], sd = sims@sigma[i])))
            }
    }

    ## get the reponse level intervals
    lwr <- apply(sim_response, 1, FUN = quantile, probs = alpha / 2)
    upr <- apply(sim_response, 1, FUN = quantile, probs = 1 - alpha / 2)
    
    ## bind to the tibble and return
    if(is.null(tb[["pred"]])) tb[["pred"]] <- out
    tb[[piNames[1]]] <- lwr
    tb[[piNames[2]]] <- upr
    tb

}

