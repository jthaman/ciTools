## TODO add pi method for GLM objects (use arm::sim)
## TODO warn about logistic or probit regression
library(arm)

add_pi.glm <- function(tb, fit, alpha = 0.05, piNames = c("LPB", "UPB"),
                       type = "response", method = "sim", nSims = 1000){
    if(fit$family$family == "binomial"){
        stop("Prediction interval for Binomial response is not supported")
    }
    if(ciType == "sim"){
        sim_pi_glm(tb, fit, alpha, piNames, type, nSims, ...)
    }
    ) else
          if(!(method %in% c("sim"))) stop("Only Simulated prediction intervals are implemented for glm objects")
}

## add_ci Wald method for glm
## TODO : hardcode more response distributions

sim_pi_glm <- function(tb, fit, alpha, piNames, type, nSims, ...){
    n_preds <- NROW(tb)
    modmat <- model.matrix(fit)
    response_distr <- fit$family$family
    inverselink <- fit$family$linkinv
    out <- predict(fit, tb, se.fit = TRUE)
    fitted_values <- fit$family$linkinv
    sims <- arm::sim(fit, n.sims = nSims)
    
    sim_response <- matrix(0, ncol = nSims, nrow = n_preds)

    ## TODO Fix this up, eliminate for loop
    for (i in 1:n_preds){
        if(response_distr == "poisson"){
            sim_response[i,] <- rpois(n = nSims, lambda = exp(rnorm(n_preds,sims@coef[i,] %*% modmat[i,], sd = sims@sigma[i])))
            }
    }

    ## get the reponse level intervals
    lwr <- apply(sim_response, 1, FUN = quantile, probs = alpha / 2)
    upr <- apply(sim_response, 1, FUN = quantile, probs = 1 - alpha / 2)
    
    
    ## bind to the tibble and return
    if(is.null(tb[["pred"]])) tb[["pred"]] <- out$fit
    tb[[piNames[1]]] <- lwr
    tb[[piNames[2]]] <- upr
    tb

}

