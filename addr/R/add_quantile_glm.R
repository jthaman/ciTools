#' @export

add_quantile.glm <- function(tb, fit, prob, quantileName = NULL, nSims = 200){
    if (prob <= 0 || prob >= 1)
        stop ("prob should be in (0,1)")
    if (is.null(quantileName))
        quantileName <- paste("quantile", prob, sep="")
    if (quantileName %in% colnames(tb)) {
        warning ("These quantiles may have already been appended to your dataframe")
        return(tb)
    }
    if ((quantileName %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe")
        return(tb)
    }
    if (fit$family$family == "binomial"){
       stop ("Quantiles for Logistic Regression don't make sense") 
    }
    if (fit$family$family == "poisson"){
        warning ("The response is not continuous, so estimated quantiles are only approximate")
        sim_quantile_pois(tb, fit, prob, quantileName, nSims)
    }
}

sim_quantile_pois <- function(tb, fit, prob, quantileName, nSims){
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

    quants <- apply(sim_response, 1, FUN = quantile, probs = prob)

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- out
    tb[[quantileName]] <- quants
    as_data_frame(tb)


}
