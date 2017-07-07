#' @export

add_quantile.glm <- function(tb, fit, prob, name = NULL, nSims = 200){
    if (prob <= 0 || prob >= 1)
        stop ("prob should be in (0,1)")
    if (is.null(name))
        name <- paste("quantile", prob, sep="")
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
        sim_quantile_pois(tb, fit, prob, name, nSims)
    }
}

sim_quantile_pois <- function(tb, fit, prob, name, nSims){
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

    quants <- apply(sim_response, 1, FUN = quantile, probs = prob, type = 1)

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- out
    tb[[name]] <- quants
    as_data_frame(tb)


}
