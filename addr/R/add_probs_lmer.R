## TODO: add_probs method for lmer objects
add_probs.merMod <- function(tb, fit, 
                          quant, probsType = "parametric", 
                          includeRanef = T, probName = NULL, ...) {
    if(probsType == "bootstrap") {
        stop ("this Type is not yet implemented")
        ##bootstrap_probs_mermod(tb, fit, quant, probName, ...)
    } else if (probType == "parametric") {
        parametric_probs_mermod(tb, fit, quant, probName, includeRanef)
    } else if (probType == "simulation") {
        simulation_probs_mermod(tb, fit, quant, probName, includeRanef, ...)
    } else if (!(probType %in% c("bootstrap", "parametric", "simulation"))) {
        stop("Incorrect type specified!")
    }
}

## could make this the default procedure with warnings
parametric_probs_merMod <- function(tb, fit, quant, probName, includeRanef){
    if (is.null(probName))
        probName <- paste("Pr(Y < ", quant, ")", sep="")
    if ((probName %in% colnames(tb))) {
        warning ("These probabilities may have already been appended to your dataframe")
        return(tb)
    }

    X <- model.matrix(reformulate(attributes(terms(fit))$term.labels), tb)
    vcovBetaHat <- vcov(fit)
    
    seFixed <- X %*% vcovBetaHat %*% t(X) %>% 
        diag() %>%
        sqrt()

    fitted <- predict(fit, newdata = tb)
    seRandom <- arm::se.ranef(fit)[[1]][1,]
    rdf <- nrow(model.matrix(fit)) - length(fixef(fit)) - (length(attributes(summary(fit)$varcor)$names) + 1)
    se_residual <- sigma(fit)
    if(includeRanef)
        seGlobal <- sqrt(seFixed^2 + seRandom^2 + se_residual^2)
    else
        seGlobal <- sqrt(seFixed^2 + se_residual^2)
    t_quantile <- (quant - fitted) / seGlobal
    t_probs <- pt(t_quantile, df = rdf)
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- fitted
    tb[[probName]] <- t_probs
    tb
    
}

## TODO : fix this up
simulation_probs_mermod <- function(tb, fit, quant, probName, includeRanef, nSims = 1000) {

    if (includeRanef) {
        which = "full"
    } else {
        which = "fixed"
    }

    probs_out <- predictInterval(fit, tb, which = which, level = 1 - alpha,
                              n.sims = nSims,
                              stat = "mean",
                              include.resid.var = TRUE)

    if(is.null(tb[["pred"]])) tb <- modelr::add_predictions(tb, fit)
    tb[[probName[1]]] <- probs_out$lwr
    tb
    
}

