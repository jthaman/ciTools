                                        #add_ci method for lme4/merMod objects
                                        #This method is a wrapper for sub-functions that will do bootstrap or parametric CIs
## TODO : include bootMer and predictInterval functions
## library(arm)
## library(lme4)
## library(merTools)

## this may be broken right now

add_ci.merMod <- function(tb, fit, 
                          alpha = 0.05, ciType = "parametric", 
                          includeRanef = T, ciNames = c("LCB", "UCB"), ...) {
    if(ciType == "bootstrap") {
        bootstrap_ci_mermod(tb, fit, alpha, ciNames, ...)
    } else if (ciType == "parametric") {
        parametric_ci_mermod(tb, fit, alpha, ciNames, includeRanef)
    } else if (ciType == "simulation") {
        simulation_ci_mermod(tb, fit, alpha, ciNames, includeRanef, ...)
    } else if (!(ciType %in% c("bootstrap", "parametric", "simulation"))) {
        stop("Incorrect type specified!")
    }
}

                                        #parametric CI for merMod objects

parametric_ci_mermod <- function(tb, fit, alpha, ciNames, includeRanef){
    
    X <- model.matrix(fit)
    vcovBetaHat <- vcov(fit)
    
    seFixed <- X %*% vcovBetaHat %*% t(X) %>% 
        diag() %>%
        sqrt()
    
    seRandom <- arm::se.ranef(fit)[[1]][1,]
    rdf <- nrow(model.matrix(fit2)) - length(fixef(fit2)) - (length(attributes(summary(fit2)$varcor)$names) + 1)
    if(includeRanef) seGlobal <- sqrt(seFixed^2 + seRandom^2) else
                                                                  seGlobal <- seFixed
    
    if(is.null(tb[["pred"]])) tb <- modelr::add_predictions(tb, fit)
    tb[[ciNames[1]]] <- tb[["pred"]] + qt(alpha/2, df = rdf) * seGlobal
    tb[[ciNames[2]]] <- tb[["pred"]] + qt(1 - alpha/2, df = rdf) * seGlobal
    tb
    
}

## simulation method using predictInterval
## Currently in Progress

simulation_ci_mermod <- function(tb, fit, alpha, ciNames, includeRanef, nSims = 1000) {

    if (includeRanef) {
        which = "full"
    } else {
        which = "fixed"
    }

    ci_out <- predictInterval(fit, tb, which = which, level = 1 - alpha,
                              n.sims = nSims,
                              stat = "mean",
                              include.resid.var = FALSE)

    if(is.null(tb[["pred"]])) tb <- modelr::add_predictions(tb, fit)
    tb[[ciNames[1]]] <- ci_out$lwr
    tb[[ciNames[2]]] <- ci_out$upr
    tb
    
}

## New Bootstrapping method with boorMer
## currently in progress
bootstrap_ci_mermod <- function(tb, fit, alpha, ciNames, includeRanef, nBS = 1000){
    

    if (includeRanef) {
        use.u = FALSE
    } else {
        use.u = TRUE
    }

    ci_out <- bootMer()

    if(is.null(tb[["pred"]])) tb <- modelr::add_predictions(tb, fit)
    tb[[ciNames[1]]] <- tb[["pred"]] + qnorm(alpha/2) * seFixed + wholePlotBS[1]
    tb[[ciNames[2]]] <- tb[["pred"]] + qnorm(1 - alpha/2) * seFixed + wholePlotBS[2]
    tb
    
}

##Bootstrap CIs for merMod objects
##Note that the bootstrapping is done on the random effects only. 

oldbootstrap_ci_mermod <- function(tb, fit, alpha, ciNames, nBS = 999){
    
    X <- model.matrix(fit)
    vcovBetaHat <- vcov(fit)
    
    seFixed <- X %*% vcovBetaHat %*% t(X) %>% 
        diag() %>%
        sqrt()
    
    wholePlotBS <- sample_ranefs(fit, nBS, alpha) 
    if(is.null(tb[["pred"]])) tb <- modelr::add_predictions(tb, fit)
    tb[[ciNames[1]]] <- tb[["pred"]] + qnorm(alpha/2) * seFixed + wholePlotBS[1]
    tb[[ciNames[2]]] <- tb[["pred"]] + qnorm(1 - alpha/2) * seFixed + wholePlotBS[2]
    tb
    
}

sample_ranefs <- function(fit, nBS, alpha){
    if(length(ranef(fit)[[1]][[1]]) < 5) warning(
                                             "Fewer than 5 whole plots; consider using parametric methods")
    ranef(fit)[[1]][[1]] %>% 
        base::sample(size = nBS, replace = T) %>%
        quantile(probs = c(alpha/2, 1-alpha/2))
}
