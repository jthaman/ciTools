#add_ci method for lme4/merMod objects
#This method is a wrapper for sub-functions that will do bootstrap or parametric CIs
## TODO : include bootMer and predictInterval functions

add_ci.lmerMod <- function(tb, fit, 
                           alpha = 0.05, ciType = "parametric", 
                           includeRanef = T, ciNames = c("LCB", "UCB"), ...){
  if(ciType == "bootstrap") return(
    bootstrap_ci_mermod(tb, fit, alpha, ciNames, ...)
  ) else
    if(ciType == "parametric") return(
      parametric_ci_mermod(tb, fit, alpha, ciNames, includeRanef)
    ) else
      if(!(ciType %in% c("bootstrap", "parametric"))) stop("Incorrect type specified!")
  
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


#Bootstrap CIs for merMod objects
#Note that the bootstrapping is done on the random effects only. 

bootstrap_ci_mermod <- function(tb, fit, alpha, ciNames, nBS = 999){
  
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

