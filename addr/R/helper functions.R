## Helper functions
make_formula <- function(fixedEffects, randomEffects, rvName = "y"){
  
  fixedPart <- paste(fixedEffects, collapse = "+")
  randomPart <- paste("+ (1|", randomEffects, ")")
  formula(paste(c(rvName, " ~ ", fixedPart, randomPart), collapse = ""))
  
}


add_predictions2 <- function (data, model, var = "pred", ...) 
{
  data[[var]] <- stats::predict(model, data, ...)
  data
}


get_prediction_se_mermod <- function(tb, fit){
  X <- model.matrix(reformulate(attributes(terms(fit))$term.labels), tb)
  vcovBetaHat <- vcov(fit)
  X %*% vcovBetaHat %*% t(X) %>% 
    diag() %>%
    sqrt()
  
}

get_resid_df_mermod <- function(fit){
  nrow(model.matrix(fit2)) - length(fixef(fit2)) - 
    (length(attributes(summary(fit2)$varcor)$names) + 1)
}


get_residual_se <- function(fit){
  fit %>%
    VarCorr %>%
    as.data.frame %>%
    last %>%
    last
}


get_pi_mermod_var <- function(tb, fit, includeRanef){
  seFixed <- get_prediction_se_mermod(tb, fit)
  seG <- arm::se.ranef(fit)[[1]][1,]
  sigmaG <- as.data.frame(VarCorr(fit2))$sdcor[1]
  sigma <- get_residual_se(fit)
  
  if(includeRanef){return(sqrt(seFixed^2 + seG^2 + seResidual^2))} else{
    return(sqrt(seFixed^2 + sigmaG^2 + sigma^2))}
  
  
}

