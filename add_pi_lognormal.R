## Add pi method for lognormal objects
add_pi_logNormal <- function(tb, fit, alpha = 0.05, piNames = c("LPB", "UPB")){
  
  out <- predict(fit, tb, interval = "prediction", level = 1 - alpha)
  if(is.null(tb[["pred"]])) tb[["pred"]] <- exp(out[, 1])
  tb[[piNames[1]]] <- exp(out[, 2])
  tb[[piNames[2]]] <- exp(out[, 3])
  tb
  
}


