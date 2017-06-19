#add_ci method for lm
add_ci.lm <- function(tb, fit, alpha = 0.05, ciNames = c("LCB", "UCB")){
  
  out <- predict(fit, tb, interval = "confidence", level = 1 - alpha)
  if(is.null(tb[["pred"]])) tb[["pred"]] <- out[, 1]
  tb[[ciNames[1]]] <- out[, 2]
  tb[[ciNames[2]]] <- out[, 3]
  tb
  
}

