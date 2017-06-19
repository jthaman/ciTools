#add_pi method for lm
add_pi.lm <- function(tb, fit, alpha = 0.05, piNames = c("LPB", "UPB")){
  
  out <- predict(fit, tb, interval = "prediction", level = 1 - alpha)
  if(is.null(tb[["pred"]])) tb[["pred"]] <- out[, 1]
  tb[[piNames[1]]] <- out[, 2]
  tb[[piNames[2]]] <- out[, 3]
  tb
  
}

