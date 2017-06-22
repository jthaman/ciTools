## add_probs generic 

add_probs <- function(tb, fit, alpha = 0.05, quantileName = probs, ...){
  UseMethod("add_probs", fit)
}


