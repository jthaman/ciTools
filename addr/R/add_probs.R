## add_probs generic 

add_probs <- function(tb, fit, quant, probName = NULL, comparison, ...){
  UseMethod("add_probs", fit)
}


