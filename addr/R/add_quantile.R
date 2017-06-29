## generic for add_quantile

add_quantile <- function(tb, fit, prob, quantileName = NULL, ...){
  UseMethod("add_quantile", fit)
}

