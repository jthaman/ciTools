## generic for add_quantile

add_quantile <- function(tb, fit, alpha = 0.05, quantileName = NULL, ...){
  UseMethod("add_ci")
}

