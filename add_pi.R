# Defining the generic for "add_pi"
add_pi <- function(tb, fit, alpha = 0.05, piNames = c("LPB", "UPB"), ...){
  UseMethod("add_pi", fit)
}


