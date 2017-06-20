# Add CI plus generics

# Defining the generic for "add_ci"
add_ci <- function(tb, fit, alpha = 0.05, ciNames = c("LCB", "UCB"), ...){
  UseMethod("add_ci")
}
