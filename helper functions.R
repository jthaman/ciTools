## Helper functions
make_formula <- function(fixedEffects, randomEffects, rvName = "y"){
  
  fixedPart <- paste(fixedEffects, collapse = "+")
  randomPart <- paste("+ (1|", randomEffects, ")")
  formula(paste(c(rvName, " ~ ", fixedPart, randomPart), collapse = ""))
  
}



