#'Estimating the Dose Response Function given the Generalized Propensity Score
#'
#'This function estimates the fitted value from a linear hierarchical mixed effects model, with the
#'generalized propensiy score and treatment occuring at the highest level
#'
#'@param formula This is a formula object specifying the dose-response model
#'@param data This is the data frame used to fit the model
#'@param intervals Specify the number of intervals to include.
#'
#'@import lme4
#'
#'@return Returns a lme4 object
#'
#'@export

continuous_nav <- function(formula, data, intervals = 50){

  mm_out <- lme4::lmer(formula, data) #running the mixed effects model using lme4

  #breaking down the outcome and the RHS
  formula_breakdown <- strsplit(formula, " ~ ")

  #breaking apart the RHS by terms
  pred_and_random <- unlist(strsplit(formula_breakdown[[1]][2], "[+]"))
  #extracting the random effects in the formula
  random <- pred_and_random[grep("\\|", pred_and_random)]
  #extracting the fixed predictors in the formula
  fixed <- pred_and_random[-grep("\\|", pred_and_random)]
  fixed_RHS <- as.formula(paste("~", paste(fixed, collapse = "+")))

  design_mat <- model.matrix(fixed_RHS, data)
}
