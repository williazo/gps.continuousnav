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
#'@export

continuous_nav <- function(formula, data, intervals = 50){

  mm_out <- lme4::lmer(formula, data = dat) #running the mixed effects model using lme4

}
