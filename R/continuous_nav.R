#'Estimating the Dose Response Function given the Generalized Propensity Score
#'
#'This function estimates the fitted value from a linear hierarchical mixed effects model, with the
#'generalized propensiy score and treatment occuring at the highest level
#'
#'@param formula This is a formula object specifying the dose-response model
#'@param data This is the data frame used to fit the model
#'@param tx This is the continuous treatment vector
#'@param intervals Specify the number of intervals to include.
#'@param covs Matrix of observed covariates
#'@param interact_vars Specifies a character subset of the variables from the matrix covs and adds in pairwise interactions between all included covariates
#'@param polynomial_vars Specifies a character subset of variables from the matrix covs to include as polynomial terms
#'@param polynomial_deg Scalar value that identifies to what power the polynomial variable should be raised
#'
#'
#'@import lme4
#'
#'@return Returns a lme4 object
#'
#'@export

continuous_nav <- function(formula, data, tx, intervals = 50, covs,
                           interact_vars, polynomial_vars, polynomial_deg){

  mm_out <- lme4::lmer(formula, data) #running the mixed effects model using lme4
  omega_hat <- summary(mm_out)$coefficients

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

  #now we want to generate the corresponding treatment and GPS values that are fixed along the observed range
  tx_values <- seq(min(tx), max(tx), length.out = intervals)
  #specifying intervals on the range of observed treatment
  r_values <- normal_gps(tx, covs, gps_val = tx_values,
                         interact_vars, polynomial_vars, polynomial_deg)
  #calculating the corresponding normal density
  r_values <- as.data.frame(r_values$gps_fit)
  #returning just the matrix of density values
  names(r_values) <- paste0("r", 1:intervals) #naming each of the score function estimates based on the treatment used

  result = list(mixed_model = mm_out, coef_est = omega_hat, fixed_tx = tx_values,
                fitted_gps = r_values)
}
