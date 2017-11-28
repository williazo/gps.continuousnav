#### GPS Score Functions: Continuous Treatment Domain
## Updated: 11/14/2017
## By: Justin Williams

#' Generalized Propensity Score for Continuous Treatment Domain
#'

# a simple function to evaluate the conditional density given a treatment value, design matrix, and coefficient estimates
R_score <- function(tx, X, beta, sigma){
  #tx     : is a continuous treatment vector
  #X      : design matrix used when estimating the beta and sigma parameters
  #beta   : fixed coefficient estimate from conditional distribution of tx given X
  #sigma  : standard deviation estimate from conditional distribution of tx given X
  (1/ sqrt(2 * pi * sigma^2)) * exp((-1/(2 * sigma^2)) * (tx - X %*% beta)^2)
}

# we want to write a function that will caclulate the GPS values at specific values of the treatment
normal_gps <- function(tx, covs, gps_val = NULL, interact_vars = NULL, polynomial_vars = NULL,
                       polynomial_deg = NULL, variable_selection = F){
  #### Arguments #####
  # tx                  : vector with the continuous treatment value, used for finding the initial parameter MLE's
  # covs                : matrix of observed covariates
  # gps_val             : scalar value or vector which contains the values to find the estimated generalized propensity score
  # polynomial_vars     : specifies a character subset of variables from the matrix covs to include as polynomial terms
  # polynomial_deg      : scalar value that identifies to what power the polynomial variable should be raised
  # interact_vars       : specifies a character subset of the variables from the matrix covs and adds in pairwise interactions between all included covariates
  # variable_selection  : indicator whether to perform variable selection using AIC backwards selection
  ####################
  library(plyr)
  library(dplyr)

  #standardizing the numeric values
  num_var <- sapply(covs, is.numeric)
  covs[ , num_var] <- lapply(covs[ , num_var], scale)

  #I want to include an option where you can specify interactions between the listed covariates as well as higher order polynomials
  if(is.null(polynomial_deg) == F & is.numeric(polynomial_deg)==F){
    stop("Degree of polynomial must be numeric")
  }
  #here I am buliding a stop point if the polynomial degrees is specified but polynomial variables are not
  if(is.null(polynomial_vars) == T & is.null(polynomial_deg) == F){
    stop("Must specify polynomial_vars with polynomial_deg")
  }
  #if polynomial variables are specified then
  if(is.null(polynomial_vars) == F & is.null(polynomial_deg) == F & polynomial_deg > 1 ){
    power_covs <- covs[ ,polynomial_vars] ^ polynomial_deg
    colnames(power_covs) <- paste0(polynomial_vars,"_",polynomial_deg)
    covs <- data.frame(covs, power_covs)
  }


  if(is.null(interact_vars) == F){
    #I want to restrict the interaction variables to only pairwise interactions as higher order interactions are often intractable
    product_vars = NULL
    final_names = NULL
    for (i in 1:(length(interact_vars)-1)){
      for(j in (i+1):length(interact_vars)){
        cross_prod = covs[, interact_vars[i]] * covs[, interact_vars[j]]
        name_vec = paste0(interact_vars[i], "X", interact_vars[j])
        final_names = c(final_names, name_vec)
        product_vars = cbind(product_vars, cross_prod)
      }
    }
    colnames(product_vars) = final_names
    covs <- data.frame(covs, product_vars)
  }

  #here I need to build in a stoping mechanism that check to make sure that the number of parameters is not larger than the number of observations
  if(ncol(covs)>=nrow(covs)){
    stop("Model is oversaturated. Reduce the number of interactions, polynomials, or covariates such that p<n")
  }

  #in this first part we specify the observed treatment values to get parameter estiamtes for the normal distribution

  #transforming the covariate matrix provided into a design matrix so that variables that are factors are converted into dummy numeric variables
  prev_na<- options('na.action')#the original na.action case
  options(na.action='na.pass') #this option allows the model matrix to be constructed with NA's allowed rather than only complete cases
  design_mat = data.frame(model.matrix( ~ . - 1, data = covs))
  options(na.action = prev_na$na.action)

  dt = data.frame(tx, design_mat)
  #combining the vector and matrices into a data.frame object. note that these need to be the same length
  lm_mod = lm(tx ~ . - 1, data = dt)
  #fittting the normal model using all covariates linearly
  beta = coef(lm_mod)
  #retrieving the beta coefficients
  sigma = summary(lm_mod)$sigma
  #retrieving sigma hat

  #adding in a stepwise selection based on AIC criterion
  if(variable_selection == T){
    var_select <- step(lm_mod, lm(tx ~ 1, data = dt), direction = "backward", trace = F)
    form_select <- formula(var_select)
    beta_select <- coef(var_select)
    num_select <- length(beta_select) #number of variables selected
    sigma_select <- summary(var_select)$sigma

    prev_na<- options('na.action')#the original na.action case
    options(na.action='na.pass') #this option allows the model matrix to be constructed with NA's allowed rather than only complete cases
    design_mat_select = data.frame(model.matrix(form_select, data = dt))
    options(na.action = prev_na$na.action)

    X_select = as.matrix(design_mat_select)
    X_select[is.na(X_select)] <- 0

    R_hat_select = R_score(tx = tx, X = X_select, beta = beta_select, sigma = sigma_select)
  }


  #now in the second part we get a vector of generalized propensity scores for each individual based on their covariates
  X = as.matrix(design_mat)

  #creating the design matrix without an intercept term, since the normal model above did not include an intercept
  X[is.na(X)] <- 0
  #assigning zero to all missing values

  #finding the observed GPS, i.e. R_hat for the full model
  R_hat = R_score(tx = tx, X = X, beta = beta, sigma = sigma)


  #if we have values that we want to evaluate the GPS at then we can specify them here
  if(is.null(gps_val)==F){
    # now we have two possible paths to evaluate the gps depending on if the function is provided as a scalar or a vector
    if(length(gps_val) == 1){ #scalar GPS
      gps_vec = rep(gps_val, times = nrow(X))
      #creating a vector that repeats the GPS value to match the length of our design matrix

      gps = R_score(tx = gps_vec, X = X, beta = beta, sigma = sigma)
      #using the normal probability density function and the parameter estimates from earlier to retrieve GPS estimate
    }

    else{
      #vector GPS
      #creating empty matrices to hold the given GPS values, and the calculated GPS scores for each individual
      val_mat = matrix(nrow=nrow(X), ncol = length(gps_val))
      #each column of the val_mat will be identical
      gps_mat = matrix(nrow=nrow(X), ncol = length(gps_val))
      colnames(gps_mat) = gps_val

      #now we want to loop over the number of different values that were given
      for(i in 1:ncol(gps_mat)){
        val_mat[, i] = rep(gps_val[i], times = nrow(gps_mat))
        #repeating the ith value to become a vector to evaluate in normal equation

        gps_mat[, i] = R_score(tx = val_mat[, i], X = X, beta = beta, sigma = sigma)
        #then for each alue we evaluate the generalized propensity score and save the corresponding vector into the ith column
      }
    }
  }

  if(variable_selection == T){
    #saving all of this information into one object
    est = list(beta = beta, sigma = sigma, lm = lm_mod, obs_GPS=R_hat, gps_fit = gps_mat,
               select_mod = var_select, obs_GPS_select = R_hat_select, num_select = num_select)
  }
  else{
    est = list(beta = beta, sigma = sigma, lm = lm_mod, obs_GPS=R_hat, gps_fit = gps_mat)
  }
  return(est)
}
