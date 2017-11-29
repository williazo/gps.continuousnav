#'t_test
#'
#'This is a function that will test balance of the covariates with respect to a
#'continuous treatment variable. In a non-randomized experiment, such as an observational experiment,
#'we expect that certain covariates will be highly correlated with the treatment and potentially the outcome.
#'We can use this function to gauge the level of imbalance in the covariates.
#'
#'The t-tests are run comparing one reference category of the continuous variable against the remaining categories.
#'By default the continuous variable will be
#'
#'@param dt Data.frame object in R. This data frame must have covariates and tx_var as variable names
#'@param covariates
#'@param tx Continuous treatment variable.
#'@param tx_cat Factor vector specifying the levels of the continuous variable being investigated.
#'
#'@return List of objects containing \code{mean_table}, \code{t_table}, and \code{obs}
#'
#'@export

t_test <- function(dt, covariates, tx = NULL, tx_cat = NULL){
  #### Arguments #####
  # dt          : data.frame object in R. This data frame must have covariates and tx_var as variable names
  # tx_var      : character value identifying the categorical values for the treatment variables
  # covariates  : chacracter value or vector which contains probable covariates to try and balance on/fit as part of GPS
  ####################

  #adding in this initial part so that this function can handle categorical variables as well
  if(is.null(tx)==T & is.null(tx_cat)==T){
    stop("Must specify either continuous treatment vector or categorical treatment vector", call. = F)
  }
  if(is.null(tx_cat)==T & is.null(tx)== F){
    warning("Setting the number of quantiles to three.", call. = F)
    tx_cut = quant_create(tx, n = 3)
    tx_var = tx_cut$quant
  }
  if(is.null(tx_cat) == F & is.null(tx) == T){
    tx_var = tx_cat
  }
  if(is.null(tx_cat) == F & is.null(tx) == F){
    warning("Continuous treatment value and categorical treatment value supplied. Using only the categorical value")
    tx_var = tx_cat
  }

  if(sum(covariates%in%names(dt))!=length(covariates)){
    stop("Covariate names not found in data.frame")
  }

  #checking to make sure our treatment groups are the same dimension as the data.frame
  if(length(tx_var)!=nrow(dt)){
    stop("Incorrect dimenstions. Treatment variable must be same length as data.frame.", call. = F)
  }

  if(sum(tx_var%in%names(dt)) == 0){
    dt <- data.frame(dt, tx_var)
  }
  #this will be a data.frame of just the covariates variable
  cov_data = dt[, covariates]
  #now I use the model.matrix to generate the design matrix without an intercept term
  #each categorical variable will have an associated dummy variable
  prev_na<- options('na.action')#the original na.action case
  options(na.action='na.pass') #this option allows the model matrix to be constructed with NA's allowed rather than only complete cases
  design_mat = data.frame(model.matrix( ~ . - 1, data = cov_data))
  options(na.action = prev_na$na.action) #changing this back for complete case analysis

  #by default the function model.matrix will drop all of the rows that are missing all of the covariates
  #therefore I adjust the specified data.frame by excluding all of the variables that are missing all of the covariates jointly
  #dt = dt[rowSums(is.na(dt[, covariates]))==0, ]
  #combining the design matrix with the original, excluding any of the duplicate variables
  dt = data.frame(dt, design_mat[, !names(design_mat)%in%names(dt)])
  lvl = levels(dt[, "tx_var"])
  #idenfitying the total number of levels for the categorical treatment variable

  t_tbl = matrix(ncol = length(lvl), nrow = length(names(design_mat)))
  mean_tbl = matrix(ncol = length(lvl), nrow = length(names(design_mat)))
  colnames(t_tbl) = lvl; colnames(mean_tbl) = lvl
  row.names(t_tbl) = names(design_mat); row.names(mean_tbl) = names(design_mat)
  #creating empty matrices to fill for the t-statistics and mean values

  # now we want to loop this over the number of covariates, and the number of distinct levels
  for(i in names(design_mat)){
    for(j in lvl){
      # running a t.test of each categorical value versus all other categories
      results = t.test(dt[dt[ , "tx_var"] == j, i], dt[dt[ , "tx_var"] != j, i])
      #saving the resulting test statistic
      t_tbl[i, j] = results$statistic
      #calculating the mean value for the ith covariate and jth categorical level
      mean_tbl[i, j] = round(mean(dt[dt[ , "tx_var"] == j, i], na.rm = T), 2)
    }
  }
  #saving both of these tables in the tbls object
  tbls = list(t_table = t_tbl, mean_table = mean_tbl, obs = nrow(design_mat))
  return(tbls)
}
