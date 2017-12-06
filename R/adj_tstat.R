#' adj_tstat
#'
#' This function allows the user to calculate the t-test statistics adjusted for the Generalized Propensity Score.
#'It follows the development of Hirano and Imbens (2004).
#'
#'@param dt Data.frame object in R. This data frame must have quant and covs as variable names
#'@param quant Character value identifying the categorical values for the treatment variables
#'@param gps_mat Matrix with the generalized propensity score values for each individual evaluated at the median of each treatment interval
#               (this can be obtained from the normal_gps() function which returns the gps object)
#'@param covs Character value or vector that contains the covariate to check for balance
#'@param block_size Scalar value that defines the number of blocks to create based on the GPS values evaluated at the median of the quantile. The default is five blocks.
#'
#'@return List of objects
#'
#'@export

adj_tstat = function(dt, quant, gps_mat, covs, block_size = 5){
  #### Arguments #####
  # dt          : data.frame object in R. This data frame must have quant and covs as variable names
  # quant       : character value identifying the categorical values for the treatment variables
  # block_size  : scalar value that defines how many quantiles to use for blocking
  # gps_mat     : matrix with the generalized propensity score values for each individual evaluated at the median of each treatment interval
  #               (this can be obtained from the normal_gps() function which returns the gps object)
  # covs        : character value or vector that contains the covariate to check for balance
  ####################

  library(ggplot2)

  if(sum(quant%in%names(dat))!=length(quant)){
    stop("Categorical treatment variable is not in data.frame", call. = F)
  }

  if(sum(covs%in%names(dat))!=length(covs)){
    stop("Balancing covariates are not in data.frame", call. = F)
  }

  #building in a simple check to make sure that the input data has the right dimensions
  if(length(levels(dt[, quant])) == ncol(gps_mat)){
    #creating an empty object that will contain the full list for each of the treatment categorical groups
    full_results = NULL

    #modifying the data set so that it will ocnvert factor/character variables to dummy variables
    prev_na<- options('na.action')#the original na.action case
    options(na.action='na.pass') #this option allows the model matrix to be constructed with NA's allowed rather than only complete cases
    design_mat = data.frame(model.matrix( ~ . - 1, data = dt[, covs]))
    options(na.action = prev_na$na.action)
    dt = data.frame(dt, design_mat) #this will create some duplicate variables for the non-factor variables but this is ok

    #looping over the number of different treatment categories
    for(j in 1:length(levels(dt[, quant]))){
      #now we want to look at the range of the GPS evaluated at the median of the corresponding treatment category
      #the block size determines how many groups we want to break this up into
      #this also requires the quant_create() function that was created above
      gps_blocks = quant_create(continuous_var =  gps_mat[which(dt[, quant]==levels(dt[, quant])[j]), j], n = block_size)

      #now I want to categorize all of the GPS values from that median value for the other treatment categories
      #then this is combined into the data.frame object that was provided as an initial argument
      dt_gps = data.frame(dt, gps_quant = cut(gps_mat[, j], breaks = gps_blocks$cut, include.lowest = T), gps = gps_mat[, j])

      #creating empty matrices to fill as I loop over the covariates and the
      mean_mat = matrix(nrow = length(names(design_mat)), ncol = length(levels(dt_gps$gps_quant)))
      se_mat = matrix(nrow = length(names(design_mat)), ncol = length(levels(dt_gps$gps_quant)))
      row.names(mean_mat) = names(design_mat); colnames(mean_mat) = levels(dt_gps$gps_quant)
      row.names(se_mat) = names(design_mat); colnames(se_mat) = levels(dt_gps$gps_quant)

      #we want to alternate which treatment group is the comparison depending on which median value of the GPS is used
      #therefore we set this treatment level to correspond the the jth column of the gps values
      #so if we are using the first treatment interval, then we will be looking at the first gps column which is evaluated for all individuals at the median of the first interval
      quant_lvl = levels(dt_gps[, quant])[j]

      #An important condition to check when using the generalized propensity score is that both groups have common support
      support_check <- ggplot(data = dt_gps, aes(x = gps, y=..density.., fill = factor(dt_gps[, quant], levels = quant_lvl)), color = "black", alpha = 0.3)+
        geom_histogram()+
        xlab(paste("Generalized Propensity Score:\n Reference Tx Group", quant_lvl))+
        ylab("Density")+
        guides(fill = F)

      #now we want to drop all of the observartions that are not in the support of the reference group
      outside_support <- which((dt_gps[dt_gps[, quant]!=quant_lvl,]$gps <
                                  min(dt_gps[dt_gps[, quant]==quant_lvl,]$gps, na.rm = T))
                               | (dt_gps[dt_gps[, quant]!=quant_lvl,]$gps >
                                    max(dt_gps[dt_gps[, quant]==quant_lvl,]$gps, na.rm = T)))

      #checking to see if we have at least one value outside the support, otherwise we don't need to make any adjustments
      if(length(outside_support) >= 1){
        dt_gps <- dt_gps[-outside_support, ] #removing the rows that are outside the support
      }

      support_adj <- ggplot(data = dt_gps, aes(x = gps, y=..density.., fill = factor(dt_gps[, quant], levels = quant_lvl)), color = "black", alpha = 0.3)+
        geom_histogram()+
        xlab(paste("Generalized Propensity Score:\n Reference Tx Group", quant_lvl))+
        ylab("Density")+
        guides(fill = F)

      #now I want to use an iteration counter that will let me fill up the corresponding matrices since the loop is over a character variable not a numeric variable
      #iteration will correspond to change the covariate of interest
      iter = 1

      #for each covariate we will have a final t-test statistic adjusted for the GPS, mean difference of the differences, and SE of difference of the differences
      #t-statistic is simply the mean_diff_diff/se_diff_diff
      t_stat_adjust = vector(length = length(covs))
      mean_diff_diff = vector(length = length(covs))
      se_diff_diff = vector(length = length(covs))


      #now we will loop over each covariate
      for(i in names(design_mat)){
        #setting k equal to the number of different blocks that were defined above (these levels will change for each GPS column that is evaluated)
        k = levels(dt_gps$gps_quant)

        #now we calcualte the mean difference between the reference treatment group with the GPS block and the remaining reference groups with the same GPS block value for covariate i
        mean_diff = lapply(k, function(k) mean(dt_gps[(dt_gps[, quant]==quant_lvl) &  (dt_gps$gps_quant==k), i], na.rm = T) - mean(dt_gps[(dt_gps[, quant]!=quant_lvl) & (dt_gps$gps_quant==k), i], na.rm = T))
        #this mean_diff object will be a list of mean differences that equal the number of blocks for the GPS
        mean_vec = unlist(mean_diff)
        #therefore to get a vector of all these differences we unlist this object
        #mean_vec will be the mean difference for covariate i at each block level between the ref treatment group and the remaining treatment groups

        #now we want to know the total number of observations within each block so that we can do a weighted mean/SE for the final value
        #uses the same lapply and unlist method above
        size = lapply(k, function(k) nrow(dt_gps[(dt_gps$gps_quant==k), ]))
        size_vec = unlist(size)

        #placing the mean difference for each block on the ith row
        #each row is a covariate and each column is a GPS block size
        mean_mat[iter, ] = mean_vec

        #standard error are calcualted as se = sqrt(s1^2/n1 + s2^2/n2)
        #the two groups are distinguished by the reference treatment interval
        #I also need to build in this check in case there is less than or equal to one observation in each group since there is no variability in that case
        se = lapply(k, function(k) ifelse(length(dt_gps[(dt_gps[, quant]==quant_lvl) &  (dt_gps$gps_quant==k), i])<=1 | length(dt_gps[(dt_gps[, quant]!=quant_lvl) & (dt_gps$gps_quant==k), i])<=1,
                                          0,
                                          sqrt(var(dt_gps[(dt_gps[, quant]==quant_lvl) & (dt_gps$gps_quant==k), i], na.rm = T) / length(dt_gps[(dt_gps[, quant]==quant_lvl) &  (dt_gps$gps_quant==k), i])
                                               + var(dt_gps[(dt_gps[, quant]!=quant_lvl) & (dt_gps$gps_quant==k), i], na.rm = T) / length(dt_gps[(dt_gps[, quant]!=quant_lvl) & (dt_gps$gps_quant==k), i]))))
        se_vec = unlist(se)
        #now for each covariate we have a standard error for each block of the GPS difference

        #storing all of these standard errors into a matrix with each row being a covariate and each column being a GPS block
        se_mat[iter, ] = se_vec

        #now for each covariate we want to calculate the overal difference of means which is weighted by the number of observations within each GPS block
        mean_diff_diff[iter] = mean(rep(mean_vec, size_vec)) #overall difference of means

        #we also calculate the total standard error based on the standard error values for each GPS block
        #SE_total = sqrt(s1^2/n1+...+sk^2/nk)
        #these values are indexed by the covariate
        se_diff_diff[iter] = sqrt(sum(se_vec^2 / size_vec)) #overall standard error for difference of difference
        t_stat_adjust[iter] = mean_diff_diff[iter] / se_diff_diff[iter]
        #increasing the iteration after each loop
        iter = iter + 1
      }
      #combining all of the values for each covariate into a list
      results = list(mean_mat = mean_mat, se_mat = se_mat, size = size_vec, mean_diff = mean_diff_diff,
                     se_diff = se_diff_diff, t_adjust = t_stat_adjust,
                     sup_gg = support_check, sup_adj = support_adj,
                     outside_count = length(outside_support),
                     var_names = names(design_mat))
      #now for each treatment interval we store the list of objects within full results, which is itself a list
      full_results[[j]] = results
    }
    return(full_results)
  }

  #return an error message if the number of treatment categories is not equal to the number of GPS columns
  else{
    stop("Number of quantiles does not match number of GPS columns", call. = F)
  }
}
