
#' Quantile Creation Function
#'
#' A simple function that allows you to create evenly spaced quantiles for a continuous varible.
#'
#' @param dt Data.fame object in R, does not need to be specified if supplying a continuous vector
#' @param continuous_var If dt is specified, then this should be a character variable specifying which variable to split on. If dt is null then this should be a continuous numeric vector.
#' @param n The number of evenly spaced quantiles to create. Default value is tertile.
#'
#' @return List containing the vector of quantile values and the continuous cutpoints used.
#'
#' @examples
#' data("iris")
#' quant_obj <- quant_create(iris$Sepal.Length)
#' table(quant_obj$quant); quant_obj$cut
#'
#' @export
#'
#building tertiles for the log density outcome measure
#quantile creation function, n specifies the number of quantiles for the outcome variable
quant_create <- function(continuous_var, n = 3, dt = NULL){
  #### Arguments #####
  # n               : the number of quantiles
  # dt              : data.frame object in R, may not need to be specified
  # continuous_var  : if dt is specified than a character variable specifying which variable to split on.
  #                   if dt is NULL then this is a continuous numeric variable to split on
  ####################

  probs = vector(length = (n-1))
  #this for loop is used to define the probabilities for cutting the continuous value
  # eg. if n = 3 then probs = c(1/3, 2/3)
  for(i in 1:(n-1)){
    probs[i] = i / n
  }
  # run this command if dt is specified
  if(is.null(dt)==F){
    if(continuous_var%in%names(dt)==T){
      #not that we include the lowest enpoint here, and that the breaks are specified starting at zero and run up to the maximum observed value
      quant = cut(dt[, continuous_var], breaks = c(min(dt[, continuous_var]), quantile(dt[ , continuous_var], probs = probs), max(dt[, continuous_var])),
                  include.lowest = T)
      #also want to save the cut point values so that these can be used as headings later
      cut_points = c(min(dt[, continuous_var]), quantile(dt[ , continuous_var], probs = probs), max(dt[, continuous_var]))
      #saving both of these objects in the quant object
      quant = list(quant = quant, cut = cut_points)
      return(quant)
    }else{
      stop("Variable specified is not in the data.frame")
    }

  } else { #if not specified then run this command
    quant = cut(continuous_var, breaks = c(min(continuous_var), quantile(continuous_var, probs = probs), max(continuous_var)),
                include.lowest = T)
    cut_points = c(min(continuous_var), quantile(continuous_var, probs = probs), max(continuous_var))
    quant = list(quant = quant, cut = cut_points)
    return(quant)
  }
}

