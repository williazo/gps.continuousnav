#' Working Directory Function
#'
#' @param comp_choice Character value specifying which of the three computers I am using. Options are: macbook, iMac, and may_lab
#' @param type_save Character value specifying which type of file to save, e.g. Image
#'
#' @return Returns the working directory of interest
#'
#' @export


#creating a function that will allow me to switch between computers easier
wd_choice <- function(comp_choice, type_save=NULL){
  macbook_wd  <- "~/Documents/UCLA/"
  iMac_wd     <- "~/Documents/Propensity_Project/"
  may_lab_wd  <- "~/PropensityScoreModeling/"
  if(is.null(type_save)==T){
    #if the type_save argument is null then we assume that we are reading in raw data
    if (comp_choice=="macbook"){
      paste0(macbook_wd, "ECOS Study/")
    } else if (comp_choice=="iMac"){
      paste0(iMac_wd, "GPS_DataFiles/")
    } else if (comp_choice=="may_lab"){
      paste0(may_lab_wd,"DataFiles/")
    }
  }
  #Image Save#
  else if (type_save == "Image"){
    if (comp_choice=="macbook"){
      paste0(macbook_wd, "GSRM/Images/")
    } else if (comp_choice=="iMac"){
      paste0(iMac_wd, "Images/")
    } else if (comp_choice=="may_lab"){
      paste0(may_lab_wd,"Images/")
    }
  }

  #Data Save#
  else if (type_save == "Data") {
    if (comp_choice == "macbook"){
      paste0(macbook_wd, "ECOS Study/")
    }else if (comp_choice == "may_lab"){
      paste0(may_lab_wd,"DataFiles/")
    }
  }

  #Dose-Response Save#
  else if (type_save == "Dose-Response"){
    if (comp_choice == "macbook"){
      paste0(macbook_wd, "GSRM/Dose_Response_Output/")
    }else if (comp_choice == "iMac"){
      paste0(iMac_wd, "Dose_Response_Output/")
    }else if (comp_choice == "may_lab"){
      paste0(may_lab_wd, "DataFiles/Dose_Response_Output/")
    }
  }

}
