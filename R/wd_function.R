#### Working Directory Function for Reading in Files from Different Computers ####

#creating a function that will allow me to switch between computers easier
wd_choice <- function(comp_choice, type_save){
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
  } else if (type_save == "Image"){
    if (comp_choice=="macbook"){
      paste0(macbook_wd, "GSRM/Images/")
    } else if (comp_choice=="iMac"){
      paste0(iMac_wd, "Images/")
    } else if (comp_choice=="may_lab"){
      paste0(may_lab_wd,"DataFiles/")
    }
  } else if (type_save == "Table"){
    ### Will come back and fill this in ###
  }
}