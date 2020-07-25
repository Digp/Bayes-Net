
probability <- function(angle = NULL, ntinfo_disc = NULL){
  
  ntinfo_m <- ntinfo_disc
  len <- max(ntinfo_m, na.rm=TRUE)
  inter_list <- 1:len
  
  ######################Run_Code###################
  
  #Calculate the number of nucleotides of our dataset
  tot <- length(ntinfo_disc$eta)
  
  ######For eta
  
  if(angle == "eta"){
    eta <- vector("list", length = length(inter_list))
      
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"eta"] == i)
      if(length(index) > 0){
        dataframe2 <- ntinfo_m[index, ]
        n <- length(index)
        eta[[i]]$total <- n / tot
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"theta"] == j)
          if (length(indim) == 0) {
            eta[[i]]$theta[j] <- 0
          } else {
            eta[[i]]$theta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"eta_menos_1"] == j)
          if (length(indim) == 0) {
            eta[[i]]$eta_menos_1[j] <- 0
          } else {
            eta[[i]]$eta_menos_1[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"theta_menos_1"] == j)
          if (length(indim) == 0) {
            eta[[i]]$theta_menos_1[j] <- 0
          } else {
            eta[[i]]$theta_menos_1[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"eta_mas_1"] == j)
          if (length(indim) == 0) {
            eta[[i]]$eta_mas_1[j] <- 0
          } else {
            eta[[i]]$eta_mas_1[j] <- length(indim) / n
          }
        }
        
      }else{
        eta[[i]]$theta <- rep(0, len)
        eta[[i]]$eta_menos_1 <- rep(0, len)
        eta[[i]]$theta_menos_1 <- rep(0, len)
        eta[[i]]$eta_mas_1 <- rep(0, len)
        eta[[i]]$total <- 0
        
      }
    }
    return(eta)
  }
  
  
  #####For beta
    
  if(angle == "theta"){
    theta <- vector("list", length = length(inter_list))
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"theta"] == i)
      if(length(index) > 0){
        dataframe2 <- ntinfo_m[index, ]
        n <- length(index)
        theta[[i]]$total <- n / tot
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"eta"] == j)
          if (length(indim) == 0) {
            theta[[i]]$eta[j] <- 0
          } else {
            theta[[i]]$eta[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"theta_menos_1"] == j)
          if (length(indim) == 0) {
            theta[[i]]$theta_menos_1[j] <- 0
          } else {
            theta[[i]]$theta_menos_1[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"theta_mas_1"] == j)
          if (length(indim) == 0) {
            theta[[i]]$theta_mas_1[j] <- 0
          } else {
            theta[[i]]$theta_mas_1[j] <- length(indim) / n
          }
        }
        for (j in 1:length(inter_list)){
          indim <- which(dataframe2[,"eta_mas_1"] == j)
          if (length(indim) == 0) {
            theta[[i]]$eta_mas_1[j] <- 0
          } else {
            theta[[i]]$eta_mas_1[j] <- length(indim) / n
          }
        }
      }else{
        theta[[i]]$eta <- rep(0, len)
        theta[[i]]$eta_mas_1 <- rep(0, len)
        theta[[i]]$theta_menos_1 <- rep(0, len)
        theta[[i]]$theta_mas_1 <- rep(0, len)
        theta[[i]]$total <- 0
      }
    }
    return(theta)
  }
  
  
  #####For eta_menos_1
  
  if(angle == "eta_menos_1"){
    eta_menos_1 <- vector("list", length = length(inter_list))
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"eta_menos_1"] == i)
      if(length(index) > 0){
        dataframe2 <- ntinfo_m[index, ]
        n <- length(index)
        eta_menos_1[[i]]$total <- n / tot
      }else{
        eta_menos_1[[i]]$total <- 0
      }
    }
    return(theta)
  }
  
  #For theta_menos_1
  
  if(angle == "theta_menos_1"){
    theta_menos_1 <- vector("list", length = length(inter_list))
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"theta_menos_1"] == i)
      if(length(index) > 0){
        dataframe2 <- ntinfo_m[index, ]
        n <- length(index)
        theta_menos_1[[i]]$total <- n / tot
      }else{
        theta_menos_1[[i]]$total <- 0
      }
    }
    return(theta_menos_1)
  }
  
  #For eta_mas_1
  
  if(angle == "eta_mas_1"){
    eta_mas_1 <- vector("list", length = length(inter_list))
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"eta_mas_1"] == i)
      if(length(index) > 0){
        dataframe2 <- ntinfo_m[index, ]
        n <- length(index)
        eta_mas_1[[i]]$total <- n / tot
      }else{
        eta_mas_1[[i]]$total <- 0
      }
    }
    return(eta_mas_1)
  }
  
  #For theta_mas_1
  
  if(angle == "theta_mas_1"){
    theta_mas_1 <- vector("list", length = length(inter_list))
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"theta_mas_1"] == i)
      if(length(index) > 0){
        dataframe2 <- ntinfo_m[index, ]
        n <- length(index)
        theta_mas_1[[i]]$total <- n / tot
      }else{
        theta_mas_1[[i]]$total <- 0
      }
    }
    return(theta_mas_1)
  }
}
