
probability <- function(angle = NULL, ntinfo_disc = NULL, ntinfo = NULL){
  
  ntinfo_m <- ntinfo_disc
  len <- max(ntinfo_m, na.rm=TRUE)
  inter_list <- 1:len
  
  ######################Run_Code###################
  
  #Calculate the number of nucleotides of our dataset
  tot <- length(ntinfo_disc$alpha)
  
  ######For alpha
  
  if(angle == "alpha"){
    alpha <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"alpha"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for beta analysis
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        alpha[[i]]$beta$prop <- gaussians_beta$prop
        alpha[[i]]$beta$mean <- gaussians_beta$mean
        alpha[[i]]$beta$sd <- gaussians_beta$sd
        alpha[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for gamma
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        alpha[[i]]$gamma$prop <- gaussians_gamma$prop
        alpha[[i]]$gamma$mean <- gaussians_gamma$mean
        alpha[[i]]$gamma$sd <- gaussians_gamma$sd
        alpha[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for delta
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        alpha[[i]]$delta$prop <- gaussians_delta$prop
        alpha[[i]]$delta$mean <- gaussians_delta$mean
        alpha[[i]]$delta$sd <- gaussians_delta$sd
        alpha[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for epsilon
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        alpha[[i]]$epsilon$prop <- gaussians_epsilon$prop
        alpha[[i]]$epsilon$mean <- gaussians_epsilon$mean
        alpha[[i]]$epsilon$sd <- gaussians_epsilon$sd
        alpha[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for zeta
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        alpha[[i]]$zeta$prop <- gaussians_zeta$prop
        alpha[[i]]$zeta$mean <- gaussians_zeta$mean
        alpha[[i]]$zeta$sd <- gaussians_zeta$sd
        alpha[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        alpha[[i]]$chi$prop <- gaussians_chi$prop
        alpha[[i]]$chi$mean <- gaussians_chi$mean
        alpha[[i]]$chi$sd <- gaussians_chi$sd
        alpha[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        alpha[[i]]$total <- n / tot
        
        
      }else{
      
        
        alpha[[i]]$beta$prop <- 0
        alpha[[i]]$beta$mean <- 0
        alpha[[i]]$beta$sd <- 0
        alpha[[i]]$beta$max_dens <- 0
        
        alpha[[i]]$gamma$prop <- 0
        alpha[[i]]$gamma$mean <- 0
        alpha[[i]]$gamma$sd <- 0
        alpha[[i]]$gamma$max_dens <- 0
        
        alpha[[i]]$delta$prop <- 0
        alpha[[i]]$delta$mean <- 0
        alpha[[i]]$delta$sd <- 0
        alpha[[i]]$delta$max_dens <- 0
        
        alpha[[i]]$epsilon$prop <- 0
        alpha[[i]]$epsilon$mean <- 0
        alpha[[i]]$epsilon$sd <- 0
        alpha[[i]]$epsilon$max_dens <- 0
        
        alpha[[i]]$zeta$prop <- 0
        alpha[[i]]$zeta$mean <- 0
        alpha[[i]]$zeta$sd <- 0
        alpha[[i]]$zeta$max_dens <- 0
        
        alpha[[i]]$chi$prop <- 0
        alpha[[i]]$chi$mean <- 0
        alpha[[i]]$chi$sd <- 0
        alpha[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        alpha[[i]]$total <- n / tot
        
        
      }
    }
    
    return(alpha)
  }
  
  
  if(angle == "beta"){
    beta <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"beta"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for beta analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        beta[[i]]$alpha$prop <- gaussians_alpha$prop
        beta[[i]]$alpha$mean <- gaussians_alpha$mean
        beta[[i]]$alpha$sd <- gaussians_alpha$sd
        beta[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for gamma
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        beta[[i]]$gamma$prop <- gaussians_gamma$prop
        beta[[i]]$gamma$mean <- gaussians_gamma$mean
        beta[[i]]$gamma$sd <- gaussians_gamma$sd
        beta[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for delta
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        beta[[i]]$delta$prop <- gaussians_delta$prop
        beta[[i]]$delta$mean <- gaussians_delta$mean
        beta[[i]]$delta$sd <- gaussians_delta$sd
        beta[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for epsilon
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        beta[[i]]$epsilon$prop <- gaussians_epsilon$prop
        beta[[i]]$epsilon$mean <- gaussians_epsilon$mean
        beta[[i]]$epsilon$sd <- gaussians_epsilon$sd
        beta[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for zeta
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        beta[[i]]$zeta$prop <- gaussians_zeta$prop
        beta[[i]]$zeta$mean <- gaussians_zeta$mean
        beta[[i]]$zeta$sd <- gaussians_zeta$sd
        beta[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        beta[[i]]$chi$prop <- gaussians_chi$prop
        beta[[i]]$chi$mean <- gaussians_chi$mean
        beta[[i]]$chi$sd <- gaussians_chi$sd
        beta[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        beta[[i]]$total <- n / tot
        
        
      }else{
        
        beta[[i]]$alpha$prop <- 0
        beta[[i]]$alpha$mean <- 0
        beta[[i]]$alpha$sd <- 0
        beta[[i]]$alpha$max_dens <- 0
        
        beta[[i]]$gamma$prop <- 0
        beta[[i]]$gamma$mean <- 0
        beta[[i]]$gamma$sd <- 0
        beta[[i]]$gamma$max_dens <- 0
        
        beta[[i]]$delta$prop <- 0
        beta[[i]]$delta$mean <- 0
        beta[[i]]$delta$sd <- 0
        beta[[i]]$delta$max_dens <- 0
        
        beta[[i]]$epsilon$prop <- 0
        beta[[i]]$epsilon$mean <- 0
        beta[[i]]$epsilon$sd <- 0
        beta[[i]]$epsilon$max_dens <- 0
        
        beta[[i]]$zeta$prop <- 0
        beta[[i]]$zeta$mean <- 0
        beta[[i]]$zeta$sd <- 0
        beta[[i]]$zeta$max_dens <- 0
        
        beta[[i]]$chi$prop <- 0
        beta[[i]]$chi$mean <- 0
        beta[[i]]$chi$sd <- 0
        beta[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        beta[[i]]$total <- n / tot
        
        
        
      }
    }
    
    return(beta)
  }
  

  #Check for gamma
  if(angle == "gamma"){
    gamma <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"gamma"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for gamma analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        gamma[[i]]$alpha$prop <- gaussians_alpha$prop
        gamma[[i]]$alpha$mean <- gaussians_alpha$mean
        gamma[[i]]$alpha$sd <- gaussians_alpha$sd
        gamma[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for gamma
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        gamma[[i]]$beta$prop <- gaussians_beta$prop
        gamma[[i]]$beta$mean <- gaussians_beta$mean
        gamma[[i]]$beta$sd <- gaussians_beta$sd
        gamma[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for delta
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        gamma[[i]]$delta$prop <- gaussians_delta$prop
        gamma[[i]]$delta$mean <- gaussians_delta$mean
        gamma[[i]]$delta$sd <- gaussians_delta$sd
        gamma[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for epsilon
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        gamma[[i]]$epsilon$prop <- gaussians_epsilon$prop
        gamma[[i]]$epsilon$mean <- gaussians_epsilon$mean
        gamma[[i]]$epsilon$sd <- gaussians_epsilon$sd
        gamma[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for zeta
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        gamma[[i]]$zeta$prop <- gaussians_zeta$prop
        gamma[[i]]$zeta$mean <- gaussians_zeta$mean
        gamma[[i]]$zeta$sd <- gaussians_zeta$sd
        gamma[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        gamma[[i]]$chi$prop <- gaussians_chi$prop
        gamma[[i]]$chi$mean <- gaussians_chi$mean
        gamma[[i]]$chi$sd <- gaussians_chi$sd
        gamma[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        gamma[[i]]$total <- n / tot
        
        
      }else{
        
        gamma[[i]]$alpha$prop <- 0
        gamma[[i]]$alpha$mean <- 0
        gamma[[i]]$alpha$sd <- 0
        gamma[[i]]$alpha$max_dens <- 0
        
        gamma[[i]]$beta$prop <- 0
        gamma[[i]]$beta$mean <- 0
        gamma[[i]]$beta$sd <- 0
        gamma[[i]]$beta$max_dens <- 0
        
        gamma[[i]]$delta$prop <- 0
        gamma[[i]]$delta$mean <- 0
        gamma[[i]]$delta$sd <- 0
        gamma[[i]]$delta$max_dens <- 0
        
        gamma[[i]]$epsilon$prop <- 0
        gamma[[i]]$epsilon$mean <- 0
        gamma[[i]]$epsilon$sd <- 0
        gamma[[i]]$epsilon$max_dens <- 0
        
        gamma[[i]]$zeta$prop <- 0
        gamma[[i]]$zeta$mean <- 0
        gamma[[i]]$zeta$sd <- 0
        gamma[[i]]$zeta$max_dens <- 0
        
        gamma[[i]]$chi$prop <- 0
        gamma[[i]]$chi$mean <- 0
        gamma[[i]]$chi$sd <- 0
        gamma[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        gamma[[i]]$total <- n / tot
        
        
        
      }
    }
    
    return(gamma)
  }
  
  #####For delta
  if(angle == "delta"){
    delta <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"delta"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for delta analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        delta[[i]]$alpha$prop <- gaussians_alpha$prop
        delta[[i]]$alpha$mean <- gaussians_alpha$mean
        delta[[i]]$alpha$sd <- gaussians_alpha$sd
        delta[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for beta
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        delta[[i]]$beta$prop <- gaussians_beta$prop
        delta[[i]]$beta$mean <- gaussians_beta$mean
        delta[[i]]$beta$sd <- gaussians_beta$sd
        delta[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for gamma
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        delta[[i]]$gamma$prop <- gaussians_gamma$prop
        delta[[i]]$gamma$mean <- gaussians_gamma$mean
        delta[[i]]$gamma$sd <- gaussians_gamma$sd
        delta[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for epsilon
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        delta[[i]]$epsilon$prop <- gaussians_epsilon$prop
        delta[[i]]$epsilon$mean <- gaussians_epsilon$mean
        delta[[i]]$epsilon$sd <- gaussians_epsilon$sd
        delta[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for zeta
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        delta[[i]]$zeta$prop <- gaussians_zeta$prop
        delta[[i]]$zeta$mean <- gaussians_zeta$mean
        delta[[i]]$zeta$sd <- gaussians_zeta$sd
        delta[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        delta[[i]]$chi$prop <- gaussians_chi$prop
        delta[[i]]$chi$mean <- gaussians_chi$mean
        delta[[i]]$chi$sd <- gaussians_chi$sd
        delta[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        delta[[i]]$total <- n / tot
        
        
      }else{
        
        delta[[i]]$alpha$prop <- 0
        delta[[i]]$alpha$mean <- 0
        delta[[i]]$alpha$sd <- 0
        delta[[i]]$alpha$max_dens <- 0
        
        delta[[i]]$beta$prop <- 0
        delta[[i]]$beta$mean <- 0
        delta[[i]]$beta$sd <- 0
        delta[[i]]$beta$max_dens <- 0
        
        delta[[i]]$gamma$prop <- 0
        delta[[i]]$gamma$mean <- 0
        delta[[i]]$gamma$sd <- 0
        delta[[i]]$gamma$max_dens <- 0
        
        delta[[i]]$epsilon$prop <- 0
        delta[[i]]$epsilon$mean <- 0
        delta[[i]]$epsilon$sd <- 0
        delta[[i]]$epsilon$max_dens <- 0
        
        delta[[i]]$zeta$prop <- 0
        delta[[i]]$zeta$mean <- 0
        delta[[i]]$zeta$sd <- 0
        delta[[i]]$zeta$max_dens <- 0
        
        delta[[i]]$chi$prop <- 0
        delta[[i]]$chi$mean <- 0
        delta[[i]]$chi$sd <- 0
        delta[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        delta[[i]]$total <- n / tot
        
        
      }
    }
    
    return(delta)
  }
  
  ##### For epsilon
  if(angle == "epsilon"){
    epsilon <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"epsilon"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for epsilon analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        epsilon[[i]]$alpha$prop <- gaussians_alpha$prop
        epsilon[[i]]$alpha$mean <- gaussians_alpha$mean
        epsilon[[i]]$alpha$sd <- gaussians_alpha$sd
        epsilon[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for epsilon
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        epsilon[[i]]$beta$prop <- gaussians_beta$prop
        epsilon[[i]]$beta$mean <- gaussians_beta$mean
        epsilon[[i]]$beta$sd <- gaussians_beta$sd
        epsilon[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for epsilon
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        epsilon[[i]]$gamma$prop <- gaussians_gamma$prop
        epsilon[[i]]$gamma$mean <- gaussians_gamma$mean
        epsilon[[i]]$gamma$sd <- gaussians_gamma$sd
        epsilon[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for epsilon
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        epsilon[[i]]$delta$prop <- gaussians_delta$prop
        epsilon[[i]]$delta$mean <- gaussians_delta$mean
        epsilon[[i]]$delta$sd <- gaussians_delta$sd
        epsilon[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for zeta
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        epsilon[[i]]$zeta$prop <- gaussians_zeta$prop
        epsilon[[i]]$zeta$mean <- gaussians_zeta$mean
        epsilon[[i]]$zeta$sd <- gaussians_zeta$sd
        epsilon[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        epsilon[[i]]$chi$prop <- gaussians_chi$prop
        epsilon[[i]]$chi$mean <- gaussians_chi$mean
        epsilon[[i]]$chi$sd <- gaussians_chi$sd
        epsilon[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        epsilon[[i]]$total <- n / tot
        
        
      }else{
        
        epsilon[[i]]$alpha$prop <- 0
        epsilon[[i]]$alpha$mean <- 0
        epsilon[[i]]$alpha$sd <- 0
        epsilon[[i]]$alpha$max_dens <- 0
        
        epsilon[[i]]$beta$prop <- 0
        epsilon[[i]]$beta$mean <- 0
        epsilon[[i]]$beta$sd <- 0
        epsilon[[i]]$beta$max_dens <- 0
        
        epsilon[[i]]$gamma$prop <- 0
        epsilon[[i]]$gamma$mean <- 0
        epsilon[[i]]$gamma$sd <- 0
        epsilon[[i]]$gamma$max_dens <- 0
        
        epsilon[[i]]$delta$prop <- 0
        epsilon[[i]]$delta$mean <- 0
        epsilon[[i]]$delta$sd <- 0
        epsilon[[i]]$delta$max_dens <- 0
        
        epsilon[[i]]$zeta$prop <- 0
        epsilon[[i]]$zeta$mean <- 0
        epsilon[[i]]$zeta$sd <- 0
        epsilon[[i]]$zeta$max_dens <- 0
        
        epsilon[[i]]$chi$prop <- 0
        epsilon[[i]]$chi$mean <- 0
        epsilon[[i]]$chi$sd <- 0
        epsilon[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        epsilon[[i]]$total <- n / tot
        
        
      }
    }
    
    return(epsilon)
  }
  
  #####For zeta
  if(angle == "zeta"){
    zeta <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"zeta"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for zeta analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        zeta[[i]]$alpha$prop <- gaussians_alpha$prop
        zeta[[i]]$alpha$mean <- gaussians_alpha$mean
        zeta[[i]]$alpha$sd <- gaussians_alpha$sd
        zeta[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for zeta
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        zeta[[i]]$beta$prop <- gaussians_beta$prop
        zeta[[i]]$beta$mean <- gaussians_beta$mean
        zeta[[i]]$beta$sd <- gaussians_beta$sd
        zeta[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for zeta
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        zeta[[i]]$gamma$prop <- gaussians_gamma$prop
        zeta[[i]]$gamma$mean <- gaussians_gamma$mean
        zeta[[i]]$gamma$sd <- gaussians_gamma$sd
        zeta[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for zeta
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        zeta[[i]]$delta$prop <- gaussians_delta$prop
        zeta[[i]]$delta$mean <- gaussians_delta$mean
        zeta[[i]]$delta$sd <- gaussians_delta$sd
        zeta[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for zeta
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        zeta[[i]]$epsilon$prop <- gaussians_epsilon$prop
        zeta[[i]]$epsilon$mean <- gaussians_epsilon$mean
        zeta[[i]]$epsilon$sd <- gaussians_epsilon$sd
        zeta[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        zeta[[i]]$chi$prop <- gaussians_chi$prop
        zeta[[i]]$chi$mean <- gaussians_chi$mean
        zeta[[i]]$chi$sd <- gaussians_chi$sd
        zeta[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        zeta[[i]]$total <- n / tot
        
        
      }else{
        
        zeta[[i]]$alpha$prop <- 0
        zeta[[i]]$alpha$mean <- 0
        zeta[[i]]$alpha$sd <- 0
        zeta[[i]]$alpha$max_dens <- 0
        
        zeta[[i]]$beta$prop <- 0
        zeta[[i]]$beta$mean <- 0
        zeta[[i]]$beta$sd <- 0
        zeta[[i]]$beta$max_dens <- 0
        
        zeta[[i]]$gamma$prop <- 0
        zeta[[i]]$gamma$mean <- 0
        zeta[[i]]$gamma$sd <- 0
        zeta[[i]]$gamma$max_dens <- 0
        
        zeta[[i]]$delta$prop <- 0
        zeta[[i]]$delta$mean <- 0
        zeta[[i]]$delta$sd <- 0
        zeta[[i]]$delta$max_dens <- 0
        
        zeta[[i]]$epsilon$prop <- 0
        zeta[[i]]$epsilon$mean <- 0
        zeta[[i]]$epsilon$sd <- 0
        zeta[[i]]$epsilon$max_dens <- 0
        
        zeta[[i]]$chi$prop <- 0
        zeta[[i]]$chi$mean <- 0
        zeta[[i]]$chi$sd <- 0
        zeta[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        zeta[[i]]$total <- n / tot
        
        
      }
    }
    
    return(zeta)
  }
  
  ####For chi
  if(angle == "chi"){
    chi <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"chi"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for chi analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        chi[[i]]$alpha$prop <- gaussians_alpha$prop
        chi[[i]]$alpha$mean <- gaussians_alpha$mean
        chi[[i]]$alpha$sd <- gaussians_alpha$sd
        chi[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for chi
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        chi[[i]]$beta$prop <- gaussians_beta$prop
        chi[[i]]$beta$mean <- gaussians_beta$mean
        chi[[i]]$beta$sd <- gaussians_beta$sd
        chi[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for chi
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        chi[[i]]$gamma$prop <- gaussians_gamma$prop
        chi[[i]]$gamma$mean <- gaussians_gamma$mean
        chi[[i]]$gamma$sd <- gaussians_gamma$sd
        chi[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for chi
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        chi[[i]]$delta$prop <- gaussians_delta$prop
        chi[[i]]$delta$mean <- gaussians_delta$mean
        chi[[i]]$delta$sd <- gaussians_delta$sd
        chi[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for chi
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        chi[[i]]$epsilon$prop <- gaussians_epsilon$prop
        chi[[i]]$epsilon$mean <- gaussians_epsilon$mean
        chi[[i]]$epsilon$sd <- gaussians_epsilon$sd
        chi[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for chi
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        chi[[i]]$zeta$prop <- gaussians_zeta$prop
        chi[[i]]$zeta$mean <- gaussians_zeta$mean
        chi[[i]]$zeta$sd <- gaussians_zeta$sd
        chi[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        n <- length(index)
        chi[[i]]$total <- n / tot
        
        
      }else{
        
        chi[[i]]$alpha$prop <- 0
        chi[[i]]$alpha$mean <- 0
        chi[[i]]$alpha$sd <- 0
        chi[[i]]$alpha$max_dens <- 0
        
        chi[[i]]$beta$prop <- 0
        chi[[i]]$beta$mean <- 0
        chi[[i]]$beta$sd <- 0
        chi[[i]]$beta$max_dens <- 0
        
        chi[[i]]$gamma$prop <- 0
        chi[[i]]$gamma$mean <- 0
        chi[[i]]$gamma$sd <- 0
        chi[[i]]$gamma$max_dens <- 0
        
        chi[[i]]$delta$prop <- 0
        chi[[i]]$delta$mean <- 0
        chi[[i]]$delta$sd <- 0
        chi[[i]]$delta$max_dens <- 0
        
        chi[[i]]$epsilon$prop <- 0
        chi[[i]]$epsilon$mean <- 0
        chi[[i]]$epsilon$sd <- 0
        chi[[i]]$epsilon$max_dens <- 0
        
        chi[[i]]$zeta$prop <- 0
        chi[[i]]$zeta$mean <- 0
        chi[[i]]$zeta$sd <- 0
        chi[[i]]$zeta$max_dens <- 0
        
        n <- length(index)
        chi[[i]]$total <- n / tot
        
        
      }
    }
    
    return(chi)
  }
}



