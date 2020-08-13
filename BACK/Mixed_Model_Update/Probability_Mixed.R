
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
  
  if(angle == "alpha_plus"){
    alpha_plus <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"alpha_plus"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for beta analysis
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        alpha_plus[[i]]$beta$prop <- gaussians_beta$prop
        alpha_plus[[i]]$beta$mean <- gaussians_beta$mean
        alpha_plus[[i]]$beta$sd <- gaussians_beta$sd
        alpha_plus[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for gamma
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        alpha_plus[[i]]$gamma$prop <- gaussians_gamma$prop
        alpha_plus[[i]]$gamma$mean <- gaussians_gamma$mean
        alpha_plus[[i]]$gamma$sd <- gaussians_gamma$sd
        alpha_plus[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for delta
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        alpha_plus[[i]]$delta$prop <- gaussians_delta$prop
        alpha_plus[[i]]$delta$mean <- gaussians_delta$mean
        alpha_plus[[i]]$delta$sd <- gaussians_delta$sd
        alpha_plus[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for epsilon
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        alpha_plus[[i]]$epsilon$prop <- gaussians_epsilon$prop
        alpha_plus[[i]]$epsilon$mean <- gaussians_epsilon$mean
        alpha_plus[[i]]$epsilon$sd <- gaussians_epsilon$sd
        alpha_plus[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for zeta
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        alpha_plus[[i]]$zeta$prop <- gaussians_zeta$prop
        alpha_plus[[i]]$zeta$mean <- gaussians_zeta$mean
        alpha_plus[[i]]$zeta$sd <- gaussians_zeta$sd
        alpha_plus[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        alpha_plus[[i]]$chi$prop <- gaussians_chi$prop
        alpha_plus[[i]]$chi$mean <- gaussians_chi$mean
        alpha_plus[[i]]$chi$sd <- gaussians_chi$sd
        alpha_plus[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        alpha_plus[[i]]$total <- n / tot
        
        
      }else{
        
        
        alpha_plus[[i]]$beta$prop <- 0
        alpha_plus[[i]]$beta$mean <- 0
        alpha_plus[[i]]$beta$sd <- 0
        alpha_plus[[i]]$beta$max_dens <- 0
        
        alpha_plus[[i]]$gamma$prop <- 0
        alpha_plus[[i]]$gamma$mean <- 0
        alpha_plus[[i]]$gamma$sd <- 0
        alpha_plus[[i]]$gamma$max_dens <- 0
        
        alpha_plus[[i]]$delta$prop <- 0
        alpha_plus[[i]]$delta$mean <- 0
        alpha_plus[[i]]$delta$sd <- 0
        alpha_plus[[i]]$delta$max_dens <- 0
        
        alpha_plus[[i]]$epsilon$prop <- 0
        alpha_plus[[i]]$epsilon$mean <- 0
        alpha_plus[[i]]$epsilon$sd <- 0
        alpha_plus[[i]]$epsilon$max_dens <- 0
        
        alpha_plus[[i]]$zeta$prop <- 0
        alpha_plus[[i]]$zeta$mean <- 0
        alpha_plus[[i]]$zeta$sd <- 0
        alpha_plus[[i]]$zeta$max_dens <- 0
        
        alpha_plus[[i]]$chi$prop <- 0
        alpha_plus[[i]]$chi$mean <- 0
        alpha_plus[[i]]$chi$sd <- 0
        alpha_plus[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        alpha_plus[[i]]$total <- n / tot
        
        
      }
    }
    
    return(alpha_plus)
  }
  
  if(angle == "beta_plus"){
    beta_plus <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"beta_plus"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for beta_plus analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        beta_plus[[i]]$alpha$prop <- gaussians_alpha$prop
        beta_plus[[i]]$alpha$mean <- gaussians_alpha$mean
        beta_plus[[i]]$alpha$sd <- gaussians_alpha$sd
        beta_plus[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for gamma
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        beta_plus[[i]]$gamma$prop <- gaussians_gamma$prop
        beta_plus[[i]]$gamma$mean <- gaussians_gamma$mean
        beta_plus[[i]]$gamma$sd <- gaussians_gamma$sd
        beta_plus[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for delta
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        beta_plus[[i]]$delta$prop <- gaussians_delta$prop
        beta_plus[[i]]$delta$mean <- gaussians_delta$mean
        beta_plus[[i]]$delta$sd <- gaussians_delta$sd
        beta_plus[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for epsilon
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        beta_plus[[i]]$epsilon$prop <- gaussians_epsilon$prop
        beta_plus[[i]]$epsilon$mean <- gaussians_epsilon$mean
        beta_plus[[i]]$epsilon$sd <- gaussians_epsilon$sd
        beta_plus[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for zeta
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        beta_plus[[i]]$zeta$prop <- gaussians_zeta$prop
        beta_plus[[i]]$zeta$mean <- gaussians_zeta$mean
        beta_plus[[i]]$zeta$sd <- gaussians_zeta$sd
        beta_plus[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        beta_plus[[i]]$chi$prop <- gaussians_chi$prop
        beta_plus[[i]]$chi$mean <- gaussians_chi$mean
        beta_plus[[i]]$chi$sd <- gaussians_chi$sd
        beta_plus[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        beta_plus[[i]]$total <- n / tot
        
        
      }else{
        
        beta_plus[[i]]$alpha$prop <- 0
        beta_plus[[i]]$alpha$mean <- 0
        beta_plus[[i]]$alpha$sd <- 0
        beta_plus[[i]]$alpha$max_dens <- 0
        
        beta_plus[[i]]$gamma$prop <- 0
        beta_plus[[i]]$gamma$mean <- 0
        beta_plus[[i]]$gamma$sd <- 0
        beta_plus[[i]]$gamma$max_dens <- 0
        
        beta_plus[[i]]$delta$prop <- 0
        beta_plus[[i]]$delta$mean <- 0
        beta_plus[[i]]$delta$sd <- 0
        beta_plus[[i]]$delta$max_dens <- 0
        
        beta_plus[[i]]$epsilon$prop <- 0
        beta_plus[[i]]$epsilon$mean <- 0
        beta_plus[[i]]$epsilon$sd <- 0
        beta_plus[[i]]$epsilon$max_dens <- 0
        
        beta_plus[[i]]$zeta$prop <- 0
        beta_plus[[i]]$zeta$mean <- 0
        beta_plus[[i]]$zeta$sd <- 0
        beta_plus[[i]]$zeta$max_dens <- 0
        
        beta_plus[[i]]$chi$prop <- 0
        beta_plus[[i]]$chi$mean <- 0
        beta_plus[[i]]$chi$sd <- 0
        beta_plus[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        beta_plus[[i]]$total <- n / tot
        
        
        
      }
    }
    
    return(beta_plus)
  }
  
  
  #Check for gamma_plus
  if(angle == "gamma_plus"){
    gamma_plus <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"gamma_plus"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for gamma_plus analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        gamma_plus[[i]]$alpha$prop <- gaussians_alpha$prop
        gamma_plus[[i]]$alpha$mean <- gaussians_alpha$mean
        gamma_plus[[i]]$alpha$sd <- gaussians_alpha$sd
        gamma_plus[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for gamma_plus
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        gamma_plus[[i]]$beta$prop <- gaussians_beta$prop
        gamma_plus[[i]]$beta$mean <- gaussians_beta$mean
        gamma_plus[[i]]$beta$sd <- gaussians_beta$sd
        gamma_plus[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for delta
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        gamma_plus[[i]]$delta$prop <- gaussians_delta$prop
        gamma_plus[[i]]$delta$mean <- gaussians_delta$mean
        gamma_plus[[i]]$delta$sd <- gaussians_delta$sd
        gamma_plus[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for epsilon
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        gamma_plus[[i]]$epsilon$prop <- gaussians_epsilon$prop
        gamma_plus[[i]]$epsilon$mean <- gaussians_epsilon$mean
        gamma_plus[[i]]$epsilon$sd <- gaussians_epsilon$sd
        gamma_plus[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for zeta
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        gamma_plus[[i]]$zeta$prop <- gaussians_zeta$prop
        gamma_plus[[i]]$zeta$mean <- gaussians_zeta$mean
        gamma_plus[[i]]$zeta$sd <- gaussians_zeta$sd
        gamma_plus[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        gamma_plus[[i]]$chi$prop <- gaussians_chi$prop
        gamma_plus[[i]]$chi$mean <- gaussians_chi$mean
        gamma_plus[[i]]$chi$sd <- gaussians_chi$sd
        gamma_plus[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        gamma_plus[[i]]$total <- n / tot
        
        
      }else{
        
        gamma_plus[[i]]$alpha$prop <- 0
        gamma_plus[[i]]$alpha$mean <- 0
        gamma_plus[[i]]$alpha$sd <- 0
        gamma_plus[[i]]$alpha$max_dens <- 0
        
        gamma_plus[[i]]$beta$prop <- 0
        gamma_plus[[i]]$beta$mean <- 0
        gamma_plus[[i]]$beta$sd <- 0
        gamma_plus[[i]]$beta$max_dens <- 0
        
        gamma_plus[[i]]$delta$prop <- 0
        gamma_plus[[i]]$delta$mean <- 0
        gamma_plus[[i]]$delta$sd <- 0
        gamma_plus[[i]]$delta$max_dens <- 0
        
        gamma_plus[[i]]$epsilon$prop <- 0
        gamma_plus[[i]]$epsilon$mean <- 0
        gamma_plus[[i]]$epsilon$sd <- 0
        gamma_plus[[i]]$epsilon$max_dens <- 0
        
        gamma_plus[[i]]$zeta$prop <- 0
        gamma_plus[[i]]$zeta$mean <- 0
        gamma_plus[[i]]$zeta$sd <- 0
        gamma_plus[[i]]$zeta$max_dens <- 0
        
        gamma_plus[[i]]$chi$prop <- 0
        gamma_plus[[i]]$chi$mean <- 0
        gamma_plus[[i]]$chi$sd <- 0
        gamma_plus[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        gamma_plus[[i]]$total <- n / tot
        
        
        
      }
    }
    
    return(gamma_plus)
  }
  
  #####For delta_plus
  if(angle == "delta_plus"){
    delta_plus <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"delta_plus"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for delta_plus analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        delta_plus[[i]]$alpha$prop <- gaussians_alpha$prop
        delta_plus[[i]]$alpha$mean <- gaussians_alpha$mean
        delta_plus[[i]]$alpha$sd <- gaussians_alpha$sd
        delta_plus[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for beta
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        delta_plus[[i]]$beta$prop <- gaussians_beta$prop
        delta_plus[[i]]$beta$mean <- gaussians_beta$mean
        delta_plus[[i]]$beta$sd <- gaussians_beta$sd
        delta_plus[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for gamma
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        delta_plus[[i]]$gamma$prop <- gaussians_gamma$prop
        delta_plus[[i]]$gamma$mean <- gaussians_gamma$mean
        delta_plus[[i]]$gamma$sd <- gaussians_gamma$sd
        delta_plus[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for epsilon
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        delta_plus[[i]]$epsilon$prop <- gaussians_epsilon$prop
        delta_plus[[i]]$epsilon$mean <- gaussians_epsilon$mean
        delta_plus[[i]]$epsilon$sd <- gaussians_epsilon$sd
        delta_plus[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for zeta
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        delta_plus[[i]]$zeta$prop <- gaussians_zeta$prop
        delta_plus[[i]]$zeta$mean <- gaussians_zeta$mean
        delta_plus[[i]]$zeta$sd <- gaussians_zeta$sd
        delta_plus[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        delta_plus[[i]]$chi$prop <- gaussians_chi$prop
        delta_plus[[i]]$chi$mean <- gaussians_chi$mean
        delta_plus[[i]]$chi$sd <- gaussians_chi$sd
        delta_plus[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        delta_plus[[i]]$total <- n / tot
        
        
      }else{
        
        delta_plus[[i]]$alpha$prop <- 0
        delta_plus[[i]]$alpha$mean <- 0
        delta_plus[[i]]$alpha$sd <- 0
        delta_plus[[i]]$alpha$max_dens <- 0
        
        delta_plus[[i]]$beta$prop <- 0
        delta_plus[[i]]$beta$mean <- 0
        delta_plus[[i]]$beta$sd <- 0
        delta_plus[[i]]$beta$max_dens <- 0
        
        delta_plus[[i]]$gamma$prop <- 0
        delta_plus[[i]]$gamma$mean <- 0
        delta_plus[[i]]$gamma$sd <- 0
        delta_plus[[i]]$gamma$max_dens <- 0
        
        delta_plus[[i]]$epsilon$prop <- 0
        delta_plus[[i]]$epsilon$mean <- 0
        delta_plus[[i]]$epsilon$sd <- 0
        delta_plus[[i]]$epsilon$max_dens <- 0
        
        delta_plus[[i]]$zeta$prop <- 0
        delta_plus[[i]]$zeta$mean <- 0
        delta_plus[[i]]$zeta$sd <- 0
        delta_plus[[i]]$zeta$max_dens <- 0
        
        delta_plus[[i]]$chi$prop <- 0
        delta_plus[[i]]$chi$mean <- 0
        delta_plus[[i]]$chi$sd <- 0
        delta_plus[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        delta_plus[[i]]$total <- n / tot
        
        
      }
    }
    
    return(delta_plus)
  }
  
  
  ####For chi_plus
  if(angle == "chi_plus"){
    chi_plus <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"chi_plus"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for chi_plus analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        chi_plus[[i]]$alpha$prop <- gaussians_alpha$prop
        chi_plus[[i]]$alpha$mean <- gaussians_alpha$mean
        chi_plus[[i]]$alpha$sd <- gaussians_alpha$sd
        chi_plus[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for chi_plus
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        chi_plus[[i]]$beta$prop <- gaussians_beta$prop
        chi_plus[[i]]$beta$mean <- gaussians_beta$mean
        chi_plus[[i]]$beta$sd <- gaussians_beta$sd
        chi_plus[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for chi_plus
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        chi_plus[[i]]$gamma$prop <- gaussians_gamma$prop
        chi_plus[[i]]$gamma$mean <- gaussians_gamma$mean
        chi_plus[[i]]$gamma$sd <- gaussians_gamma$sd
        chi_plus[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for chi_plus
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        chi_plus[[i]]$delta$prop <- gaussians_delta$prop
        chi_plus[[i]]$delta$mean <- gaussians_delta$mean
        chi_plus[[i]]$delta$sd <- gaussians_delta$sd
        chi_plus[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for chi_plus
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        chi_plus[[i]]$epsilon$prop <- gaussians_epsilon$prop
        chi_plus[[i]]$epsilon$mean <- gaussians_epsilon$mean
        chi_plus[[i]]$epsilon$sd <- gaussians_epsilon$sd
        chi_plus[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for chi_plus
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        chi_plus[[i]]$zeta$prop <- gaussians_zeta$prop
        chi_plus[[i]]$zeta$mean <- gaussians_zeta$mean
        chi_plus[[i]]$zeta$sd <- gaussians_zeta$sd
        chi_plus[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        n <- length(index)
        chi_plus[[i]]$total <- n / tot
        
        
      }else{
        
        chi_plus[[i]]$alpha$prop <- 0
        chi_plus[[i]]$alpha$mean <- 0
        chi_plus[[i]]$alpha$sd <- 0
        chi_plus[[i]]$alpha$max_dens <- 0
        
        chi_plus[[i]]$beta$prop <- 0
        chi_plus[[i]]$beta$mean <- 0
        chi_plus[[i]]$beta$sd <- 0
        chi_plus[[i]]$beta$max_dens <- 0
        
        chi_plus[[i]]$gamma$prop <- 0
        chi_plus[[i]]$gamma$mean <- 0
        chi_plus[[i]]$gamma$sd <- 0
        chi_plus[[i]]$gamma$max_dens <- 0
        
        chi_plus[[i]]$delta$prop <- 0
        chi_plus[[i]]$delta$mean <- 0
        chi_plus[[i]]$delta$sd <- 0
        chi_plus[[i]]$delta$max_dens <- 0
        
        chi_plus[[i]]$epsilon$prop <- 0
        chi_plus[[i]]$epsilon$mean <- 0
        chi_plus[[i]]$epsilon$sd <- 0
        chi_plus[[i]]$epsilon$max_dens <- 0
        
        chi_plus[[i]]$zeta$prop <- 0
        chi_plus[[i]]$zeta$mean <- 0
        chi_plus[[i]]$zeta$sd <- 0
        chi_plus[[i]]$zeta$max_dens <- 0
        
        n <- length(index)
        chi_plus[[i]]$total <- n / tot
        
        
      }
    }
    
    return(chi_plus)
  }
  
  if(angle == "beta_minus"){
    beta_minus <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"beta_minus"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for beta_minus analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        beta_minus[[i]]$alpha$prop <- gaussians_alpha$prop
        beta_minus[[i]]$alpha$mean <- gaussians_alpha$mean
        beta_minus[[i]]$alpha$sd <- gaussians_alpha$sd
        beta_minus[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for gamma
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        beta_minus[[i]]$gamma$prop <- gaussians_gamma$prop
        beta_minus[[i]]$gamma$mean <- gaussians_gamma$mean
        beta_minus[[i]]$gamma$sd <- gaussians_gamma$sd
        beta_minus[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for delta
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        beta_minus[[i]]$delta$prop <- gaussians_delta$prop
        beta_minus[[i]]$delta$mean <- gaussians_delta$mean
        beta_minus[[i]]$delta$sd <- gaussians_delta$sd
        beta_minus[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for epsilon
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        beta_minus[[i]]$epsilon$prop <- gaussians_epsilon$prop
        beta_minus[[i]]$epsilon$mean <- gaussians_epsilon$mean
        beta_minus[[i]]$epsilon$sd <- gaussians_epsilon$sd
        beta_minus[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for zeta
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        beta_minus[[i]]$zeta$prop <- gaussians_zeta$prop
        beta_minus[[i]]$zeta$mean <- gaussians_zeta$mean
        beta_minus[[i]]$zeta$sd <- gaussians_zeta$sd
        beta_minus[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        beta_minus[[i]]$chi$prop <- gaussians_chi$prop
        beta_minus[[i]]$chi$mean <- gaussians_chi$mean
        beta_minus[[i]]$chi$sd <- gaussians_chi$sd
        beta_minus[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        beta_minus[[i]]$total <- n / tot
        
        
      }else{
        
        beta_minus[[i]]$alpha$prop <- 0
        beta_minus[[i]]$alpha$mean <- 0
        beta_minus[[i]]$alpha$sd <- 0
        beta_minus[[i]]$alpha$max_dens <- 0
        
        beta_minus[[i]]$gamma$prop <- 0
        beta_minus[[i]]$gamma$mean <- 0
        beta_minus[[i]]$gamma$sd <- 0
        beta_minus[[i]]$gamma$max_dens <- 0
        
        beta_minus[[i]]$delta$prop <- 0
        beta_minus[[i]]$delta$mean <- 0
        beta_minus[[i]]$delta$sd <- 0
        beta_minus[[i]]$delta$max_dens <- 0
        
        beta_minus[[i]]$epsilon$prop <- 0
        beta_minus[[i]]$epsilon$mean <- 0
        beta_minus[[i]]$epsilon$sd <- 0
        beta_minus[[i]]$epsilon$max_dens <- 0
        
        beta_minus[[i]]$zeta$prop <- 0
        beta_minus[[i]]$zeta$mean <- 0
        beta_minus[[i]]$zeta$sd <- 0
        beta_minus[[i]]$zeta$max_dens <- 0
        
        beta_minus[[i]]$chi$prop <- 0
        beta_minus[[i]]$chi$mean <- 0
        beta_minus[[i]]$chi$sd <- 0
        beta_minus[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        beta_minus[[i]]$total <- n / tot
        
        
        
      }
    }
    
    return(beta_minus)
  }
  
  
  #Check for gamma_minus
  if(angle == "gamma_minus"){
    gamma_minus <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"gamma_minus"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for gamma_minus analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        gamma_minus[[i]]$alpha$prop <- gaussians_alpha$prop
        gamma_minus[[i]]$alpha$mean <- gaussians_alpha$mean
        gamma_minus[[i]]$alpha$sd <- gaussians_alpha$sd
        gamma_minus[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for gamma_minus
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        gamma_minus[[i]]$beta$prop <- gaussians_beta$prop
        gamma_minus[[i]]$beta$mean <- gaussians_beta$mean
        gamma_minus[[i]]$beta$sd <- gaussians_beta$sd
        gamma_minus[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for delta
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        gamma_minus[[i]]$delta$prop <- gaussians_delta$prop
        gamma_minus[[i]]$delta$mean <- gaussians_delta$mean
        gamma_minus[[i]]$delta$sd <- gaussians_delta$sd
        gamma_minus[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for epsilon
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        gamma_minus[[i]]$epsilon$prop <- gaussians_epsilon$prop
        gamma_minus[[i]]$epsilon$mean <- gaussians_epsilon$mean
        gamma_minus[[i]]$epsilon$sd <- gaussians_epsilon$sd
        gamma_minus[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for zeta
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        gamma_minus[[i]]$zeta$prop <- gaussians_zeta$prop
        gamma_minus[[i]]$zeta$mean <- gaussians_zeta$mean
        gamma_minus[[i]]$zeta$sd <- gaussians_zeta$sd
        gamma_minus[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        gamma_minus[[i]]$chi$prop <- gaussians_chi$prop
        gamma_minus[[i]]$chi$mean <- gaussians_chi$mean
        gamma_minus[[i]]$chi$sd <- gaussians_chi$sd
        gamma_minus[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        gamma_minus[[i]]$total <- n / tot
        
        
      }else{
        
        gamma_minus[[i]]$alpha$prop <- 0
        gamma_minus[[i]]$alpha$mean <- 0
        gamma_minus[[i]]$alpha$sd <- 0
        gamma_minus[[i]]$alpha$max_dens <- 0
        
        gamma_minus[[i]]$beta$prop <- 0
        gamma_minus[[i]]$beta$mean <- 0
        gamma_minus[[i]]$beta$sd <- 0
        gamma_minus[[i]]$beta$max_dens <- 0
        
        gamma_minus[[i]]$delta$prop <- 0
        gamma_minus[[i]]$delta$mean <- 0
        gamma_minus[[i]]$delta$sd <- 0
        gamma_minus[[i]]$delta$max_dens <- 0
        
        gamma_minus[[i]]$epsilon$prop <- 0
        gamma_minus[[i]]$epsilon$mean <- 0
        gamma_minus[[i]]$epsilon$sd <- 0
        gamma_minus[[i]]$epsilon$max_dens <- 0
        
        gamma_minus[[i]]$zeta$prop <- 0
        gamma_minus[[i]]$zeta$mean <- 0
        gamma_minus[[i]]$zeta$sd <- 0
        gamma_minus[[i]]$zeta$max_dens <- 0
        
        gamma_minus[[i]]$chi$prop <- 0
        gamma_minus[[i]]$chi$mean <- 0
        gamma_minus[[i]]$chi$sd <- 0
        gamma_minus[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        gamma_minus[[i]]$total <- n / tot
        
        
        
      }
    }
    
    return(gamma_minus)
  }
  
  
  #####For delta_minus
  if(angle == "delta_minus"){
    delta_minus <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"delta_minus"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for delta_minus analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        delta_minus[[i]]$alpha$prop <- gaussians_alpha$prop
        delta_minus[[i]]$alpha$mean <- gaussians_alpha$mean
        delta_minus[[i]]$alpha$sd <- gaussians_alpha$sd
        delta_minus[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for beta
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        delta_minus[[i]]$beta$prop <- gaussians_beta$prop
        delta_minus[[i]]$beta$mean <- gaussians_beta$mean
        delta_minus[[i]]$beta$sd <- gaussians_beta$sd
        delta_minus[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for gamma
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        delta_minus[[i]]$gamma$prop <- gaussians_gamma$prop
        delta_minus[[i]]$gamma$mean <- gaussians_gamma$mean
        delta_minus[[i]]$gamma$sd <- gaussians_gamma$sd
        delta_minus[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for epsilon
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        delta_minus[[i]]$epsilon$prop <- gaussians_epsilon$prop
        delta_minus[[i]]$epsilon$mean <- gaussians_epsilon$mean
        delta_minus[[i]]$epsilon$sd <- gaussians_epsilon$sd
        delta_minus[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for zeta
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        delta_minus[[i]]$zeta$prop <- gaussians_zeta$prop
        delta_minus[[i]]$zeta$mean <- gaussians_zeta$mean
        delta_minus[[i]]$zeta$sd <- gaussians_zeta$sd
        delta_minus[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        delta_minus[[i]]$chi$prop <- gaussians_chi$prop
        delta_minus[[i]]$chi$mean <- gaussians_chi$mean
        delta_minus[[i]]$chi$sd <- gaussians_chi$sd
        delta_minus[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        delta_minus[[i]]$total <- n / tot
        
        
      }else{
        
        delta_minus[[i]]$alpha$prop <- 0
        delta_minus[[i]]$alpha$mean <- 0
        delta_minus[[i]]$alpha$sd <- 0
        delta_minus[[i]]$alpha$max_dens <- 0
        
        delta_minus[[i]]$beta$prop <- 0
        delta_minus[[i]]$beta$mean <- 0
        delta_minus[[i]]$beta$sd <- 0
        delta_minus[[i]]$beta$max_dens <- 0
        
        delta_minus[[i]]$gamma$prop <- 0
        delta_minus[[i]]$gamma$mean <- 0
        delta_minus[[i]]$gamma$sd <- 0
        delta_minus[[i]]$gamma$max_dens <- 0
        
        delta_minus[[i]]$epsilon$prop <- 0
        delta_minus[[i]]$epsilon$mean <- 0
        delta_minus[[i]]$epsilon$sd <- 0
        delta_minus[[i]]$epsilon$max_dens <- 0
        
        delta_minus[[i]]$zeta$prop <- 0
        delta_minus[[i]]$zeta$mean <- 0
        delta_minus[[i]]$zeta$sd <- 0
        delta_minus[[i]]$zeta$max_dens <- 0
        
        delta_minus[[i]]$chi$prop <- 0
        delta_minus[[i]]$chi$mean <- 0
        delta_minus[[i]]$chi$sd <- 0
        delta_minus[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        delta_minus[[i]]$total <- n / tot
        
        
      }
    }
    
    return(delta_minus)
  }
  
  ##### For epsilon_minus
  if(angle == "epsilon_minus"){
    epsilon_minus <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"epsilon_minus"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for epsilon_minus analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        epsilon_minus[[i]]$alpha$prop <- gaussians_alpha$prop
        epsilon_minus[[i]]$alpha$mean <- gaussians_alpha$mean
        epsilon_minus[[i]]$alpha$sd <- gaussians_alpha$sd
        epsilon_minus[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for epsilon_minus
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        epsilon_minus[[i]]$beta$prop <- gaussians_beta$prop
        epsilon_minus[[i]]$beta$mean <- gaussians_beta$mean
        epsilon_minus[[i]]$beta$sd <- gaussians_beta$sd
        epsilon_minus[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for epsilon_minus
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        epsilon_minus[[i]]$gamma$prop <- gaussians_gamma$prop
        epsilon_minus[[i]]$gamma$mean <- gaussians_gamma$mean
        epsilon_minus[[i]]$gamma$sd <- gaussians_gamma$sd
        epsilon_minus[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for epsilon_minus
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        epsilon_minus[[i]]$delta$prop <- gaussians_delta$prop
        epsilon_minus[[i]]$delta$mean <- gaussians_delta$mean
        epsilon_minus[[i]]$delta$sd <- gaussians_delta$sd
        epsilon_minus[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for zeta
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        epsilon_minus[[i]]$zeta$prop <- gaussians_zeta$prop
        epsilon_minus[[i]]$zeta$mean <- gaussians_zeta$mean
        epsilon_minus[[i]]$zeta$sd <- gaussians_zeta$sd
        epsilon_minus[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        epsilon_minus[[i]]$chi$prop <- gaussians_chi$prop
        epsilon_minus[[i]]$chi$mean <- gaussians_chi$mean
        epsilon_minus[[i]]$chi$sd <- gaussians_chi$sd
        epsilon_minus[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        epsilon_minus[[i]]$total <- n / tot
        
        
      }else{
        
        epsilon_minus[[i]]$alpha$prop <- 0
        epsilon_minus[[i]]$alpha$mean <- 0
        epsilon_minus[[i]]$alpha$sd <- 0
        epsilon_minus[[i]]$alpha$max_dens <- 0
        
        epsilon_minus[[i]]$beta$prop <- 0
        epsilon_minus[[i]]$beta$mean <- 0
        epsilon_minus[[i]]$beta$sd <- 0
        epsilon_minus[[i]]$beta$max_dens <- 0
        
        epsilon_minus[[i]]$gamma$prop <- 0
        epsilon_minus[[i]]$gamma$mean <- 0
        epsilon_minus[[i]]$gamma$sd <- 0
        epsilon_minus[[i]]$gamma$max_dens <- 0
        
        epsilon_minus[[i]]$delta$prop <- 0
        epsilon_minus[[i]]$delta$mean <- 0
        epsilon_minus[[i]]$delta$sd <- 0
        epsilon_minus[[i]]$delta$max_dens <- 0
        
        epsilon_minus[[i]]$zeta$prop <- 0
        epsilon_minus[[i]]$zeta$mean <- 0
        epsilon_minus[[i]]$zeta$sd <- 0
        epsilon_minus[[i]]$zeta$max_dens <- 0
        
        epsilon_minus[[i]]$chi$prop <- 0
        epsilon_minus[[i]]$chi$mean <- 0
        epsilon_minus[[i]]$chi$sd <- 0
        epsilon_minus[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        epsilon_minus[[i]]$total <- n / tot
        
        
      }
    }
    
    return(epsilon_minus)
  }
  
  #####For zeta_minus
  if(angle == "zeta_minus"){
    zeta_minus <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"zeta_minus"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for zeta_minus analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        zeta_minus[[i]]$alpha$prop <- gaussians_alpha$prop
        zeta_minus[[i]]$alpha$mean <- gaussians_alpha$mean
        zeta_minus[[i]]$alpha$sd <- gaussians_alpha$sd
        zeta_minus[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for zeta_minus
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        zeta_minus[[i]]$beta$prop <- gaussians_beta$prop
        zeta_minus[[i]]$beta$mean <- gaussians_beta$mean
        zeta_minus[[i]]$beta$sd <- gaussians_beta$sd
        zeta_minus[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for zeta_minus
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        zeta_minus[[i]]$gamma$prop <- gaussians_gamma$prop
        zeta_minus[[i]]$gamma$mean <- gaussians_gamma$mean
        zeta_minus[[i]]$gamma$sd <- gaussians_gamma$sd
        zeta_minus[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for zeta_minus
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        zeta_minus[[i]]$delta$prop <- gaussians_delta$prop
        zeta_minus[[i]]$delta$mean <- gaussians_delta$mean
        zeta_minus[[i]]$delta$sd <- gaussians_delta$sd
        zeta_minus[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for zeta_minus
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        zeta_minus[[i]]$epsilon$prop <- gaussians_epsilon$prop
        zeta_minus[[i]]$epsilon$mean <- gaussians_epsilon$mean
        zeta_minus[[i]]$epsilon$sd <- gaussians_epsilon$sd
        zeta_minus[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for chi
        gaussians_chi<- fit_data(na.omit(dataframe2_no_disc$chi))
        zeta_minus[[i]]$chi$prop <- gaussians_chi$prop
        zeta_minus[[i]]$chi$mean <- gaussians_chi$mean
        zeta_minus[[i]]$chi$sd <- gaussians_chi$sd
        zeta_minus[[i]]$chi$max_dens <- gaussians_chi$max_dens
        
        n <- length(index)
        zeta_minus[[i]]$total <- n / tot
        
        
      }else{
        
        zeta_minus[[i]]$alpha$prop <- 0
        zeta_minus[[i]]$alpha$mean <- 0
        zeta_minus[[i]]$alpha$sd <- 0
        zeta_minus[[i]]$alpha$max_dens <- 0
        
        zeta_minus[[i]]$beta$prop <- 0
        zeta_minus[[i]]$beta$mean <- 0
        zeta_minus[[i]]$beta$sd <- 0
        zeta_minus[[i]]$beta$max_dens <- 0
        
        zeta_minus[[i]]$gamma$prop <- 0
        zeta_minus[[i]]$gamma$mean <- 0
        zeta_minus[[i]]$gamma$sd <- 0
        zeta_minus[[i]]$gamma$max_dens <- 0
        
        zeta_minus[[i]]$delta$prop <- 0
        zeta_minus[[i]]$delta$mean <- 0
        zeta_minus[[i]]$delta$sd <- 0
        zeta_minus[[i]]$delta$max_dens <- 0
        
        zeta_minus[[i]]$epsilon$prop <- 0
        zeta_minus[[i]]$epsilon$mean <- 0
        zeta_minus[[i]]$epsilon$sd <- 0
        zeta_minus[[i]]$epsilon$max_dens <- 0
        
        zeta_minus[[i]]$chi$prop <- 0
        zeta_minus[[i]]$chi$mean <- 0
        zeta_minus[[i]]$chi$sd <- 0
        zeta_minus[[i]]$chi$max_dens <- 0
        
        n <- length(index)
        zeta_minus[[i]]$total <- n / tot
        
        
      }
    }
    
    return(zeta_minus)
  }
  
  ####For chi_minus
  if(angle == "chi_minus"){
    chi_minus <- vector("list", length = length(inter_list))
    
    for(i in 1:length(inter_list)){
      index <- which(ntinfo_m[,"chi_minus"] == i)
      if(length(index) > 5){
        dataframe2 <- ntinfo_m[index, ]
        dataframe2_no_disc <- ntinfo[index,]
        
        #Check for chi_minus analysis
        gaussians_alpha<- fit_data(na.omit(dataframe2_no_disc$alpha))
        chi_minus[[i]]$alpha$prop <- gaussians_alpha$prop
        chi_minus[[i]]$alpha$mean <- gaussians_alpha$mean
        chi_minus[[i]]$alpha$sd <- gaussians_alpha$sd
        chi_minus[[i]]$alpha$max_dens <- gaussians_alpha$max_dens
        
        #Check for chi_minus
        gaussians_beta<- fit_data(na.omit(dataframe2_no_disc$beta))
        chi_minus[[i]]$beta$prop <- gaussians_beta$prop
        chi_minus[[i]]$beta$mean <- gaussians_beta$mean
        chi_minus[[i]]$beta$sd <- gaussians_beta$sd
        chi_minus[[i]]$beta$max_dens <- gaussians_beta$max_dens
        
        #Check for chi_minus
        gaussians_gamma<- fit_data(na.omit(dataframe2_no_disc$gamma))
        chi_minus[[i]]$gamma$prop <- gaussians_gamma$prop
        chi_minus[[i]]$gamma$mean <- gaussians_gamma$mean
        chi_minus[[i]]$gamma$sd <- gaussians_gamma$sd
        chi_minus[[i]]$gamma$max_dens <- gaussians_gamma$max_dens
        
        #Check for chi_minus
        gaussians_delta<- fit_data(na.omit(dataframe2_no_disc$delta))
        chi_minus[[i]]$delta$prop <- gaussians_delta$prop
        chi_minus[[i]]$delta$mean <- gaussians_delta$mean
        chi_minus[[i]]$delta$sd <- gaussians_delta$sd
        chi_minus[[i]]$delta$max_dens <- gaussians_delta$max_dens
        
        #Check for chi_minus
        gaussians_epsilon<- fit_data(na.omit(dataframe2_no_disc$epsilon))
        chi_minus[[i]]$epsilon$prop <- gaussians_epsilon$prop
        chi_minus[[i]]$epsilon$mean <- gaussians_epsilon$mean
        chi_minus[[i]]$epsilon$sd <- gaussians_epsilon$sd
        chi_minus[[i]]$epsilon$max_dens <- gaussians_epsilon$max_dens
        
        #Check for chi_minus
        gaussians_zeta<- fit_data(na.omit(dataframe2_no_disc$zeta))
        chi_minus[[i]]$zeta$prop <- gaussians_zeta$prop
        chi_minus[[i]]$zeta$mean <- gaussians_zeta$mean
        chi_minus[[i]]$zeta$sd <- gaussians_zeta$sd
        chi_minus[[i]]$zeta$max_dens <- gaussians_zeta$max_dens
        
        n <- length(index)
        chi_minus[[i]]$total <- n / tot
        
        
      }else{
        
        chi_minus[[i]]$alpha$prop <- 0
        chi_minus[[i]]$alpha$mean <- 0
        chi_minus[[i]]$alpha$sd <- 0
        chi_minus[[i]]$alpha$max_dens <- 0
        
        chi_minus[[i]]$beta$prop <- 0
        chi_minus[[i]]$beta$mean <- 0
        chi_minus[[i]]$beta$sd <- 0
        chi_minus[[i]]$beta$max_dens <- 0
        
        chi_minus[[i]]$gamma$prop <- 0
        chi_minus[[i]]$gamma$mean <- 0
        chi_minus[[i]]$gamma$sd <- 0
        chi_minus[[i]]$gamma$max_dens <- 0
        
        chi_minus[[i]]$delta$prop <- 0
        chi_minus[[i]]$delta$mean <- 0
        chi_minus[[i]]$delta$sd <- 0
        chi_minus[[i]]$delta$max_dens <- 0
        
        chi_minus[[i]]$epsilon$prop <- 0
        chi_minus[[i]]$epsilon$mean <- 0
        chi_minus[[i]]$epsilon$sd <- 0
        chi_minus[[i]]$epsilon$max_dens <- 0
        
        chi_minus[[i]]$zeta$prop <- 0
        chi_minus[[i]]$zeta$mean <- 0
        chi_minus[[i]]$zeta$sd <- 0
        chi_minus[[i]]$zeta$max_dens <- 0
        
        n <- length(index)
        chi_minus[[i]]$total <- n / tot
        
        
      }
    }
    
    return(chi_minus)
  }
  
}



