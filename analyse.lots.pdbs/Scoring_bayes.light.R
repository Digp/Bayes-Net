scoring_bayes <- function(ntinfo_disc= NULL, bd = bayesian_data, ntinfo = NULL){
  
  #windows <- bd$windows
  alpha <- bd$alpha
  beta <- bd$beta
  gamma <- bd$gamma
  delta <- bd$delta
  epsilon <- bd$epsilon
  zeta <- bd$zeta
  chi <- bd$chi
  
  #eval <- minipack::discretize(ntinfo_m = evaluate, windows = windows)
  eval <- ntinfo_disc
  eval2 <- ntinfo
  
  ind <-which(eval %in% 0)
  if(length(ind) > 0){
    eval[ind] <- 1
  }
  

  percent <- c(sum(c(nulling(dnorms_calc("beta",eval[,"beta"],"alpha", eval2[,"alpha"])),
                 nulling(dnorms_calc("gamma",eval[,"gamma"],"alpha",eval2[,"alpha"])),
                 nulling(dnorms_calc("delta",eval[,"delta"],"alpha",eval2[,"alpha"])),
                 nulling(dnorms_calc("epsilon",eval[,"epsilon"],"alpha",eval2[,"alpha"])),
                 nulling(dnorms_calc("zeta",eval[,"zeta"],"alpha",eval2[,"alpha"])),
                 nulling(dnorms_calc("chi",eval[,"chi"],"alpha",eval2[,"alpha"])))) / 6,
  
   sum(c(nulling(dnorms_calc("alpha",eval[,"alpha"],"beta",eval2[,"beta"])),
                nulling(dnorms_calc("gamma",eval[,"gamma"],"beta",eval2[,"beta"])),
                nulling(dnorms_calc("delta",eval[,"delta"],"beta",eval2[,"beta"])),
                nulling(dnorms_calc("epsilon",eval[,"epsilon"],"beta",eval2[,"beta"])),
                nulling(dnorms_calc("zeta",eval[,"zeta"],"beta",eval2[,"beta"])),
                nulling(dnorms_calc("chi",eval[,"chi"],"beta",eval2[,"beta"])))) / 6,
  
  sum(c(nulling(dnorms_calc("alpha",eval[,"alpha"],"gamma",eval2[,"gamma"])),
                 nulling(dnorms_calc("beta",eval[,"beta"],"gamma",eval2[,"gamma"])),
                 nulling(dnorms_calc("delta",eval[,"delta"],"gamma",eval2[,"gamma"])), 
                 nulling(dnorms_calc("epsilon",eval[,"epsilon"],"gamma",eval2[,"gamma"])),
                 nulling(dnorms_calc("zeta",eval[,"zeta"],"gamma",eval2[,"gamma"])),
                 nulling(dnorms_calc("chi",eval[,"chi"],"gamma",eval2[,"gamma"])))) / 6,
  
  sum(c(nulling(dnorms_calc("alpha",eval[,"alpha"],"delta",eval2[,"delta"])),
                 nulling(dnorms_calc("beta",eval[,"beta"],"delta",eval2[,"delta"])),
                 nulling(dnorms_calc("gamma",eval[,"gamma"],"delta",eval2[,"delta"])),
                 nulling(dnorms_calc("epsilon",eval[,"epsilon"],"delta",eval2[,"delta"])),
                 nulling(dnorms_calc("zeta",eval[,"zeta"],"delta",eval2[,"delta"])),
                 nulling(dnorms_calc("chi",eval[,"chi"],"delta",eval2[,"delta"])))) / 6,
  
  sum(c(nulling(dnorms_calc("alpha",eval[,"alpha"],"epsilon",eval2[,"epsilon"])),
                   nulling(dnorms_calc("beta",eval[,"beta"],"epsilon",eval2[,"epsilon"])),
                   nulling(dnorms_calc("gamma",eval[,"gamma"],"epsilon",eval2[,"epsilon"])),
                   nulling(dnorms_calc("delta",eval[,"delta"],"epsilon",eval2[,"epsilon"])),
                   nulling(dnorms_calc("zeta",eval[,"zeta"],"epsilon",eval2[,"epsilon"])),
                   nulling(dnorms_calc("chi",eval[,"chi"],"epsilon",eval2[,"epsilon"])))) / 6,
  
  sum(c(nulling(dnorms_calc("alpha",eval[,"alpha"],"zeta",eval2[,"zeta"])),
                nulling(dnorms_calc("beta",eval[,"beta"],"zeta",eval2[,"zeta"])),
                nulling(dnorms_calc("gamma",eval[,"gamma"],"zeta",eval2[,"zeta"])),
                nulling(dnorms_calc("delta",eval[,"delta"],"zeta",eval2[,"zeta"])),
                nulling(dnorms_calc("epsilon",eval[,"epsilon"],"zeta",eval2[,"zeta"])),
                nulling(dnorms_calc("chi",eval[,"chi"],"zeta",eval2[,"zeta"])))) / 6,
  
  sum(c(nulling(dnorms_calc("alpha",eval[,"alpha"],"chi",eval2[,"chi"])),
               nulling(dnorms_calc("beta",eval[,"beta"],"chi",eval2[,"chi"])),
               nulling(dnorms_calc("gamma",eval[,"gamma"],"chi",eval2[,"chi"])),
               nulling(dnorms_calc("delta",eval[,"delta"],"chi",eval2[,"chi"])),
               nulling(dnorms_calc("epsilon",eval[,"epsilon"],"chi",eval2[,"chi"])),
               nulling(dnorms_calc("zeta",eval[,"zeta"],"chi",eval2[,"chi"])))) / 6)
  
  
  #Find how many 0s are found. 0s mean where original data angle is a NA
  #ind_zero <- length(which(percent %in% 0))
  
  #percent$total <- sum(c(percent$alpha, percent$beta, percent$gamma, percent$delta, percent$epsilon, percent$zeta, percent$chi)) / (7 - ind_zero)
  
  return(sum(percent) / (7 -length(which(percent %in% 0))))
}
