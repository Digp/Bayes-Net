scoring_bayes <- function(ntinfo_disc= NULL, bd = bayesian_data, ntinfo = NULL){
  
  #windows <- bd$windows
  alpha <- bd$alpha
  beta <- bd$beta
  gamma <- bd$gamma
  delta <- bd$delta
  epsilon <- bd$epsilon
  zeta <- bd$zeta
  chi <- bd$chi
  alpha_plus <-bd$alpha_plus
  beta_plus <- bd$beta_plus
  gamma_plus <- bd$gamma_plus
  delta_plus <- bd$delta_plus
  chi_plus <- bd$chi_plus
  beta_minus <- bd$beta_minus
  gamma_minus <- bd$gamma_minus
  delta_minus <- bd$delta_minus
  epsilon_minus <- bd$epsilon_minus
  zeta_minus <-bd$zeta_minus
  chi_minus <- bd$chi_minus
  
  #eval <- minipack::discretize(ntinfo_m = evaluate, windows = windows)
  eval <- ntinfo_disc
  eval2 <- ntinfo
  
  bd <- bayesian_data
  
  ind <-which(eval %in% 0)
  if(length(ind) > 0){
    eval[ind] <- 1
  }
  
  percent <- c()
  
  alpha_amb <- c(nulling(dnorms_calc("beta",eval[,"beta"],"alpha", eval2[,"alpha"])),
                 nulling(dnorms_calc("gamma",eval[,"gamma"],"alpha",eval2[,"alpha"])),
                 nulling(dnorms_calc("delta_minus",eval[,"delta_minus"],"alpha",eval2[,"alpha"])),
                 nulling(dnorms_calc("delta",eval[,"delta"],"alpha",eval2[,"alpha"])),
                 nulling(dnorms_calc("epsilon_minus",eval[,"epsilon_minus"],"alpha",eval2[,"alpha"])),
                 nulling(dnorms_calc("zeta_minus",eval[,"zeta_minus"],"alpha",eval2[,"alpha"])),
                 nulling(dnorms_calc("chi",eval[,"chi"],"alpha",eval2[,"alpha"])))
  
  beta_amb <- c(nulling(dnorms_calc("alpha",eval[,"alpha"],"beta",eval2[,"beta"])),
                nulling(dnorms_calc("gamma",eval[,"gamma"],"beta",eval2[,"beta"])),
                nulling(dnorms_calc("delta_minus",eval[,"delta_minus"],"beta",eval2[,"beta"])),
                nulling(dnorms_calc("delta",eval[,"delta"],"beta",eval2[,"beta"])),
                nulling(dnorms_calc("epsilon_minus",eval[,"epsilon_minus"],"beta",eval2[,"beta"])),
                nulling(dnorms_calc("zeta_minus",eval[,"zeta_minus"],"beta",eval2[,"beta"])),
                nulling(dnorms_calc("chi",eval[,"chi"],"beta",eval2[,"beta"])))
  
  gamma_amb <- c(nulling(dnorms_calc("alpha",eval[,"alpha"],"gamma",eval2[,"gamma"])),
                 nulling(dnorms_calc("beta",eval[,"beta"],"gamma",eval2[,"gamma"])),
                 nulling(dnorms_calc("delta",eval[,"delta"],"gamma",eval2[,"gamma"])), 
                 nulling(dnorms_calc("delta_minus",eval[,"delta_minus"],"gamma",eval2[,"gamma"])),
                 nulling(dnorms_calc("epsilon",eval[,"epsilon"],"gamma",eval2[,"gamma"])),
                 nulling(dnorms_calc("epsilon_minus",eval[,"epsilon_minus"],"gamma",eval2[,"gamma"])),
                 nulling(dnorms_calc("zeta",eval[,"zeta"],"gamma",eval2[,"gamma"])),
                 nulling(dnorms_calc("zeta_minus",eval[,"zeta_minus"],"gamma",eval2[,"gamma"])),
                 nulling(dnorms_calc("chi",eval[,"chi"],"gamma",eval2[,"gamma"])))
  
  delta_amb <- c(nulling(dnorms_calc("alpha",eval[,"alpha"],"delta",eval2[,"delta"])),
                 nulling(dnorms_calc("beta",eval[,"beta"],"delta",eval2[,"delta"])),
                 nulling(dnorms_calc("gamma",eval[,"gamma"],"delta",eval2[,"delta"])),
                 nulling(dnorms_calc("epsilon",eval[,"epsilon"],"delta",eval2[,"delta"])),
                 nulling(dnorms_calc("zeta",eval[,"zeta"],"delta",eval2[,"delta"])),
                 nulling(dnorms_calc("chi",eval[,"chi"],"delta",eval2[,"delta"])))
  
  epsilon_amb <- c(nulling(dnorms_calc("alpha",eval[,"alpha"],"epsilon",eval2[,"epsilon"])),
                   nulling(dnorms_calc("alpha_plus",eval[,"alpha_plus"],"epsilon",eval2[,"epsilon"])),
                   nulling(dnorms_calc("beta",eval[,"beta"],"epsilon",eval2[,"epsilon"])),
                   nulling(dnorms_calc("beta_plus",eval[,"beta_plus"],"epsilon",eval2[,"epsilon"])),
                   nulling(dnorms_calc("gamma",eval[,"gamma"],"epsilon",eval2[,"epsilon"])),
                   nulling(dnorms_calc("gamma_plus",eval[,"gamma_plus"],"epsilon",eval2[,"epsilon"])),
                   nulling(dnorms_calc("delta",eval[,"delta"],"epsilon",eval2[,"epsilon"])),
                   nulling(dnorms_calc("delta_plus",eval[,"delta_plus"],"epsilon",eval2[,"epsilon"])),
                   nulling(dnorms_calc("zeta",eval[,"zeta"],"epsilon",eval2[,"epsilon"])))
  
  zeta_amb <- c(nulling(dnorms_calc("alpha_plus",eval[,"alpha_plus"],"zeta",eval2[,"zeta"])),
                nulling(dnorms_calc("beta_plus",eval[,"beta_plus"],"zeta",eval2[,"zeta"])),
                nulling(dnorms_calc("gamma_plus",eval[,"gamma_plus"],"zeta",eval2[,"zeta"])),
                nulling(dnorms_calc("delta",eval[,"delta"],"zeta",eval2[,"zeta"])),
                nulling(dnorms_calc("delta_plus",eval[,"delta_plus"],"zeta",eval2[,"zeta"])),
                nulling(dnorms_calc("epsilon",eval[,"epsilon"],"zeta",eval2[,"zeta"])),
                nulling(dnorms_calc("chi",eval[,"chi"],"zeta",eval2[,"zeta"])))
  
  chi_amb <- c(nulling(dnorms_calc("alpha",eval[,"alpha"],"chi",eval2[,"chi"])),
               nulling(dnorms_calc("beta",eval[,"beta"],"chi",eval2[,"chi"])),
               nulling(dnorms_calc("gamma",eval[,"gamma"],"chi",eval2[,"chi"])),
               nulling(dnorms_calc("delta",eval[,"delta"],"chi",eval2[,"chi"])),
               nulling(dnorms_calc("epsilon",eval[,"epsilon"],"chi",eval2[,"chi"])),
               nulling(dnorms_calc("zeta",eval[,"zeta"],"chi",eval2[,"chi"])))
  
  percent$alpha <- sum(alpha_amb) / length(alpha_amb)
  percent$beta <- sum(beta_amb) / length(beta_amb)
  percent$gamma <- sum(gamma_amb) / length(gamma_amb)
  percent$delta <- sum(delta_amb) / length(delta_amb)
  percent$epsilon <- sum(epsilon_amb) / length(epsilon_amb)
  percent$zeta <- sum(zeta_amb) / length(zeta_amb)
  percent$chi <- sum(chi_amb) / length(chi_amb)
  
  #Find how many 0s are found. 0s mean where original data angle is a NA
  ind_zero <- length(which(percent %in% 0))
  
  percent$total <- sum(c(percent$alpha, percent$beta, percent$gamma, percent$delta, percent$epsilon, percent$zeta, percent$chi)) / (7 - ind_zero)
  
  return(percent)
}
