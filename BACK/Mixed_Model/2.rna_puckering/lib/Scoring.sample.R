scoring_sampling <- function(ntinfo_disc= NULL, bang = bayesian_total, ntinfo = NULL){
  
  #windows <- bd$windows
  #alpha <- bd$alpha
  #beta <- bd$beta
  #gamma <- bd$gamma
  #delta <- bd$delta
  #epsilon <- bd$epsilon
  #zeta <- bd$zeta
  #chi <- bd$chi
  
  #eval <- minipack::discretize(ntinfo_m = evaluate, windows = windows)
  eval <- ntinfo_disc
  eval2 <- ntinfo
  
  ind <-which(eval %in% 0)
  if(length(ind) > 0){
    eval[ind] <- 1
  }
  
  percent <- c()
  
  alpha_amb <- c(nulling(neotest("beta",eval[,"beta"],"alpha", eval2[,"alpha"], bayesian_total = bang)),
                 nulling(neotest("gamma",eval[,"gamma"],"alpha",eval2[,"alpha"], bayesian_total = bang)),
                 nulling(neotest("delta",eval[,"delta"],"alpha",eval2[,"alpha"], bayesian_total = bang)),
                 nulling(neotest("epsilon",eval[,"epsilon"],"alpha",eval2[,"alpha"], bayesian_total = bang)),
                 nulling(neotest("zeta",eval[,"zeta"],"alpha",eval2[,"alpha"], bayesian_total = bang)),
                 nulling(neotest("chi",eval[,"chi"],"alpha",eval2[,"alpha"], bayesian_total = bang)))
  
  beta_amb <- c(nulling(neotest("alpha",eval[,"alpha"],"beta",eval2[,"beta"], bayesian_total = bang)),
                nulling(neotest("gamma",eval[,"gamma"],"beta",eval2[,"beta"], bayesian_total = bang)),
                nulling(neotest("delta",eval[,"delta"],"beta",eval2[,"beta"], bayesian_total = bang)),
                nulling(neotest("epsilon",eval[,"epsilon"],"beta",eval2[,"beta"], bayesian_total = bang)),
                nulling(neotest("zeta",eval[,"zeta"],"beta",eval2[,"beta"], bayesian_total = bang)),
                nulling(neotest("chi",eval[,"chi"],"beta",eval2[,"beta"], bayesian_total = bang)))
  
  gamma_amb <- c(nulling(neotest("alpha",eval[,"alpha"],"gamma",eval2[,"gamma"], bayesian_total = bang)),
                 nulling(neotest("beta",eval[,"beta"],"gamma",eval2[,"gamma"], bayesian_total = bang)),
                 nulling(neotest("delta",eval[,"delta"],"gamma",eval2[,"gamma"], bayesian_total = bang)), 
                 nulling(neotest("epsilon",eval[,"epsilon"],"gamma",eval2[,"gamma"], bayesian_total = bang)),
                 nulling(neotest("zeta",eval[,"zeta"],"gamma",eval2[,"gamma"], bayesian_total = bang)),
                 nulling(neotest("chi",eval[,"chi"],"gamma",eval2[,"gamma"], bayesian_total = bang)))
  
  delta_amb <- c(nulling(neotest("alpha",eval[,"alpha"],"delta",eval2[,"delta"], bayesian_total = bang)),
                 nulling(neotest("beta",eval[,"beta"],"delta",eval2[,"delta"], bayesian_total = bang)),
                 nulling(neotest("gamma",eval[,"gamma"],"delta",eval2[,"delta"], bayesian_total = bang)),
                 nulling(neotest("epsilon",eval[,"epsilon"],"delta",eval2[,"delta"], bayesian_total = bang)),
                 nulling(neotest("zeta",eval[,"zeta"],"delta",eval2[,"delta"], bayesian_total = bang)),
                 nulling(neotest("chi",eval[,"chi"],"delta",eval2[,"delta"], bayesian_total = bang)))
  
  epsilon_amb <- c(nulling(neotest("alpha",eval[,"alpha"],"epsilon",eval2[,"epsilon"], bayesian_total = bang)),
                   nulling(neotest("beta",eval[,"beta"],"epsilon",eval2[,"epsilon"], bayesian_total = bang)),
                   nulling(neotest("gamma",eval[,"gamma"],"epsilon",eval2[,"epsilon"], bayesian_total = bang)),
                   nulling(neotest("delta",eval[,"delta"],"epsilon",eval2[,"epsilon"], bayesian_total = bang)),
                   nulling(neotest("zeta",eval[,"zeta"],"epsilon",eval2[,"epsilon"], bayesian_total = bang)),
                   nulling(neotest("chi",eval[,"chi"],"epsilon",eval2[,"epsilon"], bayesian_total = bang)))
  
  zeta_amb <- c(nulling(neotest("alpha",eval[,"alpha"],"zeta",eval2[,"zeta"], bayesian_total = bang)),
                nulling(neotest("beta",eval[,"beta"],"zeta",eval2[,"zeta"], bayesian_total = bang)),
                nulling(neotest("gamma",eval[,"gamma"],"zeta",eval2[,"zeta"], bayesian_total = bang)),
                nulling(neotest("delta",eval[,"delta"],"zeta",eval2[,"zeta"], bayesian_total = bang)),
                nulling(neotest("epsilon",eval[,"epsilon"],"zeta",eval2[,"zeta"], bayesian_total = bang)),
                nulling(neotest("chi",eval[,"chi"],"zeta",eval2[,"zeta"], bayesian_total = bang)))
  
  chi_amb <- c(nulling(neotest("alpha",eval[,"alpha"],"chi",eval2[,"chi"], bayesian_total = bang)),
               nulling(neotest("beta",eval[,"beta"],"chi",eval2[,"chi"], bayesian_total = bang)),
               nulling(neotest("gamma",eval[,"gamma"],"chi",eval2[,"chi"], bayesian_total = bang)),
               nulling(neotest("delta",eval[,"delta"],"chi",eval2[,"chi"], bayesian_total = bang)),
               nulling(neotest("epsilon",eval[,"epsilon"],"chi",eval2[,"chi"], bayesian_total = bang)),
               nulling(neotest("zeta",eval[,"zeta"],"chi",eval2[,"chi"], bayesian_total = bang)))
  
  percent$alpha <- sum(alpha_amb) / length(alpha_amb)
  percent$beta <- sum(beta_amb) / length(beta_amb)
  percent$gamma <- sum(gamma_amb) / length(gamma_amb)
  percent$delta <- sum(delta_amb) / length(delta_amb)
  percent$epsilon <- sum(epsilon_amb) / length(epsilon_amb)
  percent$zeta <- sum(zeta_amb) / length(zeta_amb)
  percent$chi <- sum(chi_amb) / length(chi_amb)
  percent$total <- sum(c(percent$alpha, percent$beta, percent$gamma, percent$delta, percent$epsilon, percent$zeta, percent$chi)) / 7
  
  return(percent)
}
