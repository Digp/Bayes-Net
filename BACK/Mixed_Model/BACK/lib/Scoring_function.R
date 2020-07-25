scoring <- function(ntinfo_disc= NULL, bd = bayesian_data, ntinfo = NULL){
  
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
  
  percent <- c()

  alpha_amb <- c(nulling(nulling(dnorms_calc("beta",eval[,"beta"],"alpha", eval2[,"alpha"])) * nulling(alpha[[eval[,"alpha"]]]$total) / nulling(beta[[eval[,"beta"]]]$total)),
                 nulling(nulling(dnorms_calc("gamma",eval[,"gamma"],"alpha", eval2[,"alpha"])) * nulling(alpha[[eval[,"alpha"]]]$total) / nulling(gamma[[eval[,"gamma"]]]$total)),
                 nulling(nulling(dnorms_calc("delta",eval[,"delta"],"alpha", eval2[,"alpha"])) * nulling(alpha[[eval[,"alpha"]]]$total) / nulling(delta[[eval[,"delta"]]]$total)),
                 nulling(nulling(dnorms_calc("epsilon",eval[,"epsilon"],"alpha", eval2[,"alpha"])) * nulling(alpha[[eval[,"alpha"]]]$total) / nulling(epsilon[[eval[,"epsilon"]]]$total)),
                 nulling(nulling(dnorms_calc("zeta",eval[,"zeta"],"alpha", eval2[,"alpha"])) * nulling(alpha[[eval[,"alpha"]]]$total) / nulling(zeta[[eval[,"zeta"]]]$total)),
                 nulling(nulling(dnorms_calc("chi",eval[,"chi"],"alpha", eval2[,"alpha"])) * nulling(alpha[[eval[,"alpha"]]]$total) / nulling(chi[[eval[,"chi"]]]$total)))
  
  beta_amb <- c(nulling(nulling(dnorms_calc("alpha",eval[,"alpha"],"beta", eval2[,"beta"])) * nulling(beta[[eval[,"beta"]]]$total) / nulling(alpha[[eval[,"alpha"]]]$total)),
                nulling(nulling(dnorms_calc("gamma",eval[,"gamma"],"beta", eval2[,"beta"])) * nulling(beta[[eval[,"beta"]]]$total) / nulling(gamma[[eval[,"gamma"]]]$total)),
                nulling(nulling(dnorms_calc("delta",eval[,"delta"],"beta", eval2[,"beta"])) * nulling(beta[[eval[,"beta"]]]$total) / nulling(delta[[eval[,"delta"]]]$total)),
                nulling(nulling(dnorms_calc("epsilon",eval[,"epsilon"],"beta", eval2[,"beta"])) * nulling(beta[[eval[,"beta"]]]$total) / nulling(epsilon[[eval[,"epsilon"]]]$total)),
                nulling(nulling(dnorms_calc("zeta",eval[,"zeta"],"beta", eval2[,"beta"])) * nulling(beta[[eval[,"beta"]]]$total) / nulling(zeta[[eval[,"zeta"]]]$total)),
                nulling(nulling(dnorms_calc("chi",eval[,"chi"],"beta", eval2[,"beta"])) * nulling(beta[[eval[,"beta"]]]$total) / nulling(chi[[eval[,"chi"]]]$total)))
    
  gamma_amb <- c(nulling(nulling(dnorms_calc("alpha",eval[,"alpha"],"gamma", eval2[,"gamma"])) * nulling(gamma[[eval[,"gamma"]]]$total) / nulling(alpha[[eval[,"alpha"]]]$total)),
                 nulling(nulling(dnorms_calc("beta",eval[,"beta"],"gamma", eval2[,"gamma"])) * nulling(gamma[[eval[,"gamma"]]]$total) / nulling(beta[[eval[,"beta"]]]$total)),
                 nulling(nulling(dnorms_calc("delta",eval[,"delta"],"gamma", eval2[,"gamma"])) * nulling(gamma[[eval[,"gamma"]]]$total) / nulling(delta[[eval[,"delta"]]]$total)),
                 nulling(nulling(dnorms_calc("epsilon",eval[,"epsilon"],"gamma", eval2[,"gamma"])) * nulling(gamma[[eval[,"gamma"]]]$total) / nulling(epsilon[[eval[,"epsilon"]]]$total)),
                 nulling(nulling(dnorms_calc("zeta",eval[,"zeta"],"gamma", eval2[,"gamma"])) * nulling(gamma[[eval[,"gamma"]]]$total) / nulling(zeta[[eval[,"zeta"]]]$total)),
                 nulling(nulling(dnorms_calc("chi",eval[,"chi"],"gamma", eval2[,"gamma"])) * nulling(gamma[[eval[,"gamma"]]]$total) / nulling(chi[[eval[,"chi"]]]$total)))
    
  delta_amb <- c(nulling(nulling(dnorms_calc("alpha",eval[,"alpha"],"delta", eval2[,"delta"])) * nulling(delta[[eval[,"delta"]]]$total) / nulling(alpha[[eval[,"alpha"]]]$total)),
                 nulling(nulling(dnorms_calc("beta",eval[,"beta"],"delta", eval2[,"delta"])) * nulling(delta[[eval[,"delta"]]]$total) / nulling(beta[[eval[,"beta"]]]$total)),
                 nulling(nulling(dnorms_calc("gamma",eval[,"gamma"],"delta", eval2[,"delta"])) * nulling(delta[[eval[,"delta"]]]$total) / nulling(gamma[[eval[,"gamma"]]]$total)),
                 nulling(nulling(dnorms_calc("epsilon",eval[,"epsilon"],"delta", eval2[,"delta"])) * nulling(delta[[eval[,"delta"]]]$total) / nulling(epsilon[[eval[,"epsilon"]]]$total)),
                 nulling(nulling(dnorms_calc("zeta",eval[,"zeta"],"delta", eval2[,"delta"])) * nulling(delta[[eval[,"delta"]]]$total) / nulling(zeta[[eval[,"zeta"]]]$total)),
                 nulling(nulling(dnorms_calc("chi",eval[,"chi"],"delta", eval2[,"delta"])) * nulling(delta[[eval[,"delta"]]]$total) / nulling(chi[[eval[,"chi"]]]$total)))
    
  epsilon_amb <- c(nulling(nulling(dnorms_calc("alpha",eval[,"alpha"],"epsilon", eval2[,"epsilon"])) * nulling(epsilon[[eval[,"epsilon"]]]$total) / nulling(alpha[[eval[,"alpha"]]]$total)),
                   nulling(nulling(dnorms_calc("beta",eval[,"beta"],"epsilon", eval2[,"epsilon"])) * nulling(epsilon[[eval[,"epsilon"]]]$total) / nulling(beta[[eval[,"beta"]]]$total)),
                   nulling(nulling(dnorms_calc("gamma",eval[,"gamma"],"epsilon", eval2[,"epsilon"])) * nulling(epsilon[[eval[,"epsilon"]]]$total) / nulling(gamma[[eval[,"gamma"]]]$total)),
                   nulling(nulling(dnorms_calc("delta",eval[,"delta"],"epsilon", eval2[,"epsilon"])) * nulling(epsilon[[eval[,"epsilon"]]]$total) / nulling(delta[[eval[,"delta"]]]$total)),
                   nulling(nulling(dnorms_calc("zeta",eval[,"zeta"],"epsilon", eval2[,"epsilon"])) * nulling(epsilon[[eval[,"epsilon"]]]$total) / nulling(zeta[[eval[,"zeta"]]]$total)),
                   nulling(nulling(dnorms_calc("chi",eval[,"chi"],"epsilon", eval2[,"epsilon"])) * nulling(epsilon[[eval[,"epsilon"]]]$total) / nulling(chi[[eval[,"chi"]]]$total)))
    
  zeta_amb <- c(nulling(nulling(dnorms_calc("alpha",eval[,"alpha"],"zeta", eval2[,"zeta"])) * nulling(zeta[[eval[,"zeta"]]]$total) / nulling(alpha[[eval[,"alpha"]]]$total)),
                nulling(nulling(dnorms_calc("beta",eval[,"beta"],"zeta", eval2[,"zeta"])) * nulling(zeta[[eval[,"zeta"]]]$total) / nulling(beta[[eval[,"beta"]]]$total)),
                nulling(nulling(dnorms_calc("gamma",eval[,"gamma"],"zeta", eval2[,"zeta"])) * nulling(zeta[[eval[,"zeta"]]]$total) / nulling(gamma[[eval[,"gamma"]]]$total)),
                nulling(nulling(dnorms_calc("delta",eval[,"delta"],"zeta", eval2[,"zeta"])) * nulling(zeta[[eval[,"zeta"]]]$total) / nulling(delta[[eval[,"delta"]]]$total)),
                nulling(nulling(dnorms_calc("epsilon",eval[,"epsilon"],"zeta", eval2[,"zeta"])) * nulling(zeta[[eval[,"zeta"]]]$total) / nulling(epsilon[[eval[,"epsilon"]]]$total)),
                nulling(nulling(dnorms_calc("chi",eval[,"chi"],"zeta", eval2[,"zeta"])) * nulling(zeta[[eval[,"zeta"]]]$total) / nulling(chi[[eval[,"chi"]]]$total)))
    
  chi_amb <- c(nulling(nulling(dnorms_calc("alpha",eval[,"alpha"],"chi", eval2[,"chi"])) * nulling(chi[[eval[,"chi"]]]$total) / nulling(alpha[[eval[,"alpha"]]]$total)),
               nulling(nulling(dnorms_calc("beta",eval[,"beta"],"chi", eval2[,"chi"])) * nulling(chi[[eval[,"chi"]]]$total) / nulling(beta[[eval[,"beta"]]]$total)),
               nulling(nulling(dnorms_calc("gamma",eval[,"gamma"],"chi", eval2[,"chi"])) * nulling(chi[[eval[,"chi"]]]$total) / nulling(gamma[[eval[,"gamma"]]]$total)),
               nulling(nulling(dnorms_calc("delta",eval[,"delta"],"chi", eval2[,"chi"])) * nulling(chi[[eval[,"chi"]]]$total) / nulling(delta[[eval[,"delta"]]]$total)),
               nulling(nulling(dnorms_calc("epsilon",eval[,"epsilon"],"chi", eval2[,"chi"])) * nulling(chi[[eval[,"chi"]]]$total) / nulling(epsilon[[eval[,"epsilon"]]]$total)),
               nulling(nulling(dnorms_calc("zeta",eval[,"zeta"],"chi", eval2[,"chi"])) * nulling(chi[[eval[,"chi"]]]$total) / nulling(zeta[[eval[,"zeta"]]]$total)))
  
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
