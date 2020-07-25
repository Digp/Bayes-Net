scoring <- function(evaluate, bd = bayesian_data){
  
  windows <- bd$windows
  alpha <- bd$alpha
  beta <- bd$beta
  gamma <- bd$gamma
  delta <- bd$delta
  epsilon <- bd$epsilon
  zeta <- bd$zeta
  chi <- bd$chi
  
  eval <- discretize(ntinfo_m = evaluate, windows = windows)
  
  ind <-which(eval %in% 0)
  if(length(ind) > 0){
    eval[ind] <- 1
  }
  
  percent <- c()

  alpha_amb <- c(nulling(alpha[[eval[,"alpha"]]]$beta[eval[,"beta"]] * alpha[[eval[,"alpha"]]]$total / beta[[eval[,"beta"]]]$total),
                 nulling(alpha[[eval[,"alpha"]]]$gamma[eval[,"gamma"]] * alpha[[eval[,"alpha"]]]$total / gamma[[eval[,"gamma"]]]$total),
                 nulling(alpha[[eval[,"alpha"]]]$delta[eval[,"delta"]] * alpha[[eval[,"alpha"]]]$total / delta[[eval[,"delta"]]]$total),
                 nulling(alpha[[eval[,"alpha"]]]$epsilon[eval[,"epsilon"]] * alpha[[eval[,"alpha"]]]$total / epsilon[[eval[,"epsilon"]]]$total),
                 nulling(alpha[[eval[,"alpha"]]]$zeta[eval[,"zeta"]] * alpha[[eval[,"alpha"]]]$total / zeta[[eval[,"zeta"]]]$total),
                 nulling(alpha[[eval[,"alpha"]]]$chi[eval[,"chi"]] * alpha[[eval[,"alpha"]]]$total / chi[[eval[,"chi"]]]$total))
    
  beta_amb <- c(nulling(beta[[eval[,"beta"]]]$alpha[eval[,"alpha"]] * beta[[eval[,"beta"]]]$total / alpha[[eval[,"alpha"]]]$total),
                nulling(beta[[eval[,"beta"]]]$gamma[eval[,"gamma"]] * beta[[eval[,"beta"]]]$total / gamma[[eval[,"gamma"]]]$total),
                nulling(beta[[eval[,"beta"]]]$delta[eval[,"delta"]] * beta[[eval[,"beta"]]]$total / delta[[eval[,"delta"]]]$total),
                nulling(beta[[eval[,"beta"]]]$epsilon[eval[,"epsilon"]] * beta[[eval[,"beta"]]]$total / epsilon[[eval[,"epsilon"]]]$total),
                nulling(beta[[eval[,"beta"]]]$zeta[eval[,"zeta"]] * beta[[eval[,"beta"]]]$total / zeta[[eval[,"zeta"]]]$total),
                nulling(beta[[eval[,"beta"]]]$chi[eval[,"chi"]] * beta[[eval[,"beta"]]]$total / chi[[eval[,"chi"]]]$total))
    
  gamma_amb <- c(nulling(gamma[[eval[,"gamma"]]]$alpha[eval[,"alpha"]] * gamma[[eval[,"gamma"]]]$total / alpha[[eval[,"alpha"]]]$total),
                 nulling(gamma[[eval[,"gamma"]]]$beta[eval[,"beta"]] * gamma[[eval[,"gamma"]]]$total / beta[[eval[,"beta"]]]$total),
                 nulling(gamma[[eval[,"gamma"]]]$delta[eval[,"delta"]] * gamma[[eval[,"gamma"]]]$total / delta[[eval[,"delta"]]]$total),
                 nulling(gamma[[eval[,"gamma"]]]$epsilon[eval[,"epsilon"]] * gamma[[eval[,"gamma"]]]$total / epsilon[[eval[,"epsilon"]]]$total),
                 nulling(gamma[[eval[,"gamma"]]]$zeta[eval[,"zeta"]] * gamma[[eval[,"gamma"]]]$total / zeta[[eval[,"zeta"]]]$total),
                 nulling( gamma[[eval[,"gamma"]]]$chi[eval[,"chi"]] * gamma[[eval[,"gamma"]]]$total / chi[[eval[,"chi"]]]$total))
    
  delta_amb <- c(nulling(delta[[eval[,"delta"]]]$alpha[eval[,"alpha"]] * delta[[eval[,"delta"]]]$total / alpha[[eval[,"alpha"]]]$total),
                 nulling(delta[[eval[,"delta"]]]$beta[eval[,"beta"]] * delta[[eval[,"delta"]]]$total / beta[[eval[,"beta"]]]$total),
                 nulling(delta[[eval[,"delta"]]]$gamma[eval[,"gamma"]] * delta[[eval[,"delta"]]]$total / gamma[[eval[,"gamma"]]]$total),
                 nulling(delta[[eval[,"delta"]]]$epsilon[eval[,"epsilon"]] * delta[[eval[,"delta"]]]$total / epsilon[[eval[,"epsilon"]]]$total),
                 nulling(delta[[eval[,"delta"]]]$zeta[eval[,"zeta"]] * delta[[eval[,"delta"]]]$total / zeta[[eval[,"zeta"]]]$total),
                 nulling(delta[[eval[,"delta"]]]$chi[eval[,"chi"]] * delta[[eval[,"delta"]]]$total / chi[[eval[,"chi"]]]$total))
    
  epsilon_amb <- c(nulling(epsilon[[eval[,"epsilon"]]]$alpha[eval[,"alpha"]] * epsilon[[eval[,"epsilon"]]]$total / alpha[[eval[,"alpha"]]]$total),
                   nulling(epsilon[[eval[,"epsilon"]]]$beta[eval[,"beta"]] * epsilon[[eval[,"epsilon"]]]$total / beta[[eval[,"beta"]]]$total),
                   nulling(epsilon[[eval[,"epsilon"]]]$gamma[eval[,"gamma"]] * epsilon[[eval[,"epsilon"]]]$total / gamma[[eval[,"gamma"]]]$total),
                   nulling(epsilon[[eval[,"epsilon"]]]$delta[eval[,"delta"]] * epsilon[[eval[,"epsilon"]]]$total / delta[[eval[,"delta"]]]$total),
                   nulling(epsilon[[eval[,"epsilon"]]]$zeta[eval[,"zeta"]] * epsilon[[eval[,"epsilon"]]]$total / zeta[[eval[,"zeta"]]]$total),
                   nulling(epsilon[[eval[,"epsilon"]]]$chi[eval[,"chi"]] * epsilon[[eval[,"epsilon"]]]$total / chi[[eval[,"chi"]]]$total))
    
  zeta_amb <- c(nulling(zeta[[eval[,"zeta"]]]$alpha[eval[,"alpha"]] * zeta[[eval[,"zeta"]]]$total / alpha[[eval[,"alpha"]]]$total),
                nulling(zeta[[eval[,"zeta"]]]$beta[eval[,"beta"]] * zeta[[eval[,"zeta"]]]$total / beta[[eval[,"beta"]]]$total),
                nulling(zeta[[eval[,"zeta"]]]$gamma[eval[,"gamma"]] * zeta[[eval[,"zeta"]]]$total / gamma[[eval[,"gamma"]]]$total),
                nulling(zeta[[eval[,"zeta"]]]$delta[eval[,"delta"]] * zeta[[eval[,"zeta"]]]$total / delta[[eval[,"delta"]]]$total),
                nulling(zeta[[eval[,"zeta"]]]$epsilon[eval[,"epsilon"]] * zeta[[eval[,"zeta"]]]$total / epsilon[[eval[,"epsilon"]]]$total),
                nulling(zeta[[eval[,"zeta"]]]$chi[eval[,"chi"]] * zeta[[eval[,"zeta"]]]$total / chi[[eval[,"chi"]]]$total))
    
  chi_amb <- c(nulling(chi[[eval[,"chi"]]]$alpha[eval[,"alpha"]] * chi[[eval[,"chi"]]]$total / alpha[[eval[,"alpha"]]]$total),
               nulling( chi[[eval[,"chi"]]]$beta[eval[,"beta"]] * chi[[eval[,"chi"]]]$total / beta[[eval[,"beta"]]]$total),
               nulling(chi[[eval[,"chi"]]]$gamma[eval[,"gamma"]] * chi[[eval[,"chi"]]]$total / gamma[[eval[,"gamma"]]]$total),
               nulling(chi[[eval[,"chi"]]]$delta[eval[,"delta"]] * chi[[eval[,"chi"]]]$total / delta[[eval[,"delta"]]]$total),
               nulling(chi[[eval[,"chi"]]]$epsilon[eval[,"epsilon"]] * chi[[eval[,"chi"]]]$total / epsilon[[eval[,"epsilon"]]]$total),
               nulling(chi[[eval[,"chi"]]]$zeta[eval[,"zeta"]] * chi[[eval[,"chi"]]]$total / zeta[[eval[,"zeta"]]]$total))
  
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


#Function to print 0 when we have Inf or NaN or NA
nulling <- function(x){
  if(length(x) == 0){
    x <- 0
  }else if(length(x) == 1){
    if (is.nan(x)){
      x <- 0
    }else if (is.na(x)){
      x <- 0
    }else if (x == Inf){
      x <- 0
    }
  }
  return(x)
}

#Normalization function
normalize <- function(data = NULL){
  max_v <- max(na.omit(data))
  min_v <- min(na.omit(data))
  
  data_norm <- c()
  for(i in 1:length(data)){
    data_norm[i] <- (data[i] - min_v) / (max_v - min_v)
  }
  
  return(data_norm)
}
