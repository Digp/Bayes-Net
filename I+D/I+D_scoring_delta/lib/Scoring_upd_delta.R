##Newly updated scoring function

scoring_up_delta <- function(evaluate, bd = bayesian_data){
  
  windows <- bd$windows
  alpha <- bd$alpha
  beta <- bd$beta
  gamma <- bd$gamma
  delta <- bd$delta
  epsilon <- bd$epsilon
  zeta <- bd$zeta
  
  eval <- discretize(ntinfo_m = evaluate, windows = windows)
  
  ind <-which(eval %in% 0)
  if(length(ind) > 0){
    eval[ind] <- 1
  }
  
  amb <- c(
           delta[[eval[,"delta"]]]$alpha[eval[,"alpha"]] * delta[[eval[,"delta"]]]$total / alpha[[eval[,"alpha"]]]$total,
           delta[[eval[,"delta"]]]$beta[eval[,"beta"]] * delta[[eval[,"delta"]]]$total / beta[[eval[,"beta"]]]$total,
           delta[[eval[,"delta"]]]$gamma[eval[,"gamma"]] * delta[[eval[,"delta"]]]$total / gamma[[eval[,"gamma"]]]$total,
           delta[[eval[,"delta"]]]$epsilon[eval[,"epsilon"]] * delta[[eval[,"delta"]]]$total / epsilon[[eval[,"epsilon"]]]$total,
           delta[[eval[,"delta"]]]$zeta[eval[,"zeta"]] * delta[[eval[,"delta"]]]$total / zeta[[eval[,"zeta"]]]$total)
  
 
    
  
  #Process for delta
  
  #Sum
  if(is.na(eval[,"delta"])){
    amt <- amb
  }else if(eval[,"delta"] +1 >36){
    amt <- append(amb, c(delta[[1]]$alpha[eval[,"alpha"]] * delta[[1]]$total / alpha[[eval[,"alpha"]]]$total,
    delta[[1]]$beta[eval[,"beta"]] * delta[[1]]$total / beta[[eval[,"beta"]]]$total,
    delta[[1]]$gamma[eval[,"gamma"]] * delta[[1]]$total / gamma[[eval[,"gamma"]]]$total,
    delta[[1]]$epsilon[eval[,"epsilon"]] * delta[[1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    delta[[1]]$zeta[eval[,"zeta"]] * delta[[1]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"delta"]+1 <= 36){
    amt <- append(amb, c(delta[[eval[,"delta"]+1]]$alpha[eval[,"alpha"]] * delta[[eval[,"delta"]+1]]$total / alpha[[eval[,"alpha"]]]$total,
    delta[[eval[,"delta"]+1]]$beta[eval[,"beta"]] * delta[[eval[,"delta"]+1]]$total / beta[[eval[,"beta"]]]$total,
    delta[[eval[,"delta"]+1]]$gamma[eval[,"gamma"]] * delta[[eval[,"delta"]+1]]$total / gamma[[eval[,"gamma"]]]$total,
    delta[[eval[,"delta"]+1]]$epsilon[eval[,"epsilon"]] * delta[[eval[,"delta"]+1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    delta[[eval[,"delta"]+1]]$zeta[eval[,"zeta"]] * delta[[eval[,"delta"]+1]]$total / zeta[[eval[,"zeta"]]]$total))
  }
  
  #Substract
  
  if(is.na(eval[,"delta"])){
    amt <- amb
  }else if(eval[,"delta"]-1 < 1){
    amt <- append(amb, c(delta[[36]]$alpha[eval[,"alpha"]] * delta[[36]]$total / alpha[[eval[,"alpha"]]]$total,
    delta[[36]]$beta[eval[,"beta"]] * delta[[36]]$total / beta[[eval[,"beta"]]]$total,
    delta[[36]]$gamma[eval[,"gamma"]] * delta[[36]]$total / gamma[[eval[,"gamma"]]]$total,
    delta[[36]]$epsilon[eval[,"epsilon"]] * delta[[36]]$total / epsilon[[eval[,"epsilon"]]]$total,
    delta[[36]]$zeta[eval[,"zeta"]] * delta[[36]]$total / zeta[[eval[,"zeta"]]]$total))
  }else if(eval[,"delta"]-1 >=1){
    amt <- append(amb,  c(delta[[eval[,"delta"]-1]]$alpha[eval[,"alpha"]] * delta[[eval[,"delta"]-1]]$total / alpha[[eval[,"alpha"]]]$total,
    delta[[eval[,"delta"]-1]]$beta[eval[,"beta"]] * delta[[eval[,"delta"]-1]]$total / beta[[eval[,"beta"]]]$total,
    delta[[eval[,"delta"]-1]]$gamma[eval[,"gamma"]] * delta[[eval[,"delta"]-1]]$total / gamma[[eval[,"gamma"]]]$total,
    delta[[eval[,"delta"]-1]]$epsilon[eval[,"epsilon"]] * delta[[eval[,"delta"]-1]]$total / epsilon[[eval[,"epsilon"]]]$total,
    delta[[eval[,"delta"]-1]]$zeta[eval[,"zeta"]] * delta[[eval[,"delta"]-1]]$total / zeta[[eval[,"zeta"]]]$total))
  }
  
  return(round(sum(amt)/length(amt), 3))
}
