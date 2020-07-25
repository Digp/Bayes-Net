scoring_all <- function(evaluate, bd = bayesian_data){
  
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

  amb <- c(alpha[[eval[,"alpha"]]]$beta[eval[,"beta"]] * alpha[[eval[,"alpha"]]]$total / beta[[eval[,"beta"]]]$total,
    alpha[[eval[,"alpha"]]]$gamma[eval[,"gamma"]] * alpha[[eval[,"alpha"]]]$total / gamma[[eval[,"gamma"]]]$total,
    alpha[[eval[,"alpha"]]]$delta[eval[,"delta"]] * alpha[[eval[,"alpha"]]]$total / delta[[eval[,"delta"]]]$total,
    alpha[[eval[,"alpha"]]]$epsilon[eval[,"epsilon"]] * alpha[[eval[,"alpha"]]]$total / epsilon[[eval[,"epsilon"]]]$total,
    alpha[[eval[,"alpha"]]]$zeta[eval[,"zeta"]] * alpha[[eval[,"alpha"]]]$total / zeta[[eval[,"zeta"]]]$total,
    alpha[[eval[,"alpha"]]]$chi[eval[,"chi"]] * alpha[[eval[,"alpha"]]]$total / chi[[eval[,"chi"]]]$total,
    alpha[[eval[,"alpha"]]]$eta[eval[,"eta"]] * alpha[[eval[,"alpha"]]]$total / eta[[eval[,"eta"]]]$total,
    alpha[[eval[,"alpha"]]]$theta[eval[,"theta"]] * alpha[[eval[,"alpha"]]]$total / theta[[eval[,"theta"]]]$total,
    
    beta[[eval[,"beta"]]]$alpha[eval[,"alpha"]] * beta[[eval[,"beta"]]]$total / alpha[[eval[,"alpha"]]]$total,
    beta[[eval[,"beta"]]]$gamma[eval[,"gamma"]] * beta[[eval[,"beta"]]]$total / gamma[[eval[,"gamma"]]]$total,
    beta[[eval[,"beta"]]]$delta[eval[,"delta"]] * beta[[eval[,"beta"]]]$total / delta[[eval[,"delta"]]]$total,
    beta[[eval[,"beta"]]]$epsilon[eval[,"epsilon"]] * beta[[eval[,"beta"]]]$total / epsilon[[eval[,"epsilon"]]]$total,
    beta[[eval[,"beta"]]]$zeta[eval[,"zeta"]] * beta[[eval[,"beta"]]]$total / zeta[[eval[,"zeta"]]]$total,
    beta[[eval[,"beta"]]]$chi[eval[,"chi"]] * beta[[eval[,"beta"]]]$total / chi[[eval[,"chi"]]]$total,
    beta[[eval[,"beta"]]]$eta[eval[,"eta"]] * beta[[eval[,"beta"]]]$total / eta[[eval[,"eta"]]]$total,
    beta[[eval[,"beta"]]]$theta[eval[,"theta"]] * beta[[eval[,"beta"]]]$total / theta[[eval[,"theta"]]]$total,
    
    gamma[[eval[,"gamma"]]]$alpha[eval[,"alpha"]] * gamma[[eval[,"gamma"]]]$total / alpha[[eval[,"alpha"]]]$total,
    gamma[[eval[,"gamma"]]]$beta[eval[,"beta"]] * gamma[[eval[,"gamma"]]]$total / beta[[eval[,"beta"]]]$total,
    gamma[[eval[,"gamma"]]]$delta[eval[,"delta"]] * gamma[[eval[,"gamma"]]]$total / delta[[eval[,"delta"]]]$total,
    gamma[[eval[,"gamma"]]]$epsilon[eval[,"epsilon"]] * gamma[[eval[,"gamma"]]]$total / epsilon[[eval[,"epsilon"]]]$total,
    gamma[[eval[,"gamma"]]]$zeta[eval[,"zeta"]] * gamma[[eval[,"gamma"]]]$total / zeta[[eval[,"zeta"]]]$total,
    gamma[[eval[,"gamma"]]]$chi[eval[,"chi"]] * gamma[[eval[,"gamma"]]]$total / chi[[eval[,"chi"]]]$total,
    gamma[[eval[,"gamma"]]]$eta[eval[,"eta"]] * gamma[[eval[,"gamma"]]]$total / eta[[eval[,"eta"]]]$total,
    gamma[[eval[,"gamma"]]]$theta[eval[,"theta"]] * gamma[[eval[,"gamma"]]]$total / theta[[eval[,"theta"]]]$total,
    
    delta[[eval[,"delta"]]]$alpha[eval[,"alpha"]] * delta[[eval[,"delta"]]]$total / alpha[[eval[,"alpha"]]]$total,
    delta[[eval[,"delta"]]]$beta[eval[,"beta"]] * delta[[eval[,"delta"]]]$total / beta[[eval[,"beta"]]]$total,
    delta[[eval[,"delta"]]]$gamma[eval[,"gamma"]] * delta[[eval[,"delta"]]]$total / gamma[[eval[,"gamma"]]]$total,
    delta[[eval[,"delta"]]]$epsilon[eval[,"epsilon"]] * delta[[eval[,"delta"]]]$total / epsilon[[eval[,"epsilon"]]]$total,
    delta[[eval[,"delta"]]]$zeta[eval[,"zeta"]] * delta[[eval[,"delta"]]]$total / zeta[[eval[,"zeta"]]]$total,
    delta[[eval[,"delta"]]]$chi[eval[,"chi"]] * delta[[eval[,"delta"]]]$total / chi[[eval[,"chi"]]]$total,
    delta[[eval[,"delta"]]]$eta[eval[,"eta"]] * delta[[eval[,"delta"]]]$total / eta[[eval[,"eta"]]]$total,
    delta[[eval[,"delta"]]]$theta[eval[,"theta"]] * delta[[eval[,"delta"]]]$total / theta[[eval[,"theta"]]]$total,
    
    epsilon[[eval[,"epsilon"]]]$alpha[eval[,"alpha"]] * epsilon[[eval[,"epsilon"]]]$total / alpha[[eval[,"alpha"]]]$total,
    epsilon[[eval[,"epsilon"]]]$beta[eval[,"beta"]] * epsilon[[eval[,"epsilon"]]]$total / beta[[eval[,"beta"]]]$total,
    epsilon[[eval[,"epsilon"]]]$gamma[eval[,"gamma"]] * epsilon[[eval[,"epsilon"]]]$total / gamma[[eval[,"gamma"]]]$total,
    epsilon[[eval[,"epsilon"]]]$delta[eval[,"delta"]] * epsilon[[eval[,"epsilon"]]]$total / delta[[eval[,"delta"]]]$total,
    epsilon[[eval[,"epsilon"]]]$zeta[eval[,"zeta"]] * epsilon[[eval[,"epsilon"]]]$total / zeta[[eval[,"zeta"]]]$total,
    epsilon[[eval[,"epsilon"]]]$chi[eval[,"chi"]] * epsilon[[eval[,"epsilon"]]]$total / chi[[eval[,"chi"]]]$total,
    epsilon[[eval[,"epsilon"]]]$eta[eval[,"eta"]] * epsilon[[eval[,"epsilon"]]]$total / eta[[eval[,"eta"]]]$total,
    epsilon[[eval[,"epsilon"]]]$theta[eval[,"theta"]] * epsilon[[eval[,"epsilon"]]]$total / theta[[eval[,"theta"]]]$total,
    
    zeta[[eval[,"zeta"]]]$alpha[eval[,"alpha"]] * zeta[[eval[,"zeta"]]]$total / alpha[[eval[,"alpha"]]]$total,
    zeta[[eval[,"zeta"]]]$beta[eval[,"beta"]] * zeta[[eval[,"zeta"]]]$total / beta[[eval[,"beta"]]]$total,
    zeta[[eval[,"zeta"]]]$gamma[eval[,"gamma"]] * zeta[[eval[,"zeta"]]]$total / gamma[[eval[,"gamma"]]]$total,
    zeta[[eval[,"zeta"]]]$delta[eval[,"delta"]] * zeta[[eval[,"zeta"]]]$total / delta[[eval[,"delta"]]]$total,
    zeta[[eval[,"zeta"]]]$epsilon[eval[,"epsilon"]] * zeta[[eval[,"zeta"]]]$total / epsilon[[eval[,"epsilon"]]]$total,
    zeta[[eval[,"zeta"]]]$chi[eval[,"chi"]] * zeta[[eval[,"zeta"]]]$total / chi[[eval[,"chi"]]]$total,
    zeta[[eval[,"zeta"]]]$eta[eval[,"eta"]] * zeta[[eval[,"zeta"]]]$total / eta[[eval[,"eta"]]]$total,
    zeta[[eval[,"zeta"]]]$theta[eval[,"theta"]] * zeta[[eval[,"zeta"]]]$total / theta[[eval[,"theta"]]]$total,
    
    chi[[eval[,"chi"]]]$alpha[eval[,"alpha"]] * chi[[eval[,"chi"]]]$total / alpha[[eval[,"alpha"]]]$total,
    chi[[eval[,"chi"]]]$beta[eval[,"beta"]] * chi[[eval[,"chi"]]]$total / beta[[eval[,"beta"]]]$total,
    chi[[eval[,"chi"]]]$gamma[eval[,"gamma"]] * chi[[eval[,"chi"]]]$total / gamma[[eval[,"gamma"]]]$total,
    chi[[eval[,"chi"]]]$delta[eval[,"delta"]] * chi[[eval[,"chi"]]]$total / delta[[eval[,"delta"]]]$total,
    chi[[eval[,"chi"]]]$epsilon[eval[,"epsilon"]] * chi[[eval[,"chi"]]]$total / epsilon[[eval[,"epsilon"]]]$total,
    chi[[eval[,"chi"]]]$eta[eval[,"eta"]] * chi[[eval[,"chi"]]]$total / eta[[eval[,"eta"]]]$total,
    chi[[eval[,"chi"]]]$theta[eval[,"theta"]] * chi[[eval[,"chi"]]]$total / theta[[eval[,"theta"]]]$total,
    chi[[eval[,"chi"]]]$zeta[eval[,"zeta"]] * chi[[eval[,"chi"]]]$total / zeta[[eval[,"zeta"]]]$total)
  
  return(round(sum(amb)/length(amb), 3))
}
