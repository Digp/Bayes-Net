scoring_delta <- function(evaluate, bd = bayesian_data){
  
  windows <- bd$windows
  alpha <- bd$alpha
  beta <- bd$beta
  gamma <- bd$gamma
  delta <- bd$delta
  epsilon <- bd$epsilon
  zeta <- bd$zeta
  
  eval <- minipack::discretize(ntinfo_m = evaluate, windows = windows)
  
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
    
  
  return(round(sum(amb)/length(amb), 3))
}
