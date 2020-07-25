scoring <- function(evaluate, bd = bayesian_data){
  
  windows <- bd$windows
  eta <- bd$eta
  theta <- bd$theta
  theta_mas_1 <- bd$theta_mas_1
  eta_mas_1 <- bd$eta_mas_1
  eta_menos_1 <- bd$eta_menos_1
  theta_menos_1 <- bd$theta_menos_1
  
  eval <- minipack::discretize(ntinfo_m = evaluate, windows = windows)
  
  ind <-which(eval %in% 0)
  if(length(ind) > 0){
    eval[ind] <- 1
  }

  amb <- c(eta[[eval[,"eta"]]]$theta[eval[,"theta"]] * eta[[eval[,"eta"]]]$total / theta[[eval[,"theta"]]]$total,
    eta[[eval[,"eta"]]]$eta_menos_1[eval[,"eta_menos_1"]] * eta[[eval[,"eta"]]]$total / eta_menos_1[[eval[,"eta_menos_1"]]]$total,
    eta[[eval[,"eta"]]]$theta_menos_1[eval[,"theta_menos_1"]] * eta[[eval[,"eta"]]]$total / theta_menos_1[[eval[,"theta_menos_1"]]]$total,
    eta[[eval[,"eta"]]]$eta_mas_1[eval[,"eta_mas_1"]] * eta[[eval[,"eta"]]]$total / eta_mas_1[[eval[,"eta_mas_1"]]]$total,

    
    theta[[eval[,"theta"]]]$eta[eval[,"eta"]] * theta[[eval[,"theta"]]]$total / eta[[eval[,"eta"]]]$total,
    theta[[eval[,"theta"]]]$theta_mas_1[eval[,"theta_mas_1"]] * theta[[eval[,"theta"]]]$total / theta_mas_1[[eval[,"theta_mas_1"]]]$total,
    theta[[eval[,"theta"]]]$theta_menos_1[eval[,"theta_menos_1"]] * theta[[eval[,"theta"]]]$total / theta_menos_1[[eval[,"theta_menos_1"]]]$total,
    theta[[eval[,"theta"]]]$eta_mas_1[eval[,"eta_mas_1"]] * theta[[eval[,"theta"]]]$total / eta_mas_1[[eval[,"eta_mas_1"]]]$total)
  
  return(round(sum(amb)/length(amb), 3))
}
