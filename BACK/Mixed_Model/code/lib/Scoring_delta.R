scoring_delta <- function(ntinfo_disc= NULL, bd = bayesian_data, ntinfo = NULL){
  
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
  
  ind <-which(eval %in% 0)
  if(length(ind) > 0){
    eval[ind] <- 1
  }
  
  amb <- c(
           nulling(dnorms_calc("delta",eval[,"delta"],"alpha") * delta[[eval[,"delta"]]]$total / alpha[[eval[,"alpha"]]]$total),
           nulling(dnorms_calc("delta",eval[,"delta"],"beta") * delta[[eval[,"delta"]]]$total / beta[[eval[,"beta"]]]$total),
           nulling(dnorms_calc("delta",eval[,"delta"],"gamma") * delta[[eval[,"delta"]]]$total / gamma[[eval[,"gamma"]]]$total),
           nulling(dnorms_calc("delta",eval[,"delta"],"epsilon") * delta[[eval[,"delta"]]]$total / epsilon[[eval[,"epsilon"]]]$total),
           nulling(dnorms_calc("delta",eval[,"delta"],"zeta") * delta[[eval[,"delta"]]]$total / zeta[[eval[,"zeta"]]]$total),
           nulling(dnorms_calc("delta",eval[,"delta"],"chi") * delta[[eval[,"delta"]]]$total / chi[[eval[,"chi"]]]$total))

  
  return(sum(amb)/length(amb))
}
