
#Function to compute the dnorms
loopout <- function(ind = NULL, x = NULL,  eval2 = NULL){
  calcs <- c()
  for(i in 1:ind){
    calcs[i] <-  (dnorm(x = eval2, mean =  x$mean[i], sd = x$sd[i]) * x$prop[i]) / x$max_dens[i]
    #print(calcs)
    #print(calcs[i])
  }
  calcs <- sum(calcs)
  return(calcs)
}

#Function to check over angle1 and angle 2 as well as the angle value
dnorms_calc <- function(angle = NULL, eval = NULL, angle2 = NULL, eval2 = NULL){
  if(is.na(eval) || is.na(eval2)){
    dnorming <- 0
  }else{
    ind_1 <- which(dict$total %in% angle)
    #let <- which(names(dict) %in% angle)
    ind_2 <- which(names(bd[[ind_1]][eval][[1]]) %in% angle2)
    
    #Extract data
    ind <- length(bd[[ind_1]][eval][[1]][[ind_2]]$prop)
    x <- bd[[ind_1]][eval][[1]][[ind_2]]
    dnorming <- loopout(ind = ind, x = x, eval2 = eval2 )
  }
  return(dnorming)
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

#Create a new function to compute all dnorms for a specific file

neotest <- function(angle1 = NULL, eval = NULL, angle2=NULL, eval2 = NULL, bayesian_total = NULL){
  storing_values <- c()
  for(i in 1:100){
    bd <- bayesian_total[[i]]
    calc <- nulling(dnorms_calc(angle1, eval,angle2, eval2))
    storing_values[i] <- calc
  }
  
  #print(storing_values)
  store <- mean(storing_values)
  return(store)
}


