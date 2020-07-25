library(minipack)
library(veriNA3d)
library(mclust)
library(ggplot2)


#Read the curated ntinfo data
#path <- "/Users/eric1999/Desktop/Gitlab/trinucleotides/NonRedundantTrinucleotides/rmsdpdbs/04.torsionals/"
path <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/new_ntinfo_retrain/"
path2 <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/BACK/lib/"

ntinfo <- read.table(paste(path, "ntinfo.txt", sep=""), header= TRUE, sep=" ", dec=".")

#Extract the identifiers of interest
ntinfo_m <- c()
ntinfo_m$alpha <- ntinfo[, "alpha"]
ntinfo_m$beta <- ntinfo[, "beta"]
ntinfo_m$gamma <- ntinfo[, "gamma"]
ntinfo_m$delta <- ntinfo[, "delta"]
ntinfo_m$epsilon <- ntinfo[, "epsilon"]
ntinfo_m$zeta <- ntinfo[, "zeta"]
ntinfo_m$chi <- ntinfo[, "chi"]
ntinfo_m <- as.data.frame(ntinfo_m)
ntinfo <- ntinfo_m

#Discretize all data into a 1 degree windowing
ntinfo_disc <- minipack::discretize(ntinfo_m, windows = 1)
#save(ntinfo_disc, file = paste(path2, "Data_Discretized.RData", sep=""))
load(paste(path2, "Data_Discretized.RData", sep=""))

#Extract probabilities for each of the angles
alpha <- probability(angle = "alpha", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
beta <- probability(angle = "beta", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
gamma <- probability(angle = "gamma", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
delta <- probability(angle = "delta", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
epsilon <- probability(angle = "epsilon", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
zeta <- probability(angle = "zeta", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
chi <- probability(angle = "chi", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
bayesian_data <- list(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, zeta=zeta, chi=chi)

path3 <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/BACK/"
save(bayesian_data, file = paste(path3, "Probabilities.RData", sep=""))


#Create a dictionary type
tot_ind <- 1:7
tot_names <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi")
tot_names_alpha <- c("beta", "gamma", "delta", "epsilon", "zeta", "chi", NA)
tot_names_beta <- c("alpha", "gamma", "delta", "epsilon", "zeta", "chi", NA)
tot_names_gamma <- c("alpha", "beta", "delta", "epsilon", "zeta", "chi", NA)
tot_names_delta<- c("alpha", "beta", "gamma", "epsilon", "zeta", "chi", NA)
tot_names_epsilon <- c("alpha", "beta", "gamma", "delta", "zeta", "chi", NA)
tot_names_zeta <- c("alpha", "beta", "gamma", "delta", "epsilon", "chi", NA)
tot_names_chi <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", NA)
dict <- data.frame(tot_ind, tot_names, tot_names_alpha, tot_names_beta, tot_names_gamma, tot_names_delta, tot_names_epsilon, tot_names_zeta, tot_names_chi)
naming <- c("index", "total", "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi")
names(dict) <- naming

#######################Use the Scoring function
#Evaluate through the scoring function

evaluation_data <- c()

for(i in 1:length(ntinfo_disc[[1]])){
  tryCatch({
    evaluation_data[i] <- scoring(ntinfo_disc = ntinfo_disc[i,], bd = bayesian_data)
  }, error = function(e) {
    print(paste("Error in", i, sep=""))
  })
}


#Determine max_min to perform normalization
max_v <- max(na.omit(evaluation_data))
min_v <- min(na.omit(evaluation_data))

for(i in 1:length(evaluation_data)){
  if(is.na(evaluation_data[i])){
  }else{
    evaluation_data[i] <- (evaluation_data[i] - min_v) / (max_v - min_v)
  }
}

#Plot
hist(ntinfo$delta, breaks = 100)
par(new=T)
plot(x = ntinfo$delta, y= as.numeric(as.character(evaluation_data)), axes=F)
axis(side=4)


######Scoring using the specific delta function
#######################Use the Scoring function
#Evaluate through the scoring function

evaluation_delta <- c()

for(i in 1:length(ntinfo_disc[[1]])){
  tryCatch({
    evaluation_delta[i] <- scoring_delta(ntinfo_disc = ntinfo_disc[i,], bd = bayesian_data)
  }, error = function(e) {
    print(paste("Error in", i, sep=""))
  })
}


#Determine max_min to perform normalization
max_delta <- max(na.omit(evaluation_delta))
min_delta <- min(na.omit(evaluation_delta))

for(i in 1:length(evaluation_delta)){
  if(is.na(evaluation_delta[i])){
  }else{
    evaluation_delta[i] <- (evaluation_delta[i] - min_delta) / (max_delta - min_delta)
  }
}

#Plot
hist(ntinfo$delta, breaks = 100)
par(new=T)
plot(x = ntinfo$delta, y= as.numeric(as.character(evaluation_delta)), axes=F)
axis(side=4)

#####Perform some data exploration

ind_threshold <- which(evaluation_data > 0.6)
evaluation_data_r <- evaluation_data[-ind_threshold]

#Determine max_min to perform normalization
max_ev <- max(na.omit(evaluation_data_r))
min_ev <- min(na.omit(evaluation_data_r))

for(i in 1:length(evaluation_data_r)){
  if(is.na(evaluation_data_r[i])){
  }else{
    evaluation_data_r[i] <- (evaluation_data_r[i] - min_ev) / (max_ev - min_ev)
  }
}

ntinfo_tt <- ntinfo[-ind_threshold,]

hist(ntinfo_tt$delta, breaks = 100)
par(new=T)
plot(x = ntinfo_tt$delta, y= as.numeric(as.character(evaluation_data_r)), axes=F)
axis(side=4)
