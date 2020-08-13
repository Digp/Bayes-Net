library(minipack)
library(veriNA3d)
library(mclust)
library(ggplot2)


#Read the curated ntinfo data
#path <- "/Users/eric1999/Desktop/Gitlab/trinucleotides/NonRedundantTrinucleotides/rmsdpdbs/04.torsionals/"
path <-"/orozco/projects/PDB.datamining/trinucleotides/NonRedundantTrinucleotides/rmsdpdbs/04.torsionals/"
path2 <- "/Users/eric1999/Desktop/Bayes-Net/BACK/Mixed_Model_Update/"

ntinfo <- read.table(paste(path2, "ntinfo.txt", sep=""), header= TRUE, sep=" ", dec=".")
load(paste0(path, "ntinfo.RData"))

PDB_ID <- as.character(ntinfo$ntID)

#Extract the identifiers of interest
ntinfo_m <- c()
ntinfo_m$alpha <- ntinfo[, "alpha"]
ntinfo_m$beta <- ntinfo[, "beta"]
ntinfo_m$gamma <- ntinfo[, "gamma"]
ntinfo_m$delta <- ntinfo[, "delta"]
ntinfo_m$epsilon <- ntinfo[, "epsilon"]
ntinfo_m$zeta <- ntinfo[, "zeta"]
ntinfo_m$chi <- ntinfo[, "chi"]

ntinfo_m$alphaplus_one <- ntinfo[, "alphaplus_one"]
ntinfo_m$betaplus_one <- ntinfo[, "betaplus_one"]
ntinfo_m$gammaplus_one <- ntinfo[, "gammaplus_one"]
ntinfo_m$deltaplus_one <- ntinfo[, "deltaplus_one"]
ntinfo_m$epsilonplus_one <- ntinfo[, "epsilonplus_one"]
ntinfo_m$zetaplus_one <- ntinfo[, "zetaplus_one"]
ntinfo_m$chiplus_one <- ntinfo[, "chiplus_one"]

ntinfo_m$alphaminus_one <- ntinfo[, "alphaminus_one"]
ntinfo_m$betaminus_one <- ntinfo[, "betaminus_one"]
ntinfo_m$gammaminus_one <- ntinfo[, "gammaminus_one"]
ntinfo_m$deltaminus_one <- ntinfo[, "deltaminus_one"]
ntinfo_m$epsilonminus_one <- ntinfo[, "epsilonminus_one"]
ntinfo_m$zetaminus_one <- ntinfo[, "zetaminus_one"]
ntinfo_m$chiminus_one <- ntinfo[, "chiminus_one"]

ntinfo_m <- as.data.frame(ntinfo_m)
ntinfo <- ntinfo_m

names(ntinfo) <- c( "alpha","beta", "gamma","delta","epsilon","zeta","chi",             
                    "alpha_plus",    "beta_plus",     "gamma_plus",    "delta_plus",    "epsilon_plus",  "zeta_plus",     "chi_plus",     
                    "alpha_minus",   "beta_minus",    "gamma_minus",   "delta_minus" ,  "epsilon_minus", "zeta_minus",    "chi_minus") 

#save(ntinfo, file = paste0(path2, "ntinfo.RData"))
#save(PDB_ID, file = paste0(path2,"PDB_ID.RData"))

load( file = paste0(path2, "ntinfo.RData"))
load(file = paste0(path2, "PDB_ID.RData"))

#Remove the following columns: alpha_minus, epsilon_plus and zeta_plus
ntinfo[, "alpha_minus"] <- NULL
ntinfo[, "epsilon_plus"] <- NULL
ntinfo[, "zeta_plus"] <- NULL




#Discretize all data into a 1 degree windowing
#ntinfo_disc <- minipack::discretize(ntinfo_m, windows = 1)
ntinfo_disc <- round(ntinfo, 0)
#save(ntinfo_disc, file = paste(path2, "Data_Discretized.RData", sep=""))
#load(paste(path2, "Data_Discretized.RData", sep=""))

#Extract probabilities for each of the angles
alpha <- probability(angle = "alpha", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
beta <- probability(angle = "beta", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
gamma <- probability(angle = "gamma", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
delta <- probability(angle = "delta", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
epsilon <- probability(angle = "epsilon", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
zeta <- probability(angle = "zeta", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
chi <- probability(angle = "chi", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
alpha_plus <- probability(angle = "alpha_plus", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
beta_plus <- probability(angle = "beta_plus", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
gamma_plus <- probability(angle = "gamma_plus", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
delta_plus <- probability(angle = "delta_plus", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
chi_plus <- probability(angle = "chi_plus", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
beta_minus <- probability(angle = "beta_minus", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
gamma_minus <- probability(angle = "gamma_minus", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
delta_minus <- probability(angle = "delta_minus", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
epsilon_minus <- probability(angle = "epsilon_minus", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
zeta_minus <- probability(angle = "zeta_minus", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)
chi_minus <- probability(angle = "chi_minus", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo)

bayesian_data <- list(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, zeta=zeta, chi=chi,
                      alpha_plus=alpha_plus, beta_plus=beta_plus, gamma_plus=gamma_plus, delta_plus=delta_plus,
                      chi_plus=chi_plus, beta_minus=beta_minus, gamma_minus=gamma_minus, delta_minus=delta_minus,
                      epsilon_minus=epsilon_minus, zeta_minus=zeta_minus, chi_minus=chi_minus)

path3 <- "/Users/eric1999/Desktop/Bayes-Net/BACK/Mixed_Model_Update/"
save(bayesian_data, file = paste(path3, "Probabilities.RData", sep=""))


#Create a dictionary type
tot_ind <- 1:18
tot_names <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi", "alpha_plus", "beta_plus", "gamma_plus", "delta_plus",
               "chi_plus", "beta_minus", "gamma_minus", "delta_minus", "epsilon_minus", "zeta_minus", "chi_minus" )
tot_names_alpha <- c("beta", "gamma", "delta", "epsilon", "zeta", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_beta <- c("alpha", "gamma", "delta", "epsilon", "zeta", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_gamma <- c("alpha", "beta", "delta", "epsilon", "zeta", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_delta<- c("alpha", "beta", "gamma", "epsilon", "zeta", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_epsilon <- c("alpha", "beta", "gamma", "delta", "zeta", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_zeta <- c("alpha", "beta", "gamma", "delta", "epsilon", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_chi <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_alpha_plus <- c("beta", "gamma", "delta", "epsilon", "zeta", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_beta_plus <- c("alpha", "gamma", "delta", "epsilon", "zeta", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_gamma_plus <- c("alpha", "beta", "delta", "epsilon", "zeta", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_delta_plus<- c("alpha", "beta", "gamma", "epsilon", "zeta", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_chi_plus <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_beta_minus <- c("alpha", "gamma", "delta", "epsilon", "zeta", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_gamma_minus <- c("alpha", "beta", "delta", "epsilon", "zeta", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_delta_minus<- c("alpha", "beta", "gamma", "epsilon", "zeta", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_epsilon_minus <- c("alpha", "beta", "gamma", "delta", "zeta", "chi", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_zeta_minus <- c("alpha", "beta", "gamma", "delta", "epsilon", "chi",NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tot_names_chi_minus <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta",NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
dict <- data.frame(tot_ind, tot_names, tot_names_alpha, tot_names_beta, tot_names_gamma, tot_names_delta, tot_names_epsilon, tot_names_zeta, tot_names_chi,
                   tot_names_alpha_plus, tot_names_beta_plus, tot_names_gamma_plus,tot_names_delta_plus, tot_names_chi_plus, tot_names_beta_minus,
                   tot_names_gamma_minus, tot_names_delta_minus, tot_names_epsilon_minus, tot_names_zeta_minus, tot_names_chi_minus)
naming <- c("index", "total", "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi", "alpha_plus", "beta_plus", "gamma_plus", "delta_plus",
            "chi_plus", "beta_minus", "gamma_minus", "delta_minus", "epsilon_minus", "zeta_minus", "chi_minus" )
names(dict) <- naming

save(dict, file = paste(path3, "Dictionary.RData", sep=""))
#######################Use the Scoring function
#Evaluate through the scoring function

#####Perform scoring

evaluation_data <- scoring_bayes(ntinfo_disc = ntinfo_disc[1,], bd = bayesian_data, ntinfo = ntinfo[1,])
evaluation_data <- as.data.frame(evaluation_data)

for(i in 2:length(ntinfo[[1]])){
  amb <- scoring_bayes(ntinfo_disc = ntinfo_disc[i,], bd = bayesian_data, ntinfo = ntinfo[i,])
  amb <- as.data.frame(amb)
  evaluation_data <- rbind(evaluation_data, amb)
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
hist(ntinfo$alpha, breaks = 100)
par(new=T)
plot(x = ntinfo$alpha, y= as.numeric(as.character(evaluation_data$alpha)), axes=F)
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
