############PERFORM EXPLORATORY ANALYSIS
library(minipack)
library(veriNA3d)
library(mclust)
library(ggplot2)

######LOAD THE DATA
#Read the curated ntinfo data
path <- "/Users/eric1999/Desktop/Gitlab/trinucleotides/NonRedundantTrinucleotides/rmsdpdbs/04.torsionals/"
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

#Obtian discretized data
load(paste(path2, "Data_Discretized.RData", sep=""))

#Load the dict option
load(file = "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/BACK/lib/Dictionary.RData")
load(file = "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/BACK/Probabilities_5.RData")
bd <- bayesian_data

#####Perform scoring

evaluation_data <- scoring_bayes(ntinfo_disc = ntinfo_disc[1,], bd = bayesian_data, ntinfo = ntinfo[1,])
evaluation_data <- as.data.frame(evaluation_data)

for(i in 2:length(ntinfo[[1]])){
  amb <- scoring_bayes(ntinfo_disc = ntinfo_disc[i,], bd = bayesian_data, ntinfo = ntinfo[i,])
  amb <- as.data.frame(amb)
  evaluation_data <- rbind(evaluation_data, amb)
}

#Plot
hist(ntinfo$delta, breaks = 100)
par(new=T)
plot(x = ntinfo$delta, y= as.numeric(as.character(evaluation_data$delta)), axes=F)
axis(side=4)

#####Perform scoring

scoring_theorem <- scoring(ntinfo_disc = ntinfo_disc[1,], bd = bayesian_data, ntinfo = ntinfo[1,])
scoring_theorem <- as.data.frame(evaluation_data)

for(i in 2:length(ntinfo[[1]])){
  amb <- scoring(ntinfo_disc = ntinfo_disc[i,], bd = bayesian_data, ntinfo = ntinfo[i,])
  amb <- as.data.frame(amb)
  scoring_theorem <- rbind(scoring_theorem, amb)
}

#Plot
hist(ntinfo$chi, breaks = 100)
par(new=T)
plot(x = ntinfo$chi, y= as.numeric(as.character(scoring_theorem$chi)), axes=F)
axis(side=4)
