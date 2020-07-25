
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

#Perform scoring with obtaining the data for each single nucleotide as data_frame
evaluation_data <- scoring(ntinfo_disc = ntinfo_disc[1,], bd = bayesian_data, ntinfo = ntinfo[1,])
evaluation_data <- as.data.frame(evaluation_data)
  
for(i in 2:length(ntinfo_disc[[1]])){
  tryCatch({
    amb <- scoring(ntinfo_disc = ntinfo_disc[i,], bd = bayesian_data, ntinfo = ntinfo[i,])
    amb <- as.data.frame(amb)
    evaluation_data <- rbind(evaluation_data, amb)
  }, error = function(e) {
    #amb<- data.frame(NA, NA, NA, NA, NA, NA, NA, NA)
    print(i)
  })
}

save(evaluation_data, file = paste(path2, "Scored_values.RData", sep=""))

#Perform normalization
norm_data <- c()
norm_data$alpha <- normalize(data = evaluation_data$alpha)
norm_data$beta <- normalize(data = evaluation_data$beta)
norm_data$gamma <- normalize(data = evaluation_data$gamma)
norm_data$delta <- normalize(data = evaluation_data$delta)
norm_data$epsilon <- normalize(data = evaluation_data$epsilon)
norm_data$zeta <- normalize(data = evaluation_data$zeta)
norm_data$chi <- normalize(data = evaluation_data$chi)
norm_data$total <- normalize(data = evaluation_data$total)
norm_data <- as.data.frame(norm_data)

################################################################################################
################################################################################################
#Explore for a north configuration nucleotide
ind_north <- which(ntinfo_disc$delta %in% 70:90)
data_north <- ntinfo[ind_north, ]

#Extract the scoring of those ntds
scor_north <- evaluation_data[ind_north,]

#Select max scored value on all of them 
mas <- max(scor_north$delta)
ntd_max_north <- which(scor_north$delta == mas)   ###Nucleotide 2175 has the higher score

#Max nucleotides
max_north <- scor_north[ntd_max_north,] #3205 from original dataset
max_ntinfo <- ntinfo[3205,]

#Select another random nucleotide from the range. We picked the 89
rand_north <- scoring(ntinfo_disc = ntinfo_disc[89,], bd = bayesian_data)
rand_ntinfo <- ntinfo[89,]
rand_north <- as.data.frame(rand_north)

################################################################################################
################################################################################################


###############Perform analysis from 1:360 on delta angle
nt_data <- max_ntinfo
nt_data <- minipack::discretize(nt_data, windows = 1)
amb <- nt_data

#Generate new dataset
for(i in 1:360){
  amb$delta <- i
  nt_data <- rbind(nt_data, amb)
}

nt_data <- nt_data[-1,]



#Evaluate scoring
nt_all <- scoring(ntinfo_disc = nt_data[1,], bd = bayesian_data)
nt_all <- as.data.frame(nt_all)

for(i in 2:length(nt_data[[1]])){
    amb <- scoring(ntinfo_disc = nt_data[i,], bd = bayesian_data)
    amb <- as.data.frame(amb)
    nt_all <- rbind(nt_all, amb)
}
  

ggplot(nt_all, aes(x = 1:360, y = nt_all$delta)) + geom_point()

################################################################################################
################################################################################################

#Evaluate scoring with the scoring_bayes function without applying bayes theorem

#Perform scoring with obtaining the data for each single nucleotide as data_frame
eval_bayes <- scoring_bayes(ntinfo_disc = ntinfo_disc[1,], bd = bayesian_data, ntinfo = ntinfo[1,])
eval_bayes <- as.data.frame(eval_bayes)

for(i in 2:length(ntinfo_disc[[1]])){
  tryCatch({
    amb <- scoring_bayes(ntinfo_disc = ntinfo_disc[i,], bd = bayesian_data, ntinfo = ntinfo[i,])
    amb <- as.data.frame(amb)
    eval_bayes <- rbind(eval_bayes, amb)
  }, error = function(e) {
    #amb<- data.frame(NA, NA, NA, NA, NA, NA, NA, NA)
    print(i)
  })
}

#Perform normalization
norm_bayes <- c()
norm_bayes$alpha <- normalize(data = eval_bayes$alpha)
norm_bayes$beta <- normalize(data = eval_bayes$beta)
norm_bayes$gamma <- normalize(data = eval_bayes$gamma)
norm_bayes$delta <- normalize(data = eval_bayes$delta)
norm_bayes$epsilon <- normalize(data = eval_bayes$epsilon)
norm_bayes$zeta <- normalize(data = eval_bayes$zeta)
norm_bayes$chi <- normalize(data = eval_bayes$chi)
norm_bayes$total <- normalize(data = eval_bayes$total)
norm_bayes <- as.data.frame(norm_bayes)


#Plot
hist(ntinfo$alpha, breaks = 100)
par(new=T)
plot(x = ntinfo$alpha, y= as.numeric(as.character(norm_bayes$alpja)), axes=F)
axis(side=4)



