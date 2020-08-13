############PERFORM EXPLORATORY ANALYSIS
library(minipack)
library(veriNA3d)
library(mclust)
library(ggplot2)

######LOAD THE DATA
#Read the curated ntinfo data
path <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/new_ntinfo_retrain/"
path2 <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/BACK/lib/"

load(paste(path2, "Dictionary.RData", sep=""))

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

PDBID_id_dssr <- as.character(ntinfo$ntID)


#Obtian discretized data
#load(paste(path2, "Data_Discretized.RData", sep=""))

#Perform scoring with obtaining the data for each single nucleotide as data_frame
evaluation_data <- scoring_bayes(ntinfo_disc = ntinfo_disc[1,], bd = bayesian_data, ntinfo = ntinfo[1,])
evaluation_data <- as.data.frame(evaluation_data)

for(i in 2:length(ntinfo_disc[[1]])){
  tryCatch({
    amb <- scoring_bayes(ntinfo_disc = ntinfo_disc[i,], bd = bayesian_data, ntinfo = ntinfo[i,])
    amb <- as.data.frame(amb)
    evaluation_data <- rbind(evaluation_data, amb)
  }, error = function(e) {
    #amb<- data.frame(NA, NA, NA, NA, NA, NA, NA, NA)
    print(i)
  })
}

#Plot
hist(ntinfo$zeta, breaks = 100)
par(new=T)
plot(x = ntinfo$zeta, y= as.numeric(as.character(evaluation_data$zeta)), axes=F)
axis(side=4)

####################################################################
#Find outliers based on visualization

##########
#Alpha_outlier
ind_a <- which(ntinfo$alpha > 340)
alpha_outlier<- PDB_ID[ind_a]
alpha_out <- paste("alpha", alpha_outlier, sep=":")

##########
#Gamma_outlier
ind_g <- which(ntinfo$gamma > 245  & ntinfo$gamma < 260)
gamma_outlier<- PDB_ID[ind_g]
gamma_out <- paste("gamma", gamma_outlier, sep=":")

##########
#Delta_outlier
ind_d <- which(ntinfo$delta > 120 & ntinfo$delta < 124)
delta_outlier<- PDB_ID[ind_d]
delta_out <- paste("delta", delta_outlier, sep=":")

###Another identifier outlier in delta. When analyzing it is found because of NAs
ind_d <- which(ntinfo$delta > 80 & ntinfo$delta < 83)
ntinfo_subset <- evaluation_data[ind_d,]
ind_double <- which(ntinfo_subset$delta < 0.2)
ntinfo_subset[ind_double,]

##########
#Epsilon_outlier
ind_e <- which(ntinfo$epsilon > 330)
epsilon_outlier <- PDB_ID[ind_e]
epsilon_out <- paste("epsilon", epsilon_outlier, sep=":")


#Bind all different outliers and save
outliers <- c(alpha_out, gamma_out, delta_out, epsilon_out)


save(outliers, file = paste(path2, "Outliers.RData", sep=""))



