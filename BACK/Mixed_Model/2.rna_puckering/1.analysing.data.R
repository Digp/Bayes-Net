

################################################################################################
################################################################################################

#Load all 100 generated files
bayesian_total <- c()

for(i in 1:100){
  load(paste(lpath, "Reference_1_degree_", i, ".RData", sep=""))
  bayesian_total[[i]] <- bayesian_data
}

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


evaluation_sample <- scoring_sampling(ntinfo_disc = ntinfo_disc[1,], bang = bayesian_total, ntinfo = ntinfo[1,])
evaluation_sample <- as.data.frame(evaluation_sample)

for(m in 2:length(ntinfo[[1]])){
  tryCatch({
    amb <- scoring_sampling(ntinfo_disc = ntinfo_disc[m,], bang = bayesian_total, ntinfo = ntinfo[m,])
    amb <- as.data.frame(amb)
    evaluation_sample <- rbind(evaluation_sample, amb)
  }, error = function(e) {
    #amb<- data.frame(NA, NA, NA, NA, NA, NA, NA, NA)
    print(i)
  })
}


#Selected non-useful nucleotides
ind_chi <- which(evaluation_sample$chi > 5)
ind_beta <- which(evaluation_sample$beta > 5)

ind_tot <- unlist(list(ind_chi, ind_beta))

#Remove selected residues from the dataset
ntinfo <- ntinfo[-ind_tot,]

ntinfo_disc <- minipack::discretize(ntinfo, windows = 1)

#Plot
hist(ntinfo$gamma, breaks = 100)
par(new=T)
plot(x = ntinfo$gamma, y= as.numeric(as.character(evaluation_sample$gamma)), axes=F)
axis(side=4)


##########################################
#Let's work for the gamma angle

i_1 <- 0:110
i_2 <- 111:225
i_3 <- 226:360

nt_1 <- which(ntinfo_disc$gamma %in% i_1)
nt_2 <- which(ntinfo_disc$gamma %in% i_2)
nt_3 <- which(ntinfo_disc$gamma %in% i_3)

sc_1 <- evaluation_sample$gamma[nt_1]

sc_2 <- evaluation_sample$gamma[nt_2]
sc_3 <- evaluation_sample$gamma[nt_3]

norm_sc2 <- normalize(sc_2)
norm_sc3 <- normalize(sc_3)

evaluation_sample$gamma[nt_1] <- sc_1
evaluation_sample$gamma[nt_2] <- norm_sc2
evaluation_sample$gamma[nt_3] <- norm_sc3
