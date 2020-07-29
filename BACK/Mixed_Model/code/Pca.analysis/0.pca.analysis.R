############PERFORM EXPLORATORY ANALYSIS
library(minipack)
library(veriNA3d)
library(mclust)
library(ggplot2)
library(devtools)
library(ggplot2)
library(ggfortify)
library(cluster)
library(factoextra)

######LOAD THE DATA
#Read the curated ntinfo data
path <- "/Users/eric1999/Desktop/Gitlab/trinucleotides/NonRedundantTrinucleotides/rmsdpdbs/04.torsionals/"
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

#Remove NAs as PCA won't process them
ntinf <- na.omit(ntinfo)

#Perform PCA
pca_res <- prcomp(ntinf)

#Graph
autoplot(pca_res)
autoplot(pca_res, data = ntinf, colour = 'alpha')
autoplot(pca_res, data = ntinf, colour = 'beta')

#Clustering
autoplot(clara(ntinf, 3))
autoplot(fanny(ntinf, 3), frame = TRUE)


res.pca <- prcomp(ntinf)
