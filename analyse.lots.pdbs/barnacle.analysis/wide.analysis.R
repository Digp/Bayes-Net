############PERFORM EXPLORATORY ANALYSIS
library(minipack)
library(veriNA3d)
library(mclust)
library(ggplot2)

path_analysis <- "/Users/eric1999/Desktop/Bayessian_RNA_val/3.analyse.pdbs/barnacle.analysis/"

load(paste(path_analysis, "Scoring.RData", sep=""))
load(paste(path_analysis, "PDB_ID.RData", sep=""))


#######################DO NOT RUN. ALREADY GENERATED
scoring_total <- rbind(scoring_complete[[1]], scoring_complete[[2]])

for(i in 3:length(scoring_complete)){
  scoring_total <- rbind(scoring_total, scoring_complete[[i]])
}


scoring_total <- as.data.frame(scoring_total)
data <- cbind(files_complete, scoring_total)
names(data) <- c("PDBs", names(scoring_total))

#save(data, file= paste(path_analysis, "Data_Appended.RData", sep=""))

#################################################
#Find the mean values overall

for(i in 1:length(data$PDBs)){
  data$mean[i] <- mean(as.numeric(data[i, names(scoring_total)]))
}

#Sort in descending way
sorted.values <- sort(data$mean, decreasing = TRUE)

#Find the maximum evaluated pdbs and select the PDBs
maximums <- sorted.values[1:10]

best_evaluate <- data[which(data$mean %in% maximums),]


#Find the worst evaluated pdbs and select the PDBs
minimums <- sorted.values[(length(sorted.values)-10):length(sorted.values)]

worst_evaluate <- data[which(data$mean %in% minimums),]
##########################

#Check the folded conformation pdb
