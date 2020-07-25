###Load proper dependencies
library(minipack)
library(bio3d)
library(veriNA3d)
library(lattice)
library(rbmn)
library(BBmisc)

#Read the desired files
#path <- "/var/folders/3l/xhm3x7fs02g9ydpvmf9_h5l40000gn/T//RtmpEs3WJo/"
vpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/BACK/0.get/Data/"
mpath <- "/Users/eric1999/Desktop/Gitlab/trinucleotides/Mining/rmsdpdbs/04.torsionals/"

#wpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/Data/"
#cpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/1.explore/"
#spath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/1.explore/"

load(paste(vpath, "Reference_4_degree.RData", sep=""))           
load(paste(mpath, "PDBID_id_dssr.RData", sep=""))
#####TEST WITH A SINGLE PDB FILE

#ntinfo <- read_file(path, pdbID = "1a9n")

#evaluate <- ntinfo[,1:7]

#Scoring function over all data
#data <- c()

#for(i in 1:length(evaluate$alpha)){
#  data[i] <- scoring(evaluate = evaluate[i,], bd = bayesian_data)
#}

#Bind with the nucleotides

#bulk_score <- cbind(ntinfo$id_dssr, data)
#name <- c("Sequence", "ID_DSSR", "Scoring")
#names(bulk_score) <- name


########RUN FOR ALL DATASET
load(paste(vpath, "ntinfo_filtered.RData"))
######Extract id_dssr information unique identifiers
ntinfo_l <- ntinfo
################################
#Retreive important nucleotide angles
ntinfo_m <- ntinfo_l[, 32:37]
ntinfo_m$chi <- ntinfo_l$chi
#ntinfo_l[is.na(ntinfo_l)] <- 0
#ntinfo_m <- na.omit(ntinfo_m)
ntinfo_m$alpha <- as.numeric(as.character(ntinfo_m$alpha))
ntinfo_m$beta <- as.numeric(as.character(ntinfo_m$beta))
ntinfo_m$gamma <- as.numeric(as.character(ntinfo_m$gamma))
ntinfo_m$delta <- as.numeric(as.character(ntinfo_m$delta))
ntinfo_m$epsilon <- as.numeric(as.character(ntinfo_m$epsilon))
ntinfo_m$zeta <- as.numeric(as.character(ntinfo_m$zeta))
ntinfo_m$chi <- as.numeric(as.character(ntinfo_m$chi))


######EVALUATE SCORING WITH THE SCORING() FUNCTION

data_all <- c()

for(i in 1:length(ntinfo_m$alpha)){
  data_all[i] <- scoring_delta(evaluate = ntinfo_m[i,], bd = bayesian_data)
}

data_all <- round(data_all, 3)


#Forget about NaN to compute min and max values
ind <- which(is.nan(data_all))
data_calc <- data_all[-ind]

#Extract the min and the max values
min_score <- min(data_calc)
max_score <- max(data_calc)

min_max_scoring <- cbind(min_score, max_score)

#Normalize values
data_normalized <- c()

for(i in 1:length(data_all)){
  data_normalized[i] <- (data_all[i] - min_score) / (max_score - min_score)
}

data_normalized <- round(data_normalized, 3)

#Bind ID_DSSR with the score

bulk_score_all <- cbind(PDBID_id_dssr, data_all, data_normalized)


#Frequencies
#hist(data_all, breaks = 100)

####EVALUATING SCORE WITH THE SCORING_UPDATED() FUNCTION

data_score_upd <- c()

for(i in 1:length(ntinfo_m$alpha)){
  data_score_upd[i] <- scoring_up_delta(evaluate = ntinfo_m[i,], bd = bayesian_data)
}

data_score_upd <- round(data_score_upd, 3)

#Forget about NaN to compute min and max values
ind <- which(is.nan(data_score_upd))
data_calc_upd <- data_score_upd[-ind]

#Extract the min and the max values
min_score_upd <- min(data_calc_upd)
max_score_upd <- max(data_calc_upd)

min_max_scoring_update <- cbind(min_score_upd, max_score_upd)

#Normalize values
data_normal_upd<- c()

for(i in 1:length(data_all)){
  data_normal_upd[i] <- (data_calc_upd[i] - min_score_upd) / (max_score_upd - min_score_upd)
}

data_normal_upd <- round(data_normal_upd, 3)

bulk_score_upd <- cbind(PDBID_id_dssr, data_score_upd, data_normal_upd)


#Load results from backend cybeRNAting
load(paste(vpath, "Scoring_cybeRNAting.RData"))


#Bind ID_DSSR with the score

Data_analysis_delta <- cbind(PDBID_id_dssr, data_all, data_normalized, data_score_upd, data_normal_upd, bulk_cybernating$Scoring_cybernating,
                       bulk_cybernating$Normalized_cybernating)

Data_analysis_delta <- as.data.frame(Data_analysis_delta)
namy <- c("PDB_ID_DSSR", "Scoring", "Scoring_Normalized", "Scoring_Updated", "Normalized_Scor_Updated", "cybeRNAting_score", "Normalized_cybeRNAting")
names(Data_analysis_delta) <- namy


#Save the data
save(Data_analysis_delta, file = paste(vpath, "Scoring_Analysis_Delta.RData"))


###############MAKE THE PLOTTING

mpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/BACK/1.explore/Plots/Delta_versus_delta/"

#Plot_Score_v1_normalized
jpeg(file=paste(mpath, "Plot_Normal_v1_Delta.jpg", sep=""))
hist(ntinfo_m$delta, breaks = 100)
par(new=T)
plot(x = ntinfo_m$delta, y= data_normalized, axes=F)
axis(side=4)
dev.off()

#Plot_Score_v2_normalized
jpeg(file=paste(mpath, "Plot_Normal_v2_Delta.jpg", sep=""))
hist(ntinfo_m$delta, breaks = 100)
par(new=T)
plot(x = ntinfo_m$delta, y= data_normal_upd, axes=F)
axis(side=4)
dev.off()

#Plot_Score_cybeRNAting_normalized
jpeg(file=paste(mpath, "Plot_cybeRNAting_Delta.jpg", sep=""))
hist(ntinfo_m$delta, breaks = 100)
par(new=T)
plot(x = ntinfo_m$delta, y= bulk_cybernating$Normalized_cybernating, axes=F)
axis(side=4)
dev.off()
