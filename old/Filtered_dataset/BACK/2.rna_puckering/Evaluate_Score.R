###Load proper dependencies
library(minipack)
library(bio3d)
library(veriNA3d)
library(lattice)
library(rbmn)
library(BBmisc)


# Load reference data

mpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/BACK/2.rna_puckering/"
vpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/BACK/0.get/Data/"

load(paste(mpath, "Reference_Data_Homog.RData"))

#Read the desired files

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


bayesian_du$windows <- 4

######EVALUATE SCORING WITH THE SCORING() FUNCTION

data_all <- c()

for(i in 1:length(ntinfo_m$alpha)){
  data_all[i] <- scoring_delta(evaluate = ntinfo_m[i,], bd = bayesian_du)
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

################PLOT
#Plot_Score_v1
jpeg(file=paste(mpath, "/Images/Plot_Subsampling_Delta.jpg", sep=""))
hist(ntinfo_m$delta, breaks = 100)
par(new=T)
plot(x = ntinfo_m$delta, y= as.numeric(as.character(data_all)), axes=F)
axis(side=4)
dev.off()

#Plot_Score_v1_normalized
jpeg(file=paste(mpath, "/Images/Plot_Subsampling_Delta_Normalized.jpg", sep=""))
hist(ntinfo_m$delta, breaks = 100)
par(new=T)
plot(x = ntinfo_m$delta, y= as.numeric(as.character(data_normalized)), axes=F)
axis(side=4)
dev.off()
