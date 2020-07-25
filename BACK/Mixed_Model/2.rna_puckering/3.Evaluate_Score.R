###Load proper dependencies
library(minipack)
library(bio3d)
library(veriNA3d)
library(lattice)
library(rbmn)
library(BBmisc)


# Load reference data

mpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/2.rna_puckering/"

load(paste(mpath, "Reference_Data_Homog.RData"))

#Read the desired files
wpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/Data/"

########RUN FOR ALL DATASET
load(paste(wpath, "ntinfo_id_dssr.RData", sep=""))

#Retreive important nucleotide angles
ntinfo_m <- ntinfo[, 66:71]
ntinfo_m$chi <- ntinfo$chi
Ids <- as.character(droplevels.factor(ntinfo$ID_DSSR))


bayesian_du$windows <- 4

######EVALUATE SCORING WITH THE SCORING() FUNCTION

data <- scoring(evaluate = ntinfo_m[1,], bd = bayesian_du)
data <- as.data.frame(data)

for(i in 2:length(ntinfo_m$alpha)){
  amb <- scoring(evaluate = ntinfo_m[i,], bd = bayesian_du)
  amb <- as.data.frame(amb)
  data <- rbind(data, amb)
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
