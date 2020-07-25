###Load proper dependencies
library(minipack)
library(bio3d)
library(veriNA3d)
library(lattice)
library(rbmn)
library(BBmisc)

impath_2 <- "/Users/eric1999/Desktop/Bayessian_RNA_val/eta_theta/"
impath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/eta_theta/Data/"

source(paste(impath_2, "Scoring_function.R", sep=""))
source(paste(impath_2, "probability.R", sep=""))

load(paste(impath, "ntinfo_eta_theta.RData", sep=""))


######EVALUATE SCORING WITH THE SCORING() FUNCTION

data_all <- c()

for(i in 1:length(ntinfo_m$eta)){
  data_all[i] <- scoring(evaluate = data_select[i,], bd = bayesian_data)
}

data_all <- round(data_all, 3)


#Forget about NaN to compute min and max values
ind <- which(is.nan(data_all))
data_calc <- data_all

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

#Plot_Score_Eta
jpeg(file=paste(impath,"Eta_Scoring.jpg", sep=""))
hist(data_select$eta, breaks = 100)
par(new=T)
plot(x = data_select$eta, y= as.numeric(as.character(data_all)), axes=F)
axis(side=4)
dev.off()

#Plot_Score_Eta_Normalized
jpeg(file=paste(impath,"Eta_Scoring_Normalized.jpg", sep=""))
hist(data_select$eta, breaks = 100)
par(new=T)
plot(x = data_select$eta, y= as.numeric(as.character(data_normalized)), axes=F)
axis(side=4)
dev.off()

#Plot_Score_Theta
jpeg(file=paste(impath,"Theta_Scoring.jpg", sep=""))
hist(data_select$theta, breaks = 100)
par(new=T)
plot(x = data_select$theta, y= as.numeric(as.character(data_all)), axes=F)
axis(side=4)
dev.off()

#Plot_Score_Theta_Normalized
jpeg(file=paste(impath,"Theta_Scoring_Normalized.jpg", sep=""))
hist(data_select$theta, breaks = 100)
par(new=T)
plot(x = data_select$theta, y= as.numeric(as.character(data_normalized)), axes=F)
axis(side=4)
dev.off()


