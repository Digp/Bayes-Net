###Load proper dependencies
library(minipack)
library(bio3d)
library(veriNA3d)
library(lattice)
library(rbmn)
library(BBmisc)

impath_2 <- "/Users/eric1999/Desktop/Bayessian_RNA_val/dihedra_eta/"
#Folder of v1 code
#impath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/dihedra_eta/Evaluate_dihedra_eta_v1/Data/"
#Folder of v2 code
impath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/dihedra_eta/Evaluate_dihedra_eta_v2/Data/"
impath3 <- "/Users/eric1999/Desktop/Bayessian_RNA_val/dihedra_eta/Evaluate_dihedra_eta_v2/"

source(paste(impath3, "Scoring_dihedral_v2.R", sep=""))
source(paste(impath3, "probability_seg.R", sep=""))

load(paste(impath_2, "ntinfo_dihedrals_eta.RData", sep=""))
load(paste(impath, "Bayesian_data_4_degree_v2.RData", sep=""))

######EVALUATE SCORING WITH THE SCORING() FUNCTION

data_all <- c()

for(i in 1:length(ntinfo_m$eta)){
  data_all[i] <- scoring_seg(evaluate = ntinfo_l[i,], bd = bayesian_data)
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

#Plot_Score_Delta
#jpeg(file=paste(impath,"Eta_Scoring.jpg", sep=""))
#hist(ntinfo_l$delta, breaks = 100)
#par(new=T)
#plot(x = ntinfo_l$eta, y= as.numeric(as.character(data_all)), axes=F)
#axis(side=4)
#dev.off()