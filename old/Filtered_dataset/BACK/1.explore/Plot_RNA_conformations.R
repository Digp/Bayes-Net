###Load proper dependencies
library(minipack)
library(bio3d)
library(veriNA3d)
library(lattice)
library(rbmn)
library(BBmisc)
library(signal)
library(bnlearn)
library(scatterplot3d)
library(GGally)
library(ggplot2)
library("PerformanceAnalytics")
library("factoextra")

wpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/Data/"

########RUN FOR ALL DATASET
load(paste(wpath, "ntinfo_id_dssr.RData", sep=""))

ntinfo_m <- na.omit(ntinfo)

#Create plots of correlations

ntdata <- ntinfo_m[66:71]
ntdata$chi <- ntinfo_m$chi


##################Prepare a plot for all dihedrals of RNA
#Plot correlations for each of the pseudotorsional angles pairs
p <- ggpairs(ntdata)+ theme_bw()
p

####################Prepare a plot for North and South RNA conformations
#Segment the data between north and south nucleotides

#North angles are comprised between 0 to 36 angles
nor <- 0:36
nor <- as.numeric(nor)

#Make nucleotides rounded
ntinfo_m$phase_round <- round(ntinfo_m$pu_phase,0)

ind <- which(ntinfo_m$phase_round %in% nor)

northern_nt <- ntinfo_m[ind,]
northern <- northern_nt[66:71]
northern$chi <- northern_nt$chi

#South angles are comprised between 144 to 180 angles

sou <- 144:180
sou <- as.numeric(sou)

ind_sou <- which(ntinfo_m$phase_round %in% sou)

southern_nt <- ntinfo_m[ind_sou,]
southern <- southern_nt[66:71]
southern$chi <- southern_nt$chi

#Make the plot
jpeg(file=paste("/Users/eric1999/Desktop/Bayessian_RNA_val/Plot_Northern_Conf.jpg", sep=""))
nort_plot <- ggpairs(northern)+ theme_bw()
nort_plot
dev.off()

jpeg(file=paste("/Users/eric1999/Desktop/Bayessian_RNA_val/Plot_Southern_Conf.jpg", sep=""))
sout_plot <- ggpairs(southern)+ theme_bw()
sout_plot
dev.off()


#Plot for all different pseudo-torsionals in South puckering
#jpeg(file=paste("/Users/eric1999/Desktop/Bayessian_RNA_val/Plot_South_All_angles.jpg", sep=""))
attach(southern)
par(mfrow=c(3,2))
hist(southern$alpha, breaks = 100)
hist(southern$beta, breaks = 100)
hist(southern$gamma, breaks = 100)
hist(southern$delta, breaks = 100)
hist(southern$epsilon, breaks = 100)
hist(southern$zeta, breaks = 100)
#dev.off()

#Plot for all different pseudo-torsionals in North puckering

attach(northern)
par(mfrow=c(3,2))
hist(northern$alpha, breaks = 100)
hist(northern$beta, breaks = 100)
hist(northern$gamma, breaks = 100)
hist(northern$delta, breaks = 100)
hist(northern$epsilon, breaks = 100)
hist(northern$zeta, breaks = 100)
