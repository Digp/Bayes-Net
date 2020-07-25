library(minipack)
library(bio3d)
library(veriNA3d)
library(lattice)
library(rbmn)
library(BBmisc)


#Path to the scoring_cybeRNAting
wpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/Data/"

########RUN FOR ALL DATASET
load(paste(wpath, "ntinfo_id_dssr.RData", sep=""))

#ntinfo_m <- na.omit(ntinfo)
ntinfo_m <- ntinfo
#Create plots of correlations

ntdata <- ntinfo_m[66:71]
ntdata$chi <- ntinfo_m$chi

###########LOAD CYBERTNATING SCORE AND OBTAIN FOR SOUTH PUCKERING
load(paste(wpath, " Scoring_cybeRNAting.RData", sep=""))


#South angles are comprised between 144 to 180 angles
sou <- 144:180
sou <- as.numeric(sou)

#Make nucleotides rounded to 0
ntinfo_m$phase_round <- round(ntinfo_m$pu_phase,0)

ind_sou <- which(ntinfo_m$phase_round %in% sou)

#Extract the nucleotides from the ntinfo dataset
southern_nt <- ntinfo_m[ind_sou,]
sout_nt <- southern_nt[,66:71]
sout_nt$chi <- southern_nt$chi

#Extract south nucleotides from the scoring function

score_south_cybe <- score_cybeRNAting[ind_sou,]
score <- as.numeric(as.character(score_south_cybe$Scoring_CybeRNAting))

#Determine max and min scoring for south angles
min_sout <- min(na.omit(score))
max_sout <- max(na.omit(score))

################DETERMINE AN ANGLE THAT DOES NOT WORK
###Apply a threshold of all structures with score upper than 0.6 scoring

passed_thresh <- which(score > 0.6)

sout_nt[passed_thresh,]
