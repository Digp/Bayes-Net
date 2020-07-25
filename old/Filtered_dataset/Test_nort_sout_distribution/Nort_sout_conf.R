library(minipack)
library(veriNA3d)
library(bio3d)

###################################
#This code is used in order to analyse the distribution of a nort and south conf nucleotides when
#ranging static values for all angles except from delta ranging from 1-360
##################################3


path <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/BACK/0.get/Data/"
mpath <- "/Users/eric1999/Desktop/Gitlab/trinucleotides/Mining/rmsdpdbs/04.torsionals/"

load(paste(mpath, "PDBID_id_dssr.RData", sep=""))

########RUN FOR ALL DATASET
load(paste(vpath, "ntinfo_filtered.RData"))
######Extract id_dssr information unique identifiers
ntinfo_l <- ntinfo
################################
#Retreive important nucleotide angles
ntinfo_m <- ntinfo_l[, 32:37]
ntinfo_m$chi <- ntinfo_l$chi
ntinfo_m$phase_round <- ntinfo_l$pu_phase
#ntinfo_l[is.na(ntinfo_l)] <- 0
#ntinfo_m <- na.omit(ntinfo_m)
ntinfo_m$alpha <- as.numeric(as.character(ntinfo_m$alpha))
ntinfo_m$beta <- as.numeric(as.character(ntinfo_m$beta))
ntinfo_m$gamma <- as.numeric(as.character(ntinfo_m$gamma))
ntinfo_m$delta <- as.numeric(as.character(ntinfo_m$delta))
ntinfo_m$epsilon <- as.numeric(as.character(ntinfo_m$epsilon))
ntinfo_m$zeta <- as.numeric(as.character(ntinfo_m$zeta))
ntinfo_m$chi <- as.numeric(as.character(ntinfo_m$chi))
ntinfo_m$phase_round <- as.numeric(as.character(ntinfo_m$phase_round))
ntinfo_m$PDB_ID_DSSR <- PDBID_id_dssr



##################Prepare a plot for all dihedrals of RNA
#Plot correlations for each of the pseudotorsional angles pairs
#p <- ggpairs(southern_nt[1:7])+ theme_bw()
#p

####################Prepare a plot for North and South RNA conformations
#Segment the data between north and south nucleotides

#North angles are comprised between 0 to 36 angles
nor <- 0:36
nor <- as.numeric(nor)

#Make nucleotides rounded
ntinfo_m$phase_round <- round(ntinfo_m$phase_round,0)

ind <- which(ntinfo_m$phase_round %in% nor)

northern_nt <- ntinfo_m[ind,]
#northern <- northern_nt[66:71]
#northern$chi <- northern_nt$chi

#South angles are comprised between 144 to 180 angles

sou <- 144:180
sou <- as.numeric(sou)

ind_sou <- which(ntinfo_m$phase_round %in% sou)

southern_nt <- ntinfo_m[ind_sou,]
#southern <- southern_nt[,66:71]
#southern$chi <- southern_nt$chi


#Select one of the nucleotides to be analyzed
random <- sample (c(1:length(northern_nt$alpha)), size=1, replace =F)
nort_select <- northern_nt[random,]
nort_sel <- nort_select[1:7]

nort_sel$delta <- 1
nort_data <- nort_sel

for(i in 2:360){
  nort_sel$delta <- i
  nort_data <- rbind(nort_data, nort_sel)
}


######EVALUATE SCORING WITH THE SCORING() FUNCTION

data_all <- c()

for(i in 1:length(sout_data$alpha)){
  data_all[i] <- scoring_wochi(evaluate = sout_data[i,], bd = bayesian_du)
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

################ANALYSIS
#Bind different data
bound <- cbind(sout_data, data_normalized)
bound <- as.data.frame(bound)

#Retreive info for nort evaluation
bound_n <- bound[nor,]
bound_s <- bound[sou,]



#Select one of the nucleotides to be analyzed
random <- sample (c(1:length(southern_nt$alpha)), size=1, replace =F)
sout_select <- southern_nt[random,]
sout_sel <- sout_select[1:7]

sout_sel$delta <- 1
sout_data <- sout_sel

for(i in 2:360){
  sout_sel$delta <- i
  sout_data <- rbind(sout_data, sout_sel)
}


################PLOT
#Plot_Score_v1
hist(sout_data$delta, breaks = 100)
par(new=T)
plot(x = sout_data$delta, y= as.numeric(as.character(data_all)), axes=F)
axis(side=4)
dev.off()

