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

vpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/BACK/0.get/Data/"
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

lpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/BACK/2.rna_puckering/100_analysis/"

################################Generate sampling of the data to extract from northern
for(i in 1:100){
  random <- sample (c(1:length(northern_nt$alpha)), size=length(southern_nt$alpha), replace =F)
  
  #Select fixed ntinfo values 
  fix <- 37:143
  fix_2 <- 181:360
  tot_fix <- append(fix, fix_2)
  ind_fix <- which(ntinfo_m$phase_round %in% tot_fix)
  
  #Select the residues
  must_nt <- ntinfo_m[ind_fix,] 
  #must_nt <- must[,66:71]
  #must_nt$chi <- must$chi
  
  #Select northern sample
  nort <- ntinfo_m[random,]
  #nort_nt <- nort[,66:71]
  #nort_nt$chi <- nort$chi
  
  sout_nt <- southern_nt
  
  
  #Bind newly generated dataset
  total <- rbind(nort, must_nt, sout_nt)
  ntinfo_cat <- total
  
  ntinfo_cat <- ntinfo_cat[,1:7]
  #Define the path to save the new models
  
  windows <- 4
  ntinfo_disc <- minipack::discretize(ntinfo_m = ntinfo_cat, windows = windows)

  #Check for all different angles

  alpha <- minipack::probability(angle = "alpha", ntinfo_disc = ntinfo_disc)
  beta <- minipack::probability(angle = "beta", ntinfo_disc = ntinfo_disc)
  gamma <- minipack::probability(angle = "gamma", ntinfo_disc = ntinfo_disc)
  delta <- minipack::probability(angle = "delta", ntinfo_disc = ntinfo_disc)
  epsilon <- minipack::probability(angle = "epsilon", ntinfo_disc = ntinfo_disc)
  zeta <- minipack::probability(angle = "zeta", ntinfo_disc = ntinfo_disc)
  chi <- minipack::probability(angle = "chi", ntinfo_disc = ntinfo_disc)
  bayesian_data <- list(windows=windows, alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, zeta=zeta, chi=chi)

  save(bayesian_data, file = paste(lpath, "Reference_4_degree_", i, ".RData", sep=""))
}



####Loop over all structures and save bayes sample in a new vector

#Load one of the bayesian_data to have a scheme for later on
i <- 1
load(paste(lpath, "Reference_4_degree_", i, ".RData", sep=""))

# Looping
bayesian_total <- vector("list", 100)

for(i in 1:100){
  #Load the bayesian_data
  load(paste(lpath, "Reference_4_degree_", i, ".RData", sep=""))
  
  #Process with removing total info
  bayesian_data$windows <- NULL
  for(m in 1:length(bayesian_data)){
    for(l in 1:length(bayesian_data[[1]])){
      bayesian_data[[m]][[l]]$total <- NULL
    }
  }
  
  bayesian_total[[i]] <- bayesian_data
}


####Make the average for all 100 samples and save the output

bayes_angle <- vector("list", 100)


# All this tructure in order to analyse the bayesian_total[[1]][[1]][[1]]
for(kk in 1:7){
  for(tt in 1:90){
    for(mm in 1:6){
      for(mi in 1:100){
        for(l in 1:length(bayesian_total[[1]][[1]][[1]][[1]])){
          bayes_angle[[l]][[mi]] <- bayesian_total[[mi]][[kk]][[tt]][[mm]][l]
        }
      }
      bayes_angle <- bayes_angle[1:90]
      for(nti in 1:length(bayes_angle)){
        bayes_angle[[nti]] <- mean(as.numeric(bayes_angle[[nti]]))
        #bayesian_data[[1]][[1]][[mm]] <- unlist(bayes_angle)
      }
       bayesian_data[[kk]][[tt]][[mm]] <- unlist(bayes_angle)
    }
  }
}

#Save data with another vector name

bayesian_du <- bayesian_data

##########################Extract information on $TOTAL

# Looping
bayesian_total <- vector("list", 100)

for(i in 1:100){
  #Load the bayesian_data
  load(paste(lpath, "Reference_4_degree_", i, ".RData", sep=""))
  
  #Process with removing total info
  bayesian_data$windows <- NULL
  #for(m in 1:length(bayesian_data)){
  #  for(l in 1:length(bayesian_data[[1]])){
  #    bayesian_data[[m]][[l]]$total <- NULL
  #  }
 # }
  
  bayesian_total[[i]] <- bayesian_data
}


#Extract total angles for each of them
# All this tructure in order to analyse the bayesian_total[[1]][[1]][[1]]
for(kk in 1:7){
  for(tt in 1:90){
      for(mi in 1:100){
          bayes_angle[[mi]] <- bayesian_total[[mi]][[kk]][[tt]]$total
      }
      bayes_angle <- mean(unlist(bayes_angle))
      bayesian_data[[kk]][[tt]]$total <- bayes_angle
  }
}

#Copy data to the previous dataframe

for(kk in 1:7){
  for(tt in 1:90){
    bayesian_du[[kk]][[tt]]$total <- bayesian_data[[kk]][[tt]]$total
  }
}



# Calculate total for every angle

mpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/BACK/2.rna_puckering/"

save(bayesian_du, file = paste(mpath, "Reference_Data_Homog.RData"))
