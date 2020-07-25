library(minipack)
library(veriNA3d)
library(mclust)
library(ggplot2)


#Read the curated ntinfo data
path <- "/Users/eric1999/Desktop/Gitlab/trinucleotides/NonRedundantTrinucleotides/rmsdpdbs/04.torsionals/"
path2 <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/BACK/lib/"
lpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/2.rna_puckering/100_analysis/"

ntinfo_m <- read.table(paste(path, "ntinfo.txt", sep=""), header= TRUE, sep=" ", dec=".")

####################Prepare a plot for North and South RNA conformations
#Segment the data between north and south nucleotides

#North angles are comprised between 0 to 36 angles
nor <- 0:36
nor <- as.numeric(nor)

#Make nucleotides rounded
ntinfo_m$phase_round <- round(ntinfo_m$pu_phase,0)

ind <- which(ntinfo_m$phase_round %in% nor)

northern_nt <- ntinfo_m[ind,]
northern <- northern_nt[,32:37]
northern$chi <- northern_nt$chi

#South angles are comprised between 144 to 180 angles

sou <- 144:180
sou <- as.numeric(sou)

ind_sou <- which(ntinfo_m$phase_round %in% sou)

southern_nt <- ntinfo_m[ind_sou,]
southern <- southern_nt[,32:37]
southern$chi <- southern_nt$chi



################################Generate sampling of the data to extract from northern
for(i in 1:100){
  random <- sample (c(1:length(northern$alpha)), size=length(southern$alpha), replace =F)
  
  #Select fixed ntinfo values 
  fix <- 37:143
  fix_2 <- 181:360
  tot_fix <- append(fix, fix_2)
  ind_fix <- which(ntinfo_m$phase_round %in% tot_fix)
  
  #Select the residues
  must <- ntinfo_m[ind_fix,] 
  must_nt <- must[,32:37]
  must_nt$chi <- must$chi
  
  #Select northern sample
  nort <- ntinfo_m[random,]
  nort_nt <- nort[,32:37]
  nort_nt$chi <- nort$chi
  
  sout_nt <- southern
  
  #Bind newly generated dataset
  total <- rbind(nort_nt, must_nt, sout_nt)
  ntinfo_cat <- total
  
  
  #Define the path to save the new models
  
  windows <- 1
  ntinfo_disc <- minipack::discretize(ntinfo_m = ntinfo_cat, windows = windows)
  
  #Non-discretized: ntinfo_cat
  #Discretized: ntinfo_disc
  
  #Check for all different angles
  
  #Extract probabilities for each of the angles
  alpha <- probability(angle = "alpha", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo_cat)
  beta <- probability(angle = "beta", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo_cat)
  gamma <- probability(angle = "gamma", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo_cat)
  delta <- probability(angle = "delta", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo_cat)
  epsilon <- probability(angle = "epsilon", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo_cat)
  zeta <- probability(angle = "zeta", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo_cat)
  chi <- probability(angle = "chi", ntinfo_disc = ntinfo_disc, ntinfo = ntinfo_cat)
  bayesian_data <- list(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, zeta=zeta, chi=chi)
  
  save(bayesian_data, file = paste(lpath, "Reference_1_degree_", i, ".RData", sep=""))
}

