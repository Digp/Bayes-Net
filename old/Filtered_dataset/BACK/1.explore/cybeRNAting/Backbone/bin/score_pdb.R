#!/usr/local/bin/Rscript
## Diego Gallego
## Date: 2018-Feb-26
## This script takes a pdb file and computes the scoring funciton on it, then
## the output is saved in the pdb format and different informative .txt files
   
## Check correct usage

#CLIargs <- c("1a9n", "/Users/eric1999/Desktop/cybeRNAting/Scripts/Backbone/dat/referencedata.txt", "scoring_function.R", 0.2, "F")


## Load dependencies ...
library(bio3d)
library(veriNA3d)
library(magrittr)
library(BBmisc)
## ... and input data
pdb           <- read.pdb(CLIargs[1])
referencedata <- read.table(CLIargs[2], header=T, stringsAsFactors=F)
source(CLIargs[3])
threshold     <- CLIargs[4]
fixedname     <- CLIargs[5]

vpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/BACK/0.get/Data/"
mpath <- "/Users/eric1999/Desktop/Gitlab/trinucleotides/Mining/rmsdpdbs/04.torsionals/"
load(paste(mpath, "PDBID_id_dssr.RData", sep=""))

#pdbID <- "1a9n"
#pdb <- pipeNucData(toupper(pdbID))
referencedata <- read.table("/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/BACK/1.explore/cybeRNAting/Backbone/dat/referencedata.txt", header=T, stringsAsFactors=F)
#torsionals <- pdb

## Solve matching problem ## Temporal PATCH
referencedata$angle <- gsub("v", "nu", referencedata$angle)

## Save pdb name if necessary
#if (fixedname) {
#    out_prefix <- "out"
#} else {
#    out_prefix <- strsplit(CLIargs[1], "\\.pdb")[[1]]
#}

## Get pdb data 
#torsionals <- measureNuc(pdb, v_shifted=FALSE, distances=NA, angles=NA, pucker=F, Dp=F)


#Load data on the whole dataset that we have from all the RNA structures
load(paste(vpath, "ntinfo_filtered.RData"))

#Retreive important nucleotide angles
ntinfo_m <- ntinfo

torsionals <- ntinfo_m
torsionals$ntID <- 1:length(torsionals[[1]])
torsionals$resid <- as.character(ntinfo_m$resid)
torsionals$resno <- as.character(ntinfo_m$resno)
torsionals$model <- as.character(ntinfo_m$model)
torsionals$chain <- as.character(ntinfo_m$chain)
torsionals$alpha <- as.numeric(as.character(ntinfo_m$alpha))
torsionals$beta <- as.numeric(as.character(ntinfo_m$beta))
torsionals$gamma <- as.numeric(as.character(ntinfo_m$gamma))
torsionals$delta <- as.numeric(as.character(ntinfo_m$delta))
torsionals$epsilon <- as.numeric(as.character(ntinfo_m$epsilon))
torsionals$zeta <- as.numeric(as.character(ntinfo_m$zeta))
torsionals$chi <- as.numeric(as.character(ntinfo_m$chi))
torsionals$nu0 <- as.numeric(as.character(ntinfo_m$nu0))
torsionals$nu1 <- as.numeric(as.character(ntinfo_m$nu1))
torsionals$nu2 <- as.numeric(as.character(ntinfo_m$nu2))
torsionals$nu3 <- as.numeric(as.character(ntinfo_m$nu3))
torsionals$nu4 <- as.numeric(as.character(ntinfo_m$nu4))
torsionals$eta <- as.numeric(as.character(ntinfo_m$eta))
torsionals$theta <- as.numeric(as.character(ntinfo_m$theta))



#Ids <- as.character(droplevels.factor(ntinfo$ID_DSSR))

## Compute scores and save summary

scores <- scoring_function(torsionals, referencedata)

#Retreive data of interest
backbone_scores <- score_backbone(scores)

ind <- which(is.nan(backbone_scores))
data_calc <- backbone_scores[-ind]

#Standarize obtained values
min_cybernating <- min(data_calc)
max_cybernating <- max(data_calc)
min_max_cybeRNAting <- cbind(min_cybernating, max_cybernating)

#Normalize data
normalized_cybernating <- c()

for(i in 1:length(backbone_scores)){
  normalized_cybernating[i] <- (backbone_scores[i] - min_cybernating) / (max_cybernating - min_cybernating)
}

bulk_cybernating <- cbind(PDBID_id_dssr, backbone_scores, normalized_cybernating)
bulk_cybernating <- as.data.frame(bulk_cybernating)
namy <- c("PDB_ID_DSSR", "Scoring_cybernating", "Normalized_cybernating")
names(bulk_cybernating) <- namy

save(bulk_cybernating, file = paste(vpath, "Scoring_cybeRNAting.RData"))


s