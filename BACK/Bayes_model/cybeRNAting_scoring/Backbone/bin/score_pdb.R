#!/usr/local/bin/Rscript
## Diego Gallego
## Date: 2018-Feb-26
## This script takes a pdb file and computes the scoring funciton on it, then
## the output is saved in the pdb format and different informative .txt files
   
## Check correct usage

path <- "/Users/eric1999/Desktop/Bayessian/"

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


#pdbID <- "1a9n"
#pdb <- pipeNucData(toupper(pdbID))
referencedata <- read.table("/Users/eric1999/Desktop/cybeRNAting/Scripts/Backbone/dat/referencedata.txt", header=T, stringsAsFactors=F)
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
load("/Users/eric1999/Desktop/Bayessian/ntinfo_id_dssr.RData")

#Retreive important nucleotide angles
ntinfo_m <- ntinfo

torsionals <- ntinfo_m

Ids <- as.character(droplevels.factor(ntinfo$ID_DSSR))

## Compute scores and save summary

scores <- scoring_function(torsionals, referencedata)

#Retreive data of interest
backbone_scores <- score_backbone(scores)

#Standarize obtained values
stand <- normalize(as.numeric(bulk_score_all$Scoring), method = "standardize")

bulk_score_all <- cbind(Ids, backbone_scores, stand)
bulk_score_all <- as.data.frame(bulk_score_all)
namy <- c("ID_DSSR", "Scoring", "Normalized Val")
names(bulk_score_all) <- namy

save(bulk_score_all, file = paste(path, "Scoring_cybeRNAting.RData"))


sugar_scores    <- score_sugar(scores)
chi_scores      <- score_chi(scores)

global_score    <- score_str(scores)

write(c("globalSCORE", "backboneSCORE", "sugarSCORE", "chiSCORE"),
      file=paste(out_prefix, "_globalSCOREs.txt", sep=""),
      append=F, sep="  ", ncolumns=4)
write(c(global_score, 
        mean(backbone_scores, rm.na=T),
        mean(sugar_scores, rm.na=T),
        mean(chi_scores, rm.na=T)),
      file=paste(out_prefix, "_globalSCOREs.txt", sep=""),
      append=T, sep="  ", ncolumns=4)

## Save new pdb file
pdb <- replace_pdb(backbone_scores, pdb, "b")
pdb <- replace_pdb(sugar_scores, pdb, "o")
write.pdb(pdb, file=paste(out_prefix, "_color.pdb", sep=""))

## Find outliers and save them
outliers <- find_outliers(scores, threshold=threshold, ntinfo=torsionals)
write.table(outliers, 
            file=paste(out_prefix, "_outliers.txt",sep=""),
            row.names=F, col.names=T)
