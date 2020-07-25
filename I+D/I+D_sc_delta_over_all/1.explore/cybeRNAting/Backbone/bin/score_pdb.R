#!/usr/local/bin/Rscript
## Diego Gallego
## Date: 2018-Feb-26
## This script takes a pdb file and computes the scoring funciton on it, then
## the output is saved in the pdb format and different informative .txt files
   
## Check correct usage

path <- "/Users/eric1999/Desktop/Bayessian/"
cpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/I+D/DATA/"
#CLIargs <- c("1a9n", "/Users/eric1999/Desktop/cybeRNAting/Scripts/Backbone/dat/referencedata.txt", "scoring_function.R", 0.2, "F")


## Load dependencies ...
library(bio3d)
library(veriNA3d)
library(magrittr)
library(BBmisc)

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
wpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/Data/"

#Load data on the whole dataset that we have from all the RNA structures
load(paste(wpath, "ntinfo_id_dssr.RData", sep=""))

#Retreive important nucleotide angles
ntinfo_m <- ntinfo

torsionals <- ntinfo_m

Ids <- as.character(droplevels.factor(ntinfo$ID_DSSR))

## Compute scores and save summary

scores <- scoring_function(torsionals, referencedata)

#Retreive data of interest
delta_score <- score_delta(scores)


#Forget about NaN to compute min and max values
ind <- which(is.na(delta_score))
data_cybeRNAting <- delta_score[-ind]

#Extract the min and the max values
min_cybeRNAting <- min(data_cybeRNAting)
max_cybeRNAting <- max(data_cybeRNAting)

min_max_cybeRNAting <- cbind(min_cybeRNAting, max_cybeRNAting)


#Normalize values
data_cybeRNAting<- c()

for(i in 1:length(delta_score)){
  data_cybeRNAting[i] <- (delta_score[i] - min_cybeRNAting) / (max_cybeRNAting - min_cybeRNAting)
}

data_cybeRNAting <- round(data_cybeRNAting, 3)




#Standarize obtained values
#stand <- normalize(as.numeric(bulk_score_all$Scoring), method = "standardize")

#bulk_score_all <- cbind(data_cy)
#bulk_score_all <- as.data.frame(bulk_score_all)
#namy <- c("ID_DSSR", "Scoring", "Normalized Val")
#names(bulk_score_all) <- namy

save(min_max_cybeRNAting, data_cybeRNAting, delta_score, file = paste(cpath, "Scoring_cybeRNAting.RData"))


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
