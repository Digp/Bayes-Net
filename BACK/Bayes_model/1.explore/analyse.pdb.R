#!/usr/bin/Rscript

## ----------------------------------------------------------------------------
## Read arguments
## ----------------------------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-v", "--vpath"), type="character", default=NULL,
              help="Path to save Reference data file ",
              metavar="character"),
  make_option(c("-c", "--cpath"), type="character", default=NULL,
              help="Path to cybeRNAting scoring values",
              metavar="character"),
  make_option(c("-w", "--wpath"), type="character", default=NULL,
              help="Path with the ntinfo data to load",
              metavar="character"),
  make_option(c("-s", "--spath"), type="character", default=NULL,
              help="Path to save the resulting scoring values",
              metavar="character"),
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$vpath) | is.null(opt$wpath) | is.null(opt$cpath) | is.null(opt$spath)) {
  print_help(opt_parser)
  stop("Please, provide necessary arguments", call.=FALSE)
}

vpath <- opt$vpath
cpath <- opt$cpath
wpath <- opt$wpath
spath <- opt$spath


###Load proper dependencies
library(minipack)
library(bio3d)
library(veriNA3d)
library(lattice)
library(rbmn)
library(BBmisc)

#Read the desired files
#path <- "/var/folders/3l/xhm3x7fs02g9ydpvmf9_h5l40000gn/T//RtmpEs3WJo/"
vpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/Data/"
wpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/Data/"
cpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/1.explore/"
spath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/1.explore/"


load(paste(vpath, "Reference_10_degree.RData", sep=""))           

#####TEST WITH A SINGLE PDB FILE

#ntinfo <- read_file(path, pdbID = "1a9n")
  
#evaluate <- ntinfo[,1:7]

#Scoring function over all data
#data <- c()

#for(i in 1:length(evaluate$alpha)){
#  data[i] <- scoring(evaluate = evaluate[i,], bd = bayesian_data)
#}

#Bind with the nucleotides

#bulk_score <- cbind(ntinfo$id_dssr, data)
#name <- c("Sequence", "ID_DSSR", "Scoring")
#names(bulk_score) <- name


########RUN FOR ALL DATASET
load(paste(wpath, "ntinfo_id_dssr.RData", sep=""))

#Retreive important nucleotide angles
ntinfo_m <- ntinfo[, 66:71]
ntinfo_m$chi <- ntinfo$chi

Ids <- as.character(droplevels.factor(ntinfo$ID_DSSR))


######EVALUATE SCORING WITH THE SCORING() FUNCTION

data_all <- c()

for(i in 1:length(ntinfo_m$alpha)){
  data_all[i] <- scoring_wochi(evaluate = ntinfo_m[i,], bd = bayesian_data)
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

#Bind ID_DSSR with the score

bulk_score_all <- cbind(Ids, data_all, data_normalized)


#Frequencies
#hist(data_all, breaks = 100)

####EVALUATING SCORE WITH THE SCORING_UPDATED() FUNCTION

data_score_upd <- c()

for(i in 1:length(ntinfo_m$alpha)){
  data_score_upd[i] <- scoring_up_wochi(evaluate = ntinfo_m[i,], bd = bayesian_data)
}

data_score_upd <- round(data_score_upd, 3)

#Forget about NaN to compute min and max values
ind <- which(is.nan(data_score_upd))
data_calc_upd <- data_all[-ind]

#Extract the min and the max values
min_score_upd <- min(data_calc_upd)
max_score_upd <- max(data_calc_upd)

min_max_scoring_update <- cbind(min_score_upd, max_score_upd)

#Normalize values
data_normal_upd<- c()

for(i in 1:length(data_all)){
  data_normal_upd[i] <- (data_calc_upd[i] - min_score_upd) / (max_score_upd - min_score_upd)
}

data_normal_upd <- round(data_normal_upd, 3)


#Load results from backend cybeRNAting
load(paste(cpath, "Scoring_cybeRNAting.RData"))

norm_cybeRNAting <- round(as.numeric(as.character(score_cybeRNAting$Norm_CybeRNAting)), 3)
scoring_CybeRNAting <- round(as.numeric(as.character(score_cybeRNAting$Scoring_CybeRNAting)), 3)

#Bind ID_DSSR with the score

bulk_score_upd <- cbind(Ids, data_all, data_normalized, data_score_upd, data_normal_upd, scoring_CybeRNAting,
                        norm_cybeRNAting)

score_upd_all <- as.data.frame(bulk_score_upd)
namy <- c("ID_DSSR", "Scoring", "Scoring_Normalized", "Scoring_Updated", "Normalized_Scor_Updated", "cybeRNAting_score", "Normalized_cybeRNAting")
names(score_upd_all) <- namy


save(min_max_scoring, min_max_scoring_update, min_max_cybeRNAting, score_upd_all, file = paste(spath, "Scoring_Analysis.RData"))

