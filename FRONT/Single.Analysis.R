#!/usr/local/bin/Rscript
#!/usr/bin/Rscript

## ----------------------------------------------------------------------------
## Read arguments
## ----------------------------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-p", "--pdb"), type="character", default=NULL,
              help="Path to your pdb file",
              metavar="character"),
  make_option(c("-c", "--cif"), type="character", default=NULL,
              help="Path to your cif file",
              metavar="character"),
  make_option(c("-i", "--id"), type="character", default=NULL,
              help="Input your PDB ID ",
              metavar="character"),
  make_option(c("-r", "--rpath"), type="character", default=NULL,
              help="Path to the file Scoring_Analysis.RData ",
              metavar="character"),
  make_option(c("-d", "--dpath"), type="character", default=NULL,
              help="Path of the Reference Data file ",
              metavar="character"),
  make_option(c("-s", "--spath"), type="character", default=NULL,
              help="Path to the save the new PDB file ",
              metavar="character")

)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$rpath) | is.null(opt$dpath) | is.null(opt$spath)) {
  print_help(opt_parser)
  stop("Please, provide necessary arguments", call.=FALSE)
}

#is.null(opt$pdb) || is.null(opt$cif) || is.null(opt$id)

rpath <- opt$rpath
dpath <- opt$dpath
spath <- opt$spath


######Load libraries
library(minipack)
library(bio3d)
library(veriNA3d)
library(lattice)
library(rbmn)
library(BBmisc)
library("tools")
####################
opt <- c()
#vpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/FRONT/1a9n.pdb"
#rpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/FRONT/Data/"
#dpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/BACK/0.get/Data/"
#spath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/FRONT/New_PDBs/"
#mpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/FRONT/Good_PDBs/"
opt$id <- "1NA2"
opt$pdb <- NULL

####Read the file
if(!is.null(opt$pdb)){
  input <- read.pdb(opt$pdb)
}else if(!is.null(opt$cif)){
  input <- cifParser(opt$cif)
  input <- cifAsPDB(input)
}else if (!is.null(opt$id)){
  input <- cifParser(opt$id)
  input <- cifAsPDB(input)
}

#Get the coordinates and angles

file_measure <- measureNuc(input)
file <- file_measure[,32:37]
file$chi <- file_measure$chi



#Load the Scoring data and the Reference Data

#load(paste(rpath, " Scoring_Analysis.RData", sep=""))
#load(paste(dpath, "Reference_10_degree.RData", sep=""))
#load(rpath)
#load(dpath)

file_disc <- minipack::discretize(ntinfo_m = file, windows = 1)

#Calculate the score for the values of the PDBID

#######################

bayes_back <- scoring_bayes(ntinfo_disc = file_disc[1,], bd = bayesian_data,ntinfo = file[1,])
bayes_back <- as.data.frame(bayes_back)

for(i in 2:length(file$alpha)){
  amb <- scoring_bayes(ntinfo_disc = file_disc[i,], bd = bayesian_data,ntinfo = file[i,])
  amb <- as.data.frame(amb)
  bayes_back <- rbind(bayes_back, amb)
}

bayes_back <- round(bayes_back, 3)

## Save new pdb file
pdb_new <- replace_pdb(bayes_back$total, input, "b")

write.pdb(pdb = pdb_new, file= paste(mpath, opt$id,"_bayes.pdb", sep=""), segid = pdb_new$atom$entid)

#######################

back_score <- scoring(ntinfo_disc = file_disc[1,], bd = bayesian_data,ntinfo = file[1,])
back_score <- as.data.frame(back_score)
#Ids <- open_file$id_dssr

for(i in 2:length(file$alpha)){
  amb <- scoring(ntinfo_disc = file_disc[i,], bd = bayesian_data,ntinfo = file[i,])
  amb <- as.data.frame(amb)
  back_score <- rbind(back_score, amb)
}

#Normalize data
back_score$alpha <- normalize(data = back_score$alpha)
back_score$beta <- normalize(data = back_score$beta)
back_score$gamma <- normalize(data = back_score$gamma)
back_score$delta <- normalize(data = back_score$delta)
back_score$epsilon <- normalize(data = back_score$epsilon)
back_score$zeta <- normalize(data = back_score$zeta)
back_score$chi <- normalize(data = back_score$chi)
back_score$total <- normalize(data = back_score$total)

back_score <- round(back_score, 3)

## Save new pdb file
pdb_new <- replace_pdb(back_score$total, input, "b")

write.pdb(pdb = pdb_new, file= paste(mpath, opt$id,"_theorem.pdb", sep=""), segid = pdb_new$atom$entid)

#######BACKTESTING FOR 1BAU##############
#########################################
########################################
#Select south configuration nucleotides from the 1bau
file$pu_round <- round(file_measure$pu_phase,0)
angles_south <- 144:180

#Select those angles in the south conformation
ind_south <- which(file$pu_round %in% angles_south)

file_south <- file[ind_south,]

bayes_back[ind_south,]


angles_nor <- 0:36
ind_nort <- which(file$pu_round %in% angles_nor)
bau_nort <- bayes_back[ind_nort,]
file_nort <- file[ind_nort,]

ind_bad <- which(bayes_back[ind_nort,]$total < 0.5)
ind_good <- which(bayes_back[ind_nort,]$total > 0.75)

bau_nort[ind_bad,]
bau_nort[ind_good,]

mean_good_ang <- c(mean(file_nort[ind_good,]$alpha), mean(file_nort[ind_good,]$beta), mean(file_nort[ind_good,]$gamma), mean(file_nort[ind_good,]$delta),
                   mean(file_nort[ind_good,]$epsilon), mean(file_nort[ind_good,]$zeta), mean(file_nort[ind_good,]$chi))

mean_bad_ang <- c(mean(na.omit(file_nort[ind_bad,]$alpha)), mean(na.omit(file_nort[ind_bad,]$beta)), mean(na.omit(file_nort[ind_bad,]$gamma)),
                  mean(na.omit(file_nort[ind_bad,]$delta)), mean(na.omit(file_nort[ind_bad,]$epsilon)), mean(na.omit(file_nort[ind_bad,]$zeta)),
                  mean(na.omit(file_nort[ind_bad,]$chi)))

mean_analysis <- data.frame(mean_good_ang, mean_bad_ang)
names(mean_analysis) <- c("GOOD ANGLES", "BAD ANGLES")
diff <- mean_analysis$`GOOD ANGLES`- mean_analysis$`BAD ANGLES`
mean_analysis$DIFF <- diff
rownames(mean_analysis) <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi")

############################
#Check nice scoring south angles in the whole ntinfo
path <- "/Users/eric1999/Desktop/Gitlab/trinucleotides/NonRedundantTrinucleotides/rmsdpdbs/04.torsionals/"

ntinfo <- read.table(paste(path, "ntinfo.txt", sep=""), header= TRUE, sep=" ", dec=".")

#Extract the identifiers of interest
ntinfo_m <- c()
ntinfo_m$alpha <- ntinfo[, "alpha"]
ntinfo_m$beta <- ntinfo[, "beta"]
ntinfo_m$gamma <- ntinfo[, "gamma"]
ntinfo_m$delta <- ntinfo[, "delta"]
ntinfo_m$epsilon <- ntinfo[, "epsilon"]
ntinfo_m$zeta <- ntinfo[, "zeta"]
ntinfo_m$chi <- ntinfo[, "chi"]
ntinfo_m$pu_round <- round(ntinfo$pu_phase,0)
ntinfo_m <- as.data.frame(ntinfo_m)
ntinfo <- ntinfo_m

#Extract nucleotides in south
ind_n_so <- which(ntinfo$pu_round %in% angles_south)

ntinfo_south <- ntinfo[ind_n_so,]

#Extract the highest scoring one
ind_high <- which(evaluation_data[ind_n_so,]$total > 0.716)

ntinfo_south[ind_high,]
