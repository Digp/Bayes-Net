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
opt$id <- "1bau"
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

write.pdb(pdb = pdb_new, file= paste(spath, opt$id,"_bayes.pdb", sep=""), segid = pdb_new$atom$entid)