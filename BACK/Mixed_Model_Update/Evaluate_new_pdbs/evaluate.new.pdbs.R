#!/usr/local/bin/Rscript
#!/usr/bin/Rscript

## ----------------------------------------------------------------------------
## Read arguments
## ----------------------------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-p", "--pdb"), type="character", default=NULL,
              help="Path to your path file",
              metavar="character"),
  make_option(c("-c", "--cores"), type="character", default=NULL,
              help="Number of working cores ",
              metavar="character"),
  make_option(c("-d", "--dpath"), type="character", default=NULL,
              help="Path to the required functions ",
              metavar="character"),
  make_option(c("-s", "--spath"), type="character", default=NULL,
              help="Path to the save the new PDB file ",
              metavar="character")
  
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$pdb) | is.null(opt$dpath) | is.null(opt$spath)) {
  print_help(opt_parser)
  stop("Please, provide necessary arguments", call.=FALSE)
}

#is.null(opt$pdb) || is.null(opt$cif) || is.null(opt$id)


######Load libraries
library(minipack)
library(bio3d)
library(veriNA3d)
library(dplyr)
library("tools")
####################
#
#vpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/FRONT/1a9n.pdb"
#rpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/FRONT/Data/"
#dpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/BACK/0.get/Data/"
#spath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/FRONT/New_PDBs/"
#mpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/FRONT/Good_PDBs/"
#opt$id <- "1NA2"

#opt$pdb <- "/Users/eric1999/Desktop/Bayes-Net/analyse.lots.pdbs/barnacle.analysis/Good/"

source(paste(opt$dpath,"dependencies.R", sep="" ))
load(paste(opt$dpath,"Dictionary.RData", sep="" ))
load(paste(opt$dpath,"Bayesian_Data.RData", sep="" ))

files_pdb <- list.files(opt$pdb, pattern = ".pdb")

scoring_complete <- vector("list", length(files_pdb))
bd <- bayesian_data

#scoring_complete <- lapply(1:length(files_pdb), function(m) {
  
  
for(m in 1:length(files_pdb)){
  
  ####Read the file
  print(files_pdb[m])
  
  pdb.data <- read.pdb(paste(opt$pdb, files_pdb[m], sep = ""))
  
  file <- measureNuc(pdb.data)[, c("alpha", "beta", "gamma", "delta", "epsilon","zeta","chi")]
  
  total_complete <- vector("list", length(file[[1]]))
  
  for(t in 1:length(file[[1]])){
    t_plus <- file[t+1,]
    if(dim(t_plus)[1] == 0){
      t_plus <- c(NA,NA,NA,NA,NA,NA,NA)
    }
    t_minus <- file[t-1,]
    if(dim(t_minus)[1] == 0){
      t_minus <- c(NA,NA,NA,NA,NA,NA,NA)

    }
    names(t_plus) <- paste0(names(file[t,]), "_plus")
    names(t_minus) <- paste0(names(file[t,]), "_minus")
    t_normal <- file[t,]
    total <- c(t_minus, t_normal, t_plus)
    total_complete[[t]] <- as.data.frame(total)
  }
  
  file <- bind_rows(total_complete)
    
  #file <- measureNuc(pdb=read.pdb(paste(opt$pdb, files_pdb[m], sep = "")), torsionals=torsionals)
  
  file_disc <- round(file, 0)
  
  #Calculate the score for the values of the PDBID
  
  #######################
  #bayes_back <- lapply(1:length(file$alpha), function(i) {
  #  round(scoring_bayes(ntinfo_disc = file_disc[i,], bd = bayesian_data,ntinfo = file[i,]), 3)
  #})
  
  
  bayes_back <- vector("list", length(file[[1]]))
  for(z in 1:length(file$alpha)){
    bayes_back[[z]] <- scoring_bayes(ntinfo_disc = file_disc[z,], bd = bayesian_data,ntinfo = file[z,])
  }
  
  bayes_back <- bind_rows(bayes_back)
  bayes_back <- data.frame(round(bayes_back, 3))
  
  #scoring_complete[[m]] <- bayes_back
  #return(bayes_back)
  
  replace.data <- replace_pdb(bayes_back$total, pdb.data, col = "b")
  
  write.pdb(replace.data, file = paste0(opt$pdb, "new_", files_pdb[m], ".pdb"))
  
  
}

#library(dplyr)
#system.time(out <- bind_rows(scoring_complete))
#system.time(out <- do.call(rbind, scoring_complete))



