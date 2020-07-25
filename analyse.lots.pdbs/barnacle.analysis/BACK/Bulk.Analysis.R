#!/orozco/homes/pluto/ematamoros/Programs/bin/miniconda2/envs/r_env/bin/Rscript

# ----------------------------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-p", "--pdb"), type="character", default=NULL,
              help="Path to your path file",
              metavar="character"),
  make_option(c("-c", "--cif"), type="character", default=NULL,
              help="Path to your cif file",
              metavar="character"),
  make_option(c("-i", "--id"), type="character", default=NULL,
              help="Input your PDB ID ",
              metavar="character"),
  make_option(c("-f", "--fpath"), type="character", default=NULL,
              help="Path to the required functions ",
              metavar="character"),
  make_option(c("-s", "--spath"), type="character", default=NULL,
              help="Path to the save the new PDB file ",
              metavar="character")
  
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$pdb) | is.null(opt$fpath) | is.null(opt$spath)) {
  print_help(opt_parser)
  stop("Please, provide necessary arguments", call.=FALSE)
}

#is.null(opt$pdb) || is.null(opt$cif) || is.null(opt$id)


######Load libraries
library(minipack)
library(bio3d)
library(veriNA3d)
library("tools")
####################
#opt <- c()
#opt$pdb <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Wide.pdb.analysis/Good/"
#vpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/FRONT/1a9n.pdb"
#rpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/FRONT/Data/"
#dpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/BACK/0.get/Data/"
#spath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/FRONT/New_PDBs/"
#mpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/FRONT/Good_PDBs/"
#opt$id <- "1NA2"
#opt$pdb <- NULL

source(paste(opt$fpath,"Scoring_bayes.R", sep="" ))
source(paste(opt$fpath,"dnorm_functions.R", sep="" ))
load(paste(opt$fpath,"Bayesian_Data.RData", sep="" ))
load(paste(opt$fpath,"Dictionary.RData", sep="" ))

discretize <- function (ntinfo_m = NULL, windows = NULL) 
{
  ntinfo_m <- round(ntinfo_m, 0)
  interval <- seq(from = 0, to = 360, by = windows)
  inter_list <- c()
  for (i in 1:length(interval)) {
    tryCatch({
      vect <- (interval[i] + 1):interval[i + 1]
      inter_list[[i]] <- vect
    }, error = function(e) {
    })
  }
  indices <- length(inter_list)
  for (j in 1:length(ntinfo_m$alpha)) {
    for (i in 1:length(ntinfo_m)) {
      la <- which(NA %in% ntinfo_m[j, i])
      if (length(la) == 0) {
        for (m in 1:length(inter_list)) {
          ind <- which(ntinfo_m[j, i] %in% inter_list[m][[1]])
          if (length(ind) > 0) {
            ntinfo_m[j, i] <- m
          }
        }
      }
    }
  }
  return(ntinfo_m)
}

#Load function discretize
files_pdb <- list.files(opt$pdb, pattern = ".pdb")
scoring_complete <- vector("list", length(files_pdb))
files_complete <- c()
       
for(i in 1:length(files_pdb)){
  
  pdbs_path <- paste(opt$pdb, files_pdb[i], sep = "")
  
  ####Read the file
  if(!is.null(opt$pdb)){
    input <- read.pdb(pdbs_path)
  }else if(!is.null(opt$cif)){
    input <- cifParser(opt$cif)
    input <- cifAsPDB(input)
  }else if (!is.null(opt$id)){
    input <- cifParser(opt$id)
    input <- cifAsPDB(input)
  }
  
  ID <- strsplit(pdbs_path, split="/", fixed=T)[[1]]
  ID <- ID[length(ID)]
  ID_f <- strsplit(ID, split=".pdb", fixed="")[[1]][1]
  #Get the coordinates and angles
 
  file_measure <- measureNuc(input)
  file <- file_measure[,32:37]
  file$chi <- file_measure$chi
  
  bd = bayesian_data
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
  
  for(m in 2:length(file$alpha)){
    amb <- scoring_bayes(ntinfo_disc = file_disc[m,], bd = bayesian_data,ntinfo = file[m,])
    amb <- as.data.frame(amb)
    bayes_back <- rbind(bayes_back, amb)
  }
  
  bayes_back <- round(bayes_back, 3)
  
  ## Save new pdb file
 # pdb_new <- replace_pdb(bayes_back$total, input, "b")
  
  files_complete[i] <- files_pdb[i]
  print(files_pdb[i])
  scoring_complete[[i]] <- bayes_back$total
}

#scoring_complete <- do.call(rbind, scoring_complete)

#data <- data.frame(files_complete, scoring_complete)
#names(data) <- c("PDB ID", "Scoring")

save(scoring_complete, file = paste(opt$spath,"Scoring.RData", sep=""))
save(files_complete, file = paste(opt$spath,"PDB_ID.RData", sep=""))

#write.table(data,file = paste(opt$spath, "Analysis_ROSSETA.txt", sep=""),row.names=FALSE)
  

 # write.pdb(pdb = pdb_new, file= paste(opt$spath, ID_f,"_out.pdb", sep=""), segid = pdb_new$atom$entid)

