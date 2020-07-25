#!/usr/bin/Rscript

## ----------------------------------------------------------------------------
## Read arguments
## ----------------------------------------------------------------------------
library(optparse)

option_list = list(
    make_option(c("-pdb", "--ipath"), type="character", default=NULL,
                help="Path to load the .pdb or .cif files ",
                metavar="character"),
    make_option(c("-o", "--opath"), type="character", default=NULL,
                help="Path to save output Rdata file",
                metavar="character"),
    make_option(c("-o", "--opath"), type="character", default=NULL,
                help="Path to save output Rdata file",
                metavar="character")
    
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$ipath) | is.null(opt$opath)) {
    print_help(opt_parser)
    stop("Please, provide necessary arguments", call.=FALSE)
}

path <- opt$vpath
opath <- opt$opath
path <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/textfiles/"
# opt <- list(vpath="textfiles/", opath="Data/")

## ----------------------------------------------------------------------------
## Load dependencies and helper functions
## ----------------------------------------------------------------------------
library(bnlearn)
library(minipack)
library(bio3d)
library(veriNA3d)
library(lattice)
library(rbmn)
#load("validation.R")  ##Can be obtained from the minipack

## ----------------------------------------------------------------------------
## Get validations data
## ----------------------------------------------------------------------------

####################Data loading
#Get the RNA list with maximum threshold of 2.5A
rna_list <- getRNAList(release=3.94, threshold="2.5A", as.df=TRUE)
rna_chain <- paste(rna_list$pdb, rna_list$chain, sep=":")
rna_cm <- paste(rna_chain, rna_list$model, sep=":")


#Extract PDB IDs from the rna_list
pdbID <- rna_list$pdb
pdbID <- unique(pdbID)

##################Validation
#Validation data on the different PDBs

for (j in 1:length(pdbID)){
  pdbID2 <- pdbID[j]
  filename <- paste(path, pdbID2, ".txt", sep = "")
  if (file.exists(filename)) {
    next()
  }

  tryCatch({
    #print(validation(pdbID2))
    write.table(validation(pdbID2), filename, row.names = TRUE, col.names = TRUE)
  }, error=function(e){
    print(paste("We have a problem in", pdbID[j], sep = " "))
    e
  })
}


#Read the above mentioned files

#Make it for one
pdbIDs <- pdbID[1]
tint <- read.table(paste(path, pdbIDs, ".txt", sep=""), header = TRUE, sep = "")
tint <- as.data.frame(tint)
tint$value_rsrz <- NULL
tint$PDBID <- pdbIDs

#Loop over all structures
for(i in 1:length(pdbID)){
  nuc <- read.table(paste(path, pdbID[i], ".txt", sep=""), header = TRUE, sep = "")
  nuc <- as.data.frame(nuc)
  nuc_d <- c()
  nuc_d$seq <- nuc$seq
  nuc_d$id_dssr <- nuc$id_dssr
  nuc_d$pucker_outlier <- nuc$pucker_outlier
  nuc_d$suite_outlier <- nuc$suite_outlier
  nuc_d$clashes <- nuc$clashes
  nuc_d$bond_lengths <- nuc$bond_lengths
  nuc_d$bond_angles <- nuc$bond_angles
  nuc_d$chirals <- nuc$chirals
  nuc_d$planes <- nuc$planes
  nuc_d$rsrz <- nuc$rsrz
  #nuc_d$symm_clashes <- nuc$symm_clashes
  nuc_d <- as.data.frame(nuc_d)
  #nuc_e <- nuc_d[,-11:-1]
  nuc_e <- nuc_d
  nuc_e$PDBID <- pdbID[i]
  names(nuc_e) <- names(tint)
  tint <- rbind(tint, nuc_e)
}

#Obtain all different chains from each nucleotide
for (i in 1:length(tint$seq)){
  ant <- strsplit(as.character(tint$id_dssr[i]), split=":")
  int <- strsplit(ant[[1]][2], split=".", fixed = TRUE)
  tint$CHAIN[i] <- int[[1]][1]
}

#Obtain all different models for each nucleotide
for (i in 1:length(tint$seq)){
  ant <- strsplit(as.character(tint$id_dssr[i]), split=":")
  tint$MODEL[i] <- ant[[1]][1]
}

tint$PDB_CHAIN <- paste(tint$PDBID, tint$CHAIN, sep=":")
tint$PDB_CM <- paste(tint$PDB_CHAIN, tint$MODEL, sep=":")

#Retreive data which is 2.9 A

ind <- which(tint$PDB_CM %in% rna_cm[1])
ping <- tint[ind,]

for(i in 1:length(rna_chain)){
  index <- which(tint$PDB_CM %in% rna_cm[i])
  if(length(index >0)){
    ping <- rbind(ping, tint[index,] )
  }
}

#Validate which PDB ID do not have any miss at the validations
#pdbID_2 <- unique(ping$PDBID)

#for(i in 1:length(pdbID_2)){
#  ind <- which(ping$PDBID %in% pdbID_2[i])
#}


#Validations 

imp_suite <- which(ping$suite_outlier %in% TRUE)
imp_pucker <- which(ping$pucker_outlier %in% TRUE)
imp_clashes <- which(ping$clashes %in% TRUE)
imp_bond_l <- which(ping$bond_lengths %in% TRUE)
imp_bond_a <- which(ping$bond_angles %in% TRUE)
imp_chirals <- which(ping$chirals %in% TRUE)
imp_planes <- which(ping$planes %in% TRUE)
imp_rsrz <- which(ping$rsrz %in% TRUE)


append_list <- c(imp_suite, imp_pucker, imp_clashes, imp_bond_l, imp_bond_a, imp_chirals, imp_planes, imp_rsrz)
uni_app_l <- unique(append_list)

#Select specific rows with no TRUE
val_nuc <- ping[-uni_app_l,]

for (i in 1:length(val_nuc$seq)){
  id_split <- strsplit(as.character(val_nuc$id_dssr[i]), split=".", fixed=TRUE)
  take <- id_split[[1]][2]
  chm <- paste(val_nuc$MODEL[i], val_nuc$PDBID[i], sep=":")
  chain <- paste(chm, val_nuc$CHAIN[i], sep=":")
  all <- paste(chain, take, sep=".")
  val_nuc$IDDSSR[i] <- all
}

path2 <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/Data/"

#save(val_nuc, file = "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/Data/Validation.RData")


## ----------------------------------------------------------------------------
## Get torsional data and clean it based on validations
## ----------------------------------------------------------------------------

#################Process validation output onto main dataset
#Obtain data on all previously validated pdbID
ntinfo <- pipeNucData(pdbID = rna_list$pdb, model = rna_list$model, chain = rna_list$chain )
id <- id_dssr_extract(ntinfo)

#save(ntinfo, file = "/Users/eric1999/Desktop/Bayessian/All_nucleotide.RData")
#load("/Users/eric1999/Desktop/Bayessian/All_nucleotide.RData")


######Extract id_dssr information unique identifiers

#Obtain sequence
seq <- ntinfo$resid

# Creation of resid_resno
resid_resno <- paste(ntinfo$resid, ntinfo$resno, sep="")

# Creation of chain.residresno
chain_resid_resno <- paste(ntinfo$chain, resid_resno, sep=".")

#Creation of id_dssr <- model:chain.residresno
id_dssr <- paste(ntinfo$model, chain_resid_resno, sep=":")

count <- 0

for(k in ntinfo$insert){
  if(k =='?'){
    # do not do anything if that happens
    count = count + 1
  } else {
    # paste into position i of vector m
    id_dssr[count] <- paste(id_dssr[count], k, sep='^')
    count = count +1
  }
  data_table <- as.data.frame(cbind(seq, id_dssr))
}

#save(ntinfo, file = "/Users/eric1999/Desktop/Bayessian/id_dssr_total.RData")
#load("/Users/eric1999/Desktop/Bayessian/Bayessian/id_dssr_total.RData")
#################################

#Check the proper nucleotides in the folder
#ntinfo$chain <- paste(ntinfo$pdbID, ntinfo$chain, sep=":")
#ntinfo$cm <- paste(ntinfo$chain, ntinfo$model, sep=":")
ntinfo$ID_DSSR <- data_table$id_dssr

#save(ntinfo, file = "/Users/eric1999/Desktop/Bayessian/ntinfo_id_dssr.RData")
load("/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/Data/ntinfo_id_dssr.RData")

#Take the validation indexes
indexes <- which(val_nuc$IDDSSR %in% ntinfo$ID_DSSR)
ntinfo_l <- ntinfo[indexes,]

################################
#Retreive important nucleotide angles
ntinfo_m <- ntinfo_l[, 66:71]
ntinfo_m$chi <- ntinfo_l$chi
#ntinfo_l[is.na(ntinfo_l)] <- 0
#ntinfo_m <- na.omit(ntinfo_m)

#Try out the discrete function for bnlearn

#disc <- discretize(ntinfo_m, method = "quantile")
#pdag = iamb(gaussian.test)

#dag = pdag2dag(pdag, ordering = c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"))
#fit <- bn.fit(dag, disc)

#Process of discretize the values according to frequency

#disc <- discretize(ntinfo_m = ntinfo_m, windows = 15)   #No need to run as it is integrated in the probability() function
windows <- 3
ntinfo_disc <- discretize(ntinfo_m = ntinfo_m, windows = windows)

#Check for all different angles

alpha <- probability(angle = "alpha", ntinfo_disc = ntinfo_disc)
beta <- probability(angle = "beta", ntinfo_disc = ntinfo_disc)
gamma <- probability(angle = "gamma", ntinfo_disc = ntinfo_disc)
delta <- probability(angle = "delta", ntinfo_disc = ntinfo_disc)
epsilon <- probability(angle = "epsilon", ntinfo_disc = ntinfo_disc)
zeta <- probability(angle = "zeta", ntinfo_disc = ntinfo_disc)
chi <- probability(angle = "chi", ntinfo_disc = ntinfo_disc)
bayesian_data <- list(windows=windows, alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, zeta=zeta, chi=chi)

save(bayesian_data, file = paste(path, "Reference_3_degree.RData", sep=""))

#In case we want to range

path <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/Data/"
range <- 2:10 

for(i in range){
  #disc <- discretize(ntinfo_m = ntinfo_m, windows = 15)   #No need to run as it is integrated in the probability() function
  windows <- i
  ntinfo_disc <- discretize(ntinfo_m = ntinfo_m, windows = windows)
  
  #Check for all different angles
  
  alpha <- probability(angle = "alpha", ntinfo_disc = ntinfo_disc)
  beta <- probability(angle = "beta", ntinfo_disc = ntinfo_disc)
  gamma <- probability(angle = "gamma", ntinfo_disc = ntinfo_disc)
  delta <- probability(angle = "delta", ntinfo_disc = ntinfo_disc)
  epsilon <- probability(angle = "epsilon", ntinfo_disc = ntinfo_disc)
  zeta <- probability(angle = "zeta", ntinfo_disc = ntinfo_disc)
  chi <- probability(angle = "chi", ntinfo_disc = ntinfo_disc)
  bayesian_data <- list(windows=windows, alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, zeta=zeta, chi=chi)
  
  save(bayesian_data, file = paste(path, "Reference_", i, "_degree.RData", sep=""))
}

