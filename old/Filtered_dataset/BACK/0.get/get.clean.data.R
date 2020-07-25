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
#path <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/BACK/0.get/Data/"
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
path <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/BACK/0.get/Data/"
## ----------------------------------------------------------------------------
## Get validated and curated data from filtering in cluster
## ----------------------------------------------------------------------------

load(paste(path, "ntinfo_filtered.RData"))
######Extract id_dssr information unique identifiers
ntinfo_l <- ntinfo
################################
#Retreive important nucleotide angles
ntinfo_m <- ntinfo_l[, 32:37]
ntinfo_m$chi <- ntinfo_l$chi
#ntinfo_l[is.na(ntinfo_l)] <- 0
#ntinfo_m <- na.omit(ntinfo_m)
ntinfo_m$alpha <- as.numeric(as.character(ntinfo_m$alpha))
ntinfo_m$beta <- as.numeric(as.character(ntinfo_m$beta))
ntinfo_m$gamma <- as.numeric(as.character(ntinfo_m$gamma))
ntinfo_m$delta <- as.numeric(as.character(ntinfo_m$delta))
ntinfo_m$epsilon <- as.numeric(as.character(ntinfo_m$epsilon))
ntinfo_m$zeta <- as.numeric(as.character(ntinfo_m$zeta))
ntinfo_m$chi <- as.numeric(as.character(ntinfo_m$chi))

#Try out the discrete function for bnlearn

#disc <- discretize(ntinfo_m, method = "quantile")
#pdag = iamb(gaussian.test)

#dag = pdag2dag(pdag, ordering = c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"))
#fit <- bn.fit(dag, disc)

#Process of discretize the values according to frequency

#disc <- discretize(ntinfo_m = ntinfo_m, windows = 15)   #No need to run as it is integrated in the probability() function
windows <- 4
ntinfo_disc <- minipack::discretize(ntinfo_m = ntinfo_m, windows = windows)

#Check for all different angles

alpha <- minipack::probability(angle = "alpha", ntinfo_disc = ntinfo_disc)
beta <- minipack::probability(angle = "beta", ntinfo_disc = ntinfo_disc)
gamma <- minipack::probability(angle = "gamma", ntinfo_disc = ntinfo_disc)
delta <- minipack::probability(angle = "delta", ntinfo_disc = ntinfo_disc)
epsilon <- minipack::probability(angle = "epsilon", ntinfo_disc = ntinfo_disc)
zeta <- minipack::probability(angle = "zeta", ntinfo_disc = ntinfo_disc)
chi <- minipack::probability(angle = "chi", ntinfo_disc = ntinfo_disc)
bayesian_data <- list(windows=windows, alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, zeta=zeta, chi=chi)

save(bayesian_data, file = paste(path, "Reference_4_degree.RData", sep=""))


###############CODE FROM HERE BELOW WAS NOT RUN
#In case we want to range

#path <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/Data/"
#range <- 2:10 

#for(i in range){
#  #disc <- discretize(ntinfo_m = ntinfo_m, windows = 15)   #No need to run as it is integrated in the probability() function
#  windows <- i
#  ntinfo_disc <- discretize(ntinfo_m = ntinfo_m, windows = windows)
  
  #Check for all different angles
  
#  alpha <- probability(angle = "alpha", ntinfo_disc = ntinfo_disc)
#  beta <- probability(angle = "beta", ntinfo_disc = ntinfo_disc)
#  gamma <- probability(angle = "gamma", ntinfo_disc = ntinfo_disc)
#  delta <- probability(angle = "delta", ntinfo_disc = ntinfo_disc)
#  epsilon <- probability(angle = "epsilon", ntinfo_disc = ntinfo_disc)
#  zeta <- probability(angle = "zeta", ntinfo_disc = ntinfo_disc)
#  chi <- probability(angle = "chi", ntinfo_disc = ntinfo_disc)
#  bayesian_data <- list(windows=windows, alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, zeta=zeta, chi=chi)
#  
#  save(bayesian_data, file = paste(path, "Reference_", i, "_degree.RData", sep=""))
#}

