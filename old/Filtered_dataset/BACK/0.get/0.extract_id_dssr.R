###Load proper dependencies
library(minipack)
library(bio3d)
library(veriNA3d)

path <- "/Users/eric1999/Desktop/Gitlab/trinucleotides/Mining/rmsdpdbs/04.torsionals/"
vpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/BACK/0.get/Data/"
all_files <- list.files(path, pattern="repres")

#Load the first ntinfo
## No em funciona Eric, disenyo estrategia alternativa.
#load(paste(path, "ntinfo_1.RData", sep=""))

data <- vector("list", length(all_files))

path_mod <- "/orozco/projects/PDB.datamining/rmsd_data2/"
ntinfo <- list()
count <- 1
for(i in 1:length(all_files)){
  init <- strsplit(all_files[i], split="representatives")
  init <- init[[1]][1]
  
  #Check all data in the file
  file <- readLines(paste(path, all_files[i], sep=""))
  data[[i]] <- file
}

PDBID_id_dssr <- unlist(data)

save(PDBID_id_dssr, file = paste(vpath, "PDBID_id_dssr.RData", sep=""))

