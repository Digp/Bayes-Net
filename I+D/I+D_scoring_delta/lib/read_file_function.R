#Extract the PDBID and all the id_dssr and pseudo-torsional angles

read_file <- function(path, pdbID = NULL, pdb = NULL, cif = NULL){
  if(!is.null(pdb)){
    file_pdb <- read.pdb(paste(path, tolower(pdb), ".pdb", sep=""))
    file_measure <- measureNuc(file_pdb)
    file <- file[,32:37]
    file$chi <- file_measure$chi
    file$id_dssr <- id_dssr_extract(file_pdb)
    
    
  }else if(!is.null(cif)){
    file_cif <- read.cif(paste(path, tolower(pdbID), ".cif.gz", sep="" ))
    file_pdb <- measureNuc(file_cif) #Buggy, not working 
    file$id_dssr <- id_dssr_extract(file_pdb)
    file <- file_pdb[, 66:71]
    file$chi <- file_pdb$chi
    
  }else if(!is.null(pdbID)){
    file_pdb <- pipeNucData(toupper(pdbID))
    file <- file_pdb[, 66:71]
    file$chi <- file_pdb$chi
    file$id_dssr <- id_dssr_extract(ntinfo = file_pdb)
  }
  return(file)
}
