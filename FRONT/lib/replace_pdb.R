## Replace pdb column o or b by given scores
replace_pdb <-
  function(score_vec, pdb, col="b") {
    nt_vec <- paste(pdb$atom$chain, pdb$atom$resno, pdb$atom$insert, sep="_")
    ntID_equivalences <- unique(nt_vec)
    
    ntID_vec <- unlist(lapply(1:nrow(pdb$atom), 
                              function(i) which(nt_vec[i] == ntID_equivalences)))
    
    replacement <- score_vec[ntID_vec]
    
    pdb$atom[, col] <- replacement
    
    if (any(is.na(pdb$atom[, col]))) {
      ind <- which(is.na(pdb$atom[, col]))
      pdb$atom[ind, col] <- 0
    }
    
    return(pdb)
  }


#Normalization function
normalize <- function(data = NULL){
  max_v <- max(na.omit(data))
  min_v <- min(na.omit(data))
  
  data_norm <- c()
  for(i in 1:length(data)){
    data_norm[i] <- (data[i] - min_v) / (max_v - min_v)
  }
  
  return(data_norm)
}
