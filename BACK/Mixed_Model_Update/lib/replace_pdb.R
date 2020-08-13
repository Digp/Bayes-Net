replace_pdb <- function(score_vec, pdb, col="b") {
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
