###Load proper dependencies
library(minipack)
library(bio3d)
library(veriNA3d)
library(lattice)
library(rbmn)
library(BBmisc)


#Read the desired files
#path <- "/var/folders/3l/xhm3x7fs02g9ydpvmf9_h5l40000gn/T//RtmpEs3WJo/"
vpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/Data/"
wpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/0.get/Data/"
cpath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/1.explore/"
spath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/BACK/1.explore/"
impath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/eta_theta/Data/"


########LOAD NTINFO
load(paste(wpath, "ntinfo_id_dssr.RData", sep=""))

####EXTRACT VALIDATION FROM ETA AND THETA ANGLES

eta_val <- which(grepl("TRUE",ntinfo$eta_valid))
theta_val <- which(grepl("TRUE",ntinfo$theta_valid))

et_th_val <- append(eta_val, theta_val)
et_th <- unique(et_th_val)

ntinfo_m <- ntinfo[et_th,]

#Process data to extract specific data

PDB_ID <- c()

for(i in 1:length(ntinfo_m[[1]])){
  amb <- strsplit(as.character(ntinfo_m$ID_DSSR[i]), split=":")
  pdb <- amb[[1]][2]
  PDB_ID[i] <- pdb
}

#Substrat eta, theta and ID_DSSR information and apply PDB_ID
ntinfo_l <- c()
ntinfo_l$eta <- ntinfo_m$eta
ntinfo_l$theta <- ntinfo_m$theta
ntinfo_l$PDB_ID <- ntinfo_m$pdbID
ntinfo_l <- as.data.frame(ntinfo_l)
#Separate eta-theta angles for each PDB 

pdbs <- as.character(unique(ntinfo_l$PDB_ID))

#################################
### Extract information for every PDB
seg <- vector("list", length(pdbs))

for(i in 1:length(pdbs)){
  ind <- which(grepl(pdbs[i], ntinfo_l$PDB_ID))
  new_dat <- ntinfo_m[ind,]
  
  #Extract data of interest
  new_file <- c()
  new_file$ID_DSSR <- new_dat$ID_DSSR
  new_file$PDB_IS <- new_dat$pdbID
  new_file$eta <- new_dat$eta
  new_file$theta <- new_dat$theta
  new_file <- as.data.frame(new_file)
  
  
  
  seg[[i]] <- new_file
}


#Select specific nucleotides

for(i in 1:length(pdbs)){
  #n_inicial <- 3
  #nfinal <- length(seg[[i]][[1]]) - 3
  for(mm in 1:length(seg[[i]][[1]])){
    
  #Analyse for the initial angles
    if (mm < 2){
      seg[[i]]$eta_menos_1[mm] <- NA
      seg[[i]]$theta_menos_1[mm] <- NA
      seg[[i]]$eta_mas_1[mm] <- NA
      seg[[i]]$theta_mas_1[mm] <- NA
    }else if (mm >=2){
      seg[[i]]$eta_menos_1[mm] <- seg[[i]]$eta[mm-1]
      seg[[i]]$theta_menos_1[mm] <- seg[[i]]$theta[mm-1]
      seg[[i]]$eta_mas_1[mm] <- seg[[i]]$eta[mm+1]
      seg[[i]]$theta_mas_1[mm] <- seg[[i]]$theta[mm+1]
    }
  }
}

#Remove the first row as everything will be NA
for(i in 1:length(pdbs)){
  seg[[i]] <- seg[[i]][-1,]
}


#######################APPEND DATA TO GENERATE DATAFRAME
data_angles <- c()
data_angles <- seg[[1]]

#Append all
for (i in 2:length(pdbs)){
  data_angles <- rbind(data_angles, seg[[i]])
}

#Round the data
data_select <- data_angles[,3:length(data_angles)]

#Discretize the data
windows <- 4

ntinfo_disc <- minipack::discretize(ntinfo_m = data_select[1,], windows = windows) 

for(i in 2:length(data_select[[1]])){
  a <- minipack::discretize(ntinfo_m = data_select[i,], windows = windows)
  ntinfo_disc <- rbind(ntinfo_disc, a)
}

#Evaluate the angles
eta <- probability(angle = "eta", ntinfo_disc = ntinfo_disc)
theta <- probability(angle = "theta", ntinfo_disc = ntinfo_disc)
theta_mas_1 <- probability(angle = "theta_mas_1", ntinfo_disc = ntinfo_disc)
eta_mas_1 <- probability(angle = "eta_mas_1", ntinfo_disc = ntinfo_disc)
theta_menos_1 <- probability(angle = "theta_menos_1", ntinfo_disc = ntinfo_disc)
eta_menos_1 <- probability(angle = "eta_menos_1", ntinfo_disc = ntinfo_disc)
bayesian_data <- list(windows=windows, eta=eta, theta=theta, theta_mas_1=theta_mas_1, eta_mas_1=eta_mas_1,
                      theta_menos_1=theta_menos_1, eta_menos_1=eta_menos_1)

impath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/eta_theta/Data/"
save(bayesian_data, file = paste(impath, "Bayesian_data.RData", sep=""))
save(data_select,file= paste(impath, "ntinfo_eta_theta.RData", sep=""))
  





  