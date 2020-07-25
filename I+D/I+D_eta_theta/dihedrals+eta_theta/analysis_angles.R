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
ntinfo_l$alpha <- ntinfo_m$alpha
ntinfo_l$beta <- ntinfo_m$beta
ntinfo_l$gamma <- ntinfo_m$gamma
ntinfo_l$delta <- ntinfo_m$delta
ntinfo_l$epsilon <- ntinfo_m$epsilon
ntinfo_l$zeta <- ntinfo_m$zeta
ntinfo_l$chi <- ntinfo_m$chi
ntinfo_l$PDB_ID <- ntinfo_m$pdbID
ntinfo_l <- as.data.frame(ntinfo_l)

#Round the data
ntinfo_l <- ntinfo_l[,1:length(ntinfo_l)-1]
ntinfo_l <- round(ntinfo_l, 0)

save(ntinfo_l, file= "/Users/eric1999/Desktop/Bayessian_RNA_val/dihedra_eta/ntinfo_dihedrals_eta.RData")


#Discretize the data
windows <- 4

ntinfo_disc <- discretize(ntinfo_m = ntinfo_l[1,], windows = windows) 


for(i in 2:length(ntinfo_l[[1]])){
  a <- discretize(ntinfo_m = ntinfo_l[i,], windows = windows)
  ntinfo_disc <- rbind(ntinfo_disc, a)
}

#Evaluate the angles for the code v1
alpha <- probability_all(angle = "alpha", ntinfo_disc = ntinfo_disc)
beta <- probability_all(angle = "beta", ntinfo_disc = ntinfo_disc)
gamma <- probability_all(angle = "gamma", ntinfo_disc = ntinfo_disc)
delta <- probability_all(angle = "delta", ntinfo_disc = ntinfo_disc)
epsilon <- probability_all(angle = "epsilon", ntinfo_disc = ntinfo_disc)
zeta <- probability_all(angle = "zeta", ntinfo_disc = ntinfo_disc)
chi <- probability_all(angle = "chi", ntinfo_disc = ntinfo_disc)
eta <- probability_all(angle = "eta", ntinfo_disc = ntinfo_disc)
theta <- probability_all(angle = "theta", ntinfo_disc = ntinfo_disc)

bayesian_data <- list(windows=windows, alpha=alpha, beta=beta, gamma=gamma,delta=delta, epsilon=epsilon,
                      zeta=zeta, chi=chi, eta=eta, theta=theta)

impath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/dihedra_eta/Evaluate_dihedra_eta_v1/Data/"
save(bayesian_data, file = paste(impath, "Bayesian_data_4_degree_v1.RData", sep=""))

  
#Evaluate the angles for the code v2
alpha <- probability_seg(angle = "alpha", ntinfo_disc = ntinfo_disc)
beta <- probability_seg(angle = "beta", ntinfo_disc = ntinfo_disc)
gamma <- probability_seg(angle = "gamma", ntinfo_disc = ntinfo_disc)
delta <- probability_seg(angle = "delta", ntinfo_disc = ntinfo_disc)
epsilon <- probability_seg(angle = "epsilon", ntinfo_disc = ntinfo_disc)
zeta <- probability_seg(angle = "zeta", ntinfo_disc = ntinfo_disc)
chi <- probability_seg(angle = "chi", ntinfo_disc = ntinfo_disc)
eta <- probability_seg(angle = "eta", ntinfo_disc = ntinfo_disc)
theta <- probability_seg(angle = "tjeta", ntinfo_disc = ntinfo_disc)

bayesian_data_v2 <- list(windows=windows, alpha=alpha, beta=beta, gamma=gamma,delta=delta, epsilon=epsilon,
                      zeta=zeta, chi=chi, eta=eta, theta=theta)

ampath <- "/Users/eric1999/Desktop/Bayessian_RNA_val/dihedra_eta/Evaluate_dihedra_eta_v2/Data/"
save(bayesian_data_v2, file = paste(ampath, "Bayesian_data_4_degree_v2.RData", sep=""))

  