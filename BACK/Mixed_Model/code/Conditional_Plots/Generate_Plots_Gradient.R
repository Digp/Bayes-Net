############PERFORM EXPLORATORY ANALYSIS
library(minipack)
library(veriNA3d)
library(mclust)
library(ggplot2)

######LOAD THE DATA
#Read the curated ntinfo data
path <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/new_ntinfo_retrain/"
path2 <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/BACK/lib/"

load(paste(path2, "Dictionary.RData", sep=""))

ntinfo <- read.table(paste(path, "ntinfo.txt", sep=""), header= TRUE, sep=" ", dec=".")

#Extract the identifiers of interest
ntinfo_m <- c()
ntinfo_m$alpha <- ntinfo[, "alpha"]
ntinfo_m$beta <- ntinfo[, "beta"]
ntinfo_m$gamma <- ntinfo[, "gamma"]
ntinfo_m$delta <- ntinfo[, "delta"]
ntinfo_m$epsilon <- ntinfo[, "epsilon"]
ntinfo_m$zeta <- ntinfo[, "zeta"]
ntinfo_m$chi <- ntinfo[, "chi"]
ntinfo_m <- as.data.frame(ntinfo_m)
ntinfo <- ntinfo_m


#Obtian discretized data
#load(paste(path2, "Data_Discretized.RData", sep=""))

ntinfo_disc <- minipack::discretize(ntinfo_m = ntinfo, windows = 1)

#Perform scoring with obtaining the data for each single nucleotide as data_frame
evaluation_data <- scoring_bayes(ntinfo_disc = ntinfo_disc[1,], bd = bayesian_data, ntinfo = ntinfo[1,])
evaluation_data <- as.data.frame(evaluation_data)

for(i in 2:length(ntinfo_disc[[1]])){
  tryCatch({
    amb <- scoring_bayes(ntinfo_disc = ntinfo_disc[i,], bd = bayesian_data, ntinfo = ntinfo[i,])
    amb <- as.data.frame(amb)
    evaluation_data <- rbind(evaluation_data, amb)
  }, error = function(e) {
    #amb<- data.frame(NA, NA, NA, NA, NA, NA, NA, NA)
    print(i)
  })
}

#Plot
hist(ntinfo$epsilon, breaks = 100)
par(new=T)
plot(x = ntinfo$epsilon, y= as.numeric(as.character(evaluation_data$epsilon)), axes=F)
axis(side=4)

total <- cbind(ntinfo, evaluation_data$delta)

#Discretize data with windowing n = 10

dictio_p <- dict[1:6,3:9]


for(i in 1:length(ntinfo)){
  index <- which(names(dictio_p) %in% names(ntinfo)[i])
  for(m in 1:length(dictio_p[,index])){
    
    tot <- as.character(dictio_p[,index])
    
    #Set the jpeg filepath
    #jpeg(filename = paste("/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/BACK/Conditional_Plots/",names(ntinfo)[i], "/", names(ntinfo)[i], "_", tot[m], sep="" ))
    
    #Create ggplot charts and save
    #
    ggplot(data = ntinfo, aes(x = ntinfo[,i], y = evaluation_data[,i])) + geom_point(aes(colour=ntinfo[,tot[m]]), size=0.1) +
      scale_colour_gradient(low = "blue", high = "red", limits = c(0, 360)) + theme_bw() + xlim(0, 360) + theme(aspect.ratio=1) +
      xlab(paste(names(ntinfo)[i], "angle (Degree)", sep=" ")) +  ylab("Score (a.u)") + 
      labs(title = paste(names(ntinfo)[i], "conditioned to", tot[m], sep=" "), color = tot[m])
    
    ggsave(filename = paste("/Users/eric1999/Desktop/Bayessian_RNA_val/Mixed_Model/BACK/Conditional_Plots/",names(ntinfo)[i], "/", names(ntinfo)[i], "_", tot[m], ".png", sep="" ),
           plot = last_plot())
    
  }
}
