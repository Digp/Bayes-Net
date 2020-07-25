library(minipack)
library(bio3d)
library(veriNA3d)

path <- "/Users/eric1999/Desktop/Bayessian_RNA_val/Filtered_dataset/"

###############ntinfo reading and processing 
ntinfo <- read.table(paste(path, "ntinfo.txt", sep=""))
names <- ntinfo[1,]

#Extract names from the first row
data_names <- c()
for(i in 1:length(names)){
  data_names[i] <- as.character(droplevels(names[1,i]))
}


#Set names
ntinfo <- ntinfo[-1,]
names(ntinfo) <- data_names

save(ntinfo, file = paste(path, "ntinfo_filtered.RData"))
