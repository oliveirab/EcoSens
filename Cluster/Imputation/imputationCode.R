### R code for trait imputation using Rphylopars 

# library(devtools)
# install_github("ericgoolsby/Rphylopars@devel", build_vignette = TRUE)

list.of.packages <- c("raster","maps","nlme","spdep","rgdal","viridis","ggplot2","rworldmap","gridExtra", "cowplot","effects","PerformanceAnalytics","dplyr","Rphylopars","pbapply","parallel","doParallel","reshape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

phy.imput <- readRDS("phy_plants.RDS")
traits.2imput <- readRDS("traits_plants.RDS")

# species names in trait data should match tips on the phylogeny
rownames(traits.2imput) <- traits.2imput[,1]
traits.2imput <- traits.2imput[phy.imput[[1]]$tip.label,]
names(traits.2imput)[1] <- "species"

# log-transform traits
traits.2imput[,-1] <- log(traits.2imput[,-1])

# 3 decimal places for trait values
# traits.2imput[,-1] <- round(traits.2imput[,-1],3)

# run imputation
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)

foreach(i=1:length(phy.imput), .packages=c('Rphylopars')) %dopar% { 
  phy.tmp <- phy.imput[[i]]
  phy.tmp$edge.length[which(phy.tmp$edge.length<0.0001)]=0.0001
    
  tmp.traitsimput <- phylopars(traits.2imput, phy.tmp)
  
  saveRDS(tmp.traitsimput, paste("Results/traitsimput_", i, ".RDS", sep = ""))
}
stopCluster(cl)

###############
