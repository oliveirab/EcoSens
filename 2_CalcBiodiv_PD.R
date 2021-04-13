
# Load packages

rm(list=ls())
gc()

list.of.packages <- c("raster","maps","nlme","spdep","rgdal","viridis","ggplot2","rworldmap","gridExtra","cowplot","effects","PerformanceAnalytics","dplyr","Rphylopars","pbapply","parallel","doParallel","reshape","picante","FD","corrplot","ggvenn","WorldFlora")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

# Source code

# Calc PD function
calc_PD <- function(i, occu.long, phylo){
  # species that occur in the cell i
  sps.in.cell <- occu.long$sps[which(occu.long$cell==cells[i])] 
  # construct tree of all species within a cell
  trx <- keep.tip(phylo, sps.in.cell) 
  # calc PD
  sum(trx$edge.length)
}

# Calc MPD function
# cophe_phyd = cophenetic distance matrix
calc_MPD <- function(i, occu.long, cophe_phyd){
  # species that occur in the cell i
  sps.in.cell <- occu.long$sps[which(occu.long$cell==cells[i])] 
  # crop
  dist_trx <- cophe_phyd[sps.in.cell, sps.in.cell]
  # MPD
  MPD_ <- mean(dist_trx[lower.tri(dist_trx)])
  return(MPD_)
}

# Load data

plants <- readRDS("Data/WorkData/plants.RDS")
phy <- readRDS("Data/WorkData/newphy.RDS")


# Calc richness

richness <- plants %>%
  group_by(cell) %>% 
  summarise(rich = length(unique(sps)))

rownames(richness) <- richness$cell

template_raster <- raster(xmn = -180, xmx = 180, 
                          ymn = -60, ymx = 85, 
                          resolution = .5, 
                          crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")

cells <- richness$cell
XY <- data.frame(xyFromCell(template_raster, as.numeric(cells)))
names(XY) <- c("x","y")

plot(rasterFromXYZ(data.frame(XY, richness$rich)))


# Calc PD

# Run with the Maximum Clade Credibility Tree
phy.mcc <- phangorn::maxCladeCred(phy)
# Compute the cophenetic distances for calculating MPD
cophe_tre <- cophenetic(phy.mcc)
# delete 100 trees for saving space
rm(phy)

##################
# Calc PD
ncores <- detectCores()
cl <- makeCluster(ncores)
clusterExport(cl, c("cells", "phy.mcc", "plants","calc_PD","keep.tip"))

PD <- pblapply(1:length(cells), function(i) calc_PD(i, plants, phy.mcc),
               cl = cl)

stopCluster(cl)
# Took 20min
PD <- unlist(PD)

saveRDS(PD, "Data/WorkData/PD.RDS")

##################
# Calc MPD
ncores <- detectCores()
cl <- makeCluster(ncores)
clusterExport(cl, c("cells", "cophe_tre", "plants","calc_MPD"))

MPD <- pblapply(1:length(cells), function(i) calc_MPD(i, plants, cophe_tre),
                cl = cl)

stopCluster(cl)
# Took 25min
MPD <- unlist(MPD)

saveRDS(MPD, "Data/WorkData/MPD.RDS") 

# no longer need the cophenetic tree. delete for saving space
rm(cophe_tre)

gc()


############ 
# Run FD in the Cluster

# Calc FD

# from skinny to occr
occr <- reshape2::dcast(plants, cell ~ sps, length)
rownames(occr) <- occr[,1]
occr <- occr[,-1]

# calc FDis

# Run PCA on traits to reduce dimensionality
traits.pc <- prcomp(traits, center = TRUE, scale = TRUE)
summary(traits.pc)
# library(devtools)
# install_github("vqv/ggbiplot")
# library(ggbiplot)
# ggbiplot(traits.pc, alpha = 0.001)
# select the 5 first PCs, explaining 88% of the complete variance in the trait dataset
traits.pc <- traits.pc$x[,1:5]

TraitDis <- dist(traits.pc[names(occr),]) # Species labels in the same order as in occr
fdis_prepared <- fdis_prep(TraitDis, as.matrix(occr))

ncores <- detectCores()
cl <- makeCluster(ncores)
clusterExport(cl, c("TraitDis", "cells", "plants","fdis_prepared","fdis_calc"))

FDis <- pblapply(1:length(cells), 
                 FUN = function(i){
                   # species that occur in the cell i
                   sps.in.cell <- plants$sps[which(plants$cell==cells[i])] 
                   if(length(sps.in.cell) > 5){
                     fdis_calc(fdis_prepared, sps.in.cell)$FDis
                   }
                   if(length(sps.in.cell) <= 5){
                     return(NA)
                   }
                 },
                 cl = cl)

stopCluster(cl)

FDis <- unlist(FDis)

saveRDS(FD, "FDis.RDS") 

# calc FRich

FRichall <- convhulln(traits, "FA")$vol

ncores <- detectCores()
cl <- makeCluster(ncores)
clusterExport(cl, c("traits.pc", "cells", "plants"))

FRich <- pblapply(1:length(cells), 
                 FUN = function(i){
                   # species that occur in the cell i
                   sps.in.cell <- plants$sps[which(plants$cell==cells[i])] 
                   if(length(sps.in.cell) > 5){
                     FRich_ <- convhulln(traits.pc[sps.in.cell,], "FA")$vol
                     return(FRich_/FRichall)
                   }
                   if(length(sps.in.cell) <= 5){
                     return(NA)
                   }
                 },
                 cl = cl)

stopCluster(cl)

FRich <- unlist(FRich)


# Calc Mean and SD traits /  cell

ncores <- detectCores()
cl <- makeCluster(ncores)
clusterExport(cl, c("cells", "plants","traitsmean"))

TMSD <- pblapply(1:length(cells), FUN = function(i){
  # species that occur in the cell i
  sps.in.cell <- plants$sps[which(plants$cell==cells[i])] 
  data.cell <- log1p(data.frame(traitsmean[sps.in.cell,])) # log-transform trait values
  
  SLA.m <- mean(data.cell$SLA)
  N.m <- mean(data.cell$N)
  P.m <- mean(data.cell$P)
  Longevity.m <- mean(data.cell$Longevity)
  SeedMass.m <- mean(data.cell$SeedMass)
  Height.m <- mean(data.cell$Height)
  RootDepth.m <- mean(data.cell$RootDept)
  WoodDens.m <- mean(data.cell$WoodDens)
  
  SLA.sd <- sd(data.cell$SLA)
  N.sd <- sd(data.cell$N)
  P.sd <- sd(data.cell$P)
  Longevity.sd <- sd(data.cell$Longevity)
  SeedMass.sd <- sd(data.cell$SeedMass)
  Height.sd <- sd(data.cell$Height)
  RootDepth.sd <- sd(data.cell$RootDept)
  WoodDens.sd <- sd(data.cell$WoodDens)
  
  return(data.frame(SLA.m,N.m,P.m,Longevity.m,SeedMass.m,Height.m,RootDepth.m,WoodDens.m,
                    SLA.sd,N.sd,P.sd,Longevity.sd,SeedMass.sd,Height.sd,RootDepth.sd,WoodDens.sd))
},
cl = cl)

stopCluster(cl)

TMSD <- do.call(rbind, TMSD)

saveRDS(TMSD, "TMSD.RDS") 

