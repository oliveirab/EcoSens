
# Load packages

rm(list=ls())
gc()

list.of.packages <- c("raster","dplyr","pbapply","parallel","doParallel","reshape","FD","reshape2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

# Source code

source("fdisp_loop.R")

# Calc FDis function
calc_FDis <- function(i, sp.names, sps.in.cell){
  sps.in.cell <- sp.names[cells.sp[[i]]]
  if(length(sps.in.cell) > 5){
    fdis_calc(fdis_prepared, sps.in.cell)$FDis
  }
  else{
    NA
  }
}

# Calc FRic function
# species that occur in the cell i
calc_FRich <- function(i, sp.names, sps.in.cell){
  sps.in.cell <- sp.names[cells.sp[[i]]]
  if(length(sps.in.cell) > 5){
    FRich_ <- convhulln(traits.pc[sps.in.cell,], "FA")$vol
    return(FRich_/FRichall)
  }
  else{
    NA
  }
}

# Load data

plants <- readRDS("Data/WorkData/plants.RDS")
traits <- readRDS("Data/WorkData/traitsmean.RDS")

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

############ 
# Run FD in the Cluster

# Calc FD

# from skinny to occr
# occr <- reshape2::dcast(plants, cell ~ sps, length)
# rownames(occr) <- occr[,1]
# occr <- occr[,-1]
# 
# saveRDS(occr, "Data/WorkData/occr.RDS")

occr <- readRDS("Data/WorkData/occr.RDS")

cells.sp <- pbapply(occr, 1, function(x) which(x==1))
cells.sp <- lapply(cells.sp, as.numeric)

sp.names <- names(occr)

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

# fdis_prepared <- fdis_prep(TraitDis, as.matrix(occr))
# saveRDS(fdis_prepared, "Data/WorkData/fdis_prepared.RDS")
fdis_prepared <- readRDS("Data/WorkData/fdis_prepared.RDS")

ncores <- detectCores()
cl <- makeCluster(ncores)
clusterExport(cl, c("cells.sp","fdis_prepared","fdis_calc","sp.names","calc_FDis"))

FDis <- pbsapply(1:nrow(occr), 
                 FUN = function(i){
                   calc_FDis(i, sp.names, sps.in.cell)
                 },
                 cl = cl)

stopCluster(cl)

saveRDS(FDis, "Data/WorkData/FDis.RDS") 

# calc FRich

FRichall <- convhulln(traits, "FA")$vol

ncores <- detectCores()
cl <- makeCluster(ncores)
clusterExport(cl, c("traits.pc", "sp.names","cells.sp", "convhulln", "FRichall", "calc_FRich"))

FRich <- pbsapply(1:nrow(occr), 
                  FUN = function(i){
                    calc_FRich(i, sp.names, sps.in.cell)
                  },
                  cl = cl)
stopCluster(cl)

saveRDS(FRich, "Data/WorkData/FRich.RDS") 

# Calc Mean and SD traits /  cell

ncores <- detectCores()
cl <- makeCluster(ncores)
clusterExport(cl, c("sp.names", "cells.sp","traits"))

TMSD <- pblapply(1:length(cells), FUN = function(i){
  # species that occur in the cell i
  sps.in.cell <- sp.names[cells.sp[[i]]]
  data.cell <- data.frame(exp(traits[sps.in.cell,])) # exp for removing log
  
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

saveRDS(TMSD, "Data/WorkData/TMSD.RDS") 

# test
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

plot(rasterFromXYZ(data.frame(XY, FDis)))
plot(rasterFromXYZ(data.frame(XY, FRich)))

plot(rasterFromXYZ(data.frame(XY, TMSD$WoodDens.m)))

