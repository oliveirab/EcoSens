library(spdep)
library(spatialreg)

setwd("/home/brolivei/SensProd/SARmodels")

# Models Rich
models2go <- readRDS("models2goRich.RDS")
datarun <- readRDS("datarun.RDS")
weightpar <- readRDS("weightpar.RDS")
nlw <- readRDS("nlw.RDS")

command_args <- commandArgs(trailingOnly = TRUE)

i <- as.numeric(command_args[1])

model.i <- spatialreg::errorsarlm(formula(models2go[i]),
                                  weights = datarun[,weightpar[i]],
                                  listw = nlw,
                                  data = datarun,
                                  zero.policy=TRUE)

saveRDS(model.i,
        paste("Results/","MD",i,".RDS",sep = ""))

# Models Disp
models2go <- readRDS("models2goDisp.RDS")
datarun <- readRDS("datarun.RDS")
weightpar <- readRDS("weightpar.RDS")
nlw <- readRDS("nlw.RDS")

command_args <- commandArgs(trailingOnly = TRUE)

i <- as.numeric(command_args[1])

model.i <- spatialreg::errorsarlm(formula(models2go[i]),
                                  weights = datarun[,weightpar[i]],
                                  listw = nlw,
                                  data = datarun,
                                  zero.policy=TRUE)

saveRDS(model.i,
        paste("Results/","MR",i,".RDS",sep = ""))
