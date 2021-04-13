
traits.2imput<- readRDS("traits_plants.RDS")

# 1) Load imputed files
my.list <- lapply(files.dir, function(x) readRDS(x)$anc_recon[traits.2imput$sps,])

# 2) Calc mean and sd trait values
traitsmean <- apply(simplify2array(my.list), 1:2, mean)
traitssd <- apply(simplify2array(my.list), 1:2, sd)

# 3) Save
saveRDS(traitsmean, "Results/traitsmean.RDS")
saveRDS(traitssd, "Results/traitssd.RDS")



