library(rethinking)
library(phytools)
source("R/functions.R")

trees <- read.nexus("data/authority_97_20_09_21.trees")
modelSimStan <- stan_model("stan/OU_sim_phy_ordinal.stan")

# Get parameter values from the prior
prior_sim <- simOUData(modelSimStan, phylo = trees, iter=1, output="parameters")

time_depth <- max(node.depth.edgelength(trees[[1]])) * 1000 # depth of phylogeny, which we'll need to scale our dt values

# Median parameter values
A <- apply(prior_sim$A, 2:3, median)
b <- apply(prior_sim$b, 2, median)
G <- t(chol(apply(prior_sim$Q, 2:3, median)))

# Organize
parms <- list(
  A = A,
  b = b,
  G = G
)

# Simulate co-evolution over a 500 year period, given some starting values for eta
sim_coev <- SDE_sim(init_eta=c(0,0), parms, time_depth=time_depth, time=500)

plot(NA, ylim=c(min(sim_coev), max(sim_coev)), xlim=c(0,nrow(sim_coev)), xlab="Time", ylab="eta")
lines(sim_coev[,1], col="darkred")
lines(sim_coev[,2], col="skyblue")


