library(targets)
library(tarchetypes)
library(future)
library(future.callr)
source("R/functions.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("ape", "bayestestR", "colorspace", "cowplot", "ggrepel", 
                            "ggtree", "phangorn", "phaseR", "phytools", "readxl", 
                            "rethinking", "rstan", "tidyverse"))
plan(callr)
# full workflow
list(
  # load phylogeny
  tar_target(filePhylo, "data/authority_97_20_09_21.trees", format = "file"),
  tar_target(phylo, read.nexus(filePhylo)),
  # load data
  tar_target(filePol, "data/political_authority_20_09_21.txt", format = "file"),
  tar_target(fileRel, "data/religious_authority_20_09_21.txt", format = "file"),
  tar_target(d, loadData(filePol, fileRel, phylo)),
  # load society names
  tar_target(fileSociety, "data/taxon_names.xlsx", format = "file"),
  tar_target(socNames, read_xlsx(fileSociety)),
  # get posterior trees for modelling
  tar_target(numTrees, 100),
  tar_target(iter, getPostSamples(numTrees)),
  # simulate OU stan model
  tar_target(fileSimStan, "stan/OU_sim_phy_ordinal.stan", format = "file"),
  tar_target(modelSimStan, stan_model(fileSimStan)),
  tar_target(dSim, simOUData(modelSimStan, phylo, iter),
             pattern = map(iter), iteration = "list"),
  tar_target(simModel, fitOUModel(modelStan, dSim, phylo, iter),
             pattern = map(dSim, iter), iteration = "list"),
  tar_target(simPost, extract.samples(simModel, pars = c("A", "Q")), 
             pattern = map(simModel), iteration = "list"),
  # fit OU stan model
  tar_target(fileStan, "stan/OU_fit_phy_ordinal.stan", format = "file"),
  tar_target(modelStan, stan_model(fileStan)),
  tar_target(ouModelInd, fitOUModel(modelStan, d, phylo, iter), 
             pattern = map(iter), iteration = "list"),
  tar_target(ouModel, sflist2stanfit(ouModelInd)),
  tar_target(post, extract.samples(ouModel)),
  # bayes factors for coevolution parameters
  tar_target(prior1, getPrior(1)),
  tar_target(prior2, getPrior(2)),
  tar_target(bf1, bayesfactor_parameters(posterior = post$A[,1,2], prior = prior1, direction = ">")),
  tar_target(bf2, bayesfactor_parameters(posterior = post$A[,2,1], prior = prior2, direction = ">")),
  tar_target(bf3, bayesfactor_parameters(posterior = post$A[,1,2] - post$A[,2,1], prior = prior1 - prior2)),
  # plot model results
  tar_target(deltaTheta, getEquilibriumChange(post)),
  tar_target(plotOU1, plotEquilibriumChange(deltaTheta)),
  tar_target(plotOU2, plotMedianTraitValues(d, post, phylo)),
  tar_target(plotOU3, plotFlowField(post)),
  tar_target(plotOU4, plotPredManifest(post)),
  tar_target(plotOU5, plotSelectionGradient(post)),
  tar_target(plotOU6, plotButterfly(phylo, iter, post, socNames)),
  tar_target(plotOU7, plotTrace(post, numTrees, numChains = 4)),
  tar_target(plotOU8, plotSim(simPost))
)