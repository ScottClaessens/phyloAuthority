library(targets)
library(tarchetypes)
library(future)
library(future.callr)
source("R/functions.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("ape", "bayestestR", "brms", "colorspace", "cowplot", 
                            "geosphere", "ggrepel", "ggtree", "HDInterval", "phangorn", 
                            "phaseR", "phytools", "readxl", "rethinking", "rstan", 
                            "tidybayes", "tidyverse"))
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
  # load longitude and latitude points
  tar_target(fileLonLat, "data/authority_coordinates.txt", format = "file"),
  tar_target(lonLat, read_tsv(fileLonLat, col_types = "cdd")),
  # plot linguistic vs. geographic distance
  tar_target(plotDist, plotDistance(d, phylo, lonLat)),
  # get posterior trees for modelling
  tar_target(numTrees, 100),
  tar_target(iter, getPostSamples(numTrees)),
  # fit phylogenetic GLMM in brms
  tar_target(modelPhyloGLMM, getPhyloGLMMinitial(d, phylo)),
  tar_target(phyloGLMMind, fitPhyloGLMM(modelPhyloGLMM, phylo, iter),
             pattern = map(iter), iteration = "list"),
  tar_target(phyloGLMM, combine_models(mlist = phyloGLMMind, check_data = FALSE)),
  tar_target(postPhyloGLMM, as_draws_array(phyloGLMM, variable = "^cor_", regex = TRUE)),
  tar_target(plotPhyloCor, plotPhyloGLMM(postPhyloGLMM)),
  # simulate OU stan model
  tar_target(fileSimStan, "stan/OU_sim_phy_ordinal.stan", format = "file"),
  tar_target(modelSimStan, stan_model(fileSimStan)),
  tar_target(dSim, simOUData(modelSimStan, phylo, iter),
             pattern = map(iter), iteration = "list"),
  tar_target(simModel, fitOUModel(modelStan, dSim, phylo, iter),
             pattern = map(dSim, iter), iteration = "list"),
  tar_target(simPost, extract.samples(simModel, pars = c("A", "Q")), 
             pattern = map(simModel), iteration = "list"),
  tar_target(plotSimulation, plotSim(simPost)),
  # fit OU stan model
  tar_target(fileStan, "stan/OU_fit_phy_ordinal.stan", format = "file"),
  tar_target(modelStan, stan_model(fileStan)),
  tar_target(ouModelInd, fitOUModel(modelStan, d, phylo, iter), 
             pattern = map(iter), iteration = "list"),
  tar_target(ouModel, sflist2stanfit(ouModelInd)),
  tar_target(post, extract.samples(ouModel)),
  # fit OU stan model with geographic control
  tar_target(fileStanGeo, "stan/OU_fit_phy_ordinal_geographic.stan", format = "file"),
  tar_target(modelStanGeo, stan_model(fileStanGeo)),
  tar_target(ouModelGeoInd, fitOUModelGeo(modelStanGeo, d, phylo, iter, lonLat),
             pattern = map(iter), iteration = "list"),
  tar_target(ouModelGeo, sflist2stanfit(ouModelGeoInd)),
  tar_target(postGeo, extract.samples(ouModelGeo)),
  # bayes factors for coevolution parameters
  tar_target(prior1, getPrior(1)),
  tar_target(prior2, getPrior(2)),
  tar_target(bf1, bayesfactor_parameters(posterior = post$A[,1,2], prior = prior1, direction = ">")),
  tar_target(bf2, bayesfactor_parameters(posterior = post$A[,2,1], prior = prior2, direction = ">")),
  tar_target(bf3, bayesfactor_parameters(posterior = post$A[,1,2] - post$A[,2,1], prior = prior1 - prior2)),
  tar_target(bf1Geo, bayesfactor_parameters(posterior = postGeo$A[,1,2], prior = prior1, direction = ">")),
  tar_target(bf2Geo, bayesfactor_parameters(posterior = postGeo$A[,2,1], prior = prior2, direction = ">")),
  tar_target(bf3Geo, bayesfactor_parameters(posterior = postGeo$A[,1,2] - postGeo$A[,2,1], prior = prior1 - prior2)),
  # plot model results
  tar_target(deltaTheta, getEquilibriumChange(post)),
  tar_target(plotOU1, plotEquilibriumChange(deltaTheta)),
  tar_target(plotOU2, plotMedianTraitValues(d, post, phylo)),
  tar_target(plotOU3, plotFlowField(post)),
  tar_target(plotOU4, plotPredManifest(post)),
  tar_target(plotOU5, plotSelectionGradient(post)),
  tar_target(plotOU6, plotButterfly(phylo, iter, post, socNames, d)),
  tar_target(plotOU7, plotTrace(post, numTrees, numChains = 4)),
  # plot model results with geographic control
  tar_target(deltaThetaGeo, getEquilibriumChange(postGeo)),
  tar_target(plotOUGeo1, plotEquilibriumChange(deltaThetaGeo, geo = TRUE)),
  tar_target(plotOUGeo2, plotMedianTraitValues(d, postGeo, phylo, geo = TRUE)),
  tar_target(plotOUGeo3, plotFlowField(postGeo, geo = TRUE)),
  tar_target(plotOUGeo4, plotPredManifest(postGeo, geo = TRUE)),
  tar_target(plotOUGeo5, plotSelectionGradient(postGeo, geo = TRUE)),
  tar_target(plotOUGeo6, plotButterfly(phylo, iter, postGeo, socNames, d, geo = TRUE)),
  tar_target(plotOUGeo7, plotTrace(postGeo, numTrees, numChains = 4, geo = TRUE))
)
