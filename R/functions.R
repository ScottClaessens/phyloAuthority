# custom functions

# load data
loadData <- function(filePol, fileRel, phylo) {
  read_tsv(filePol, col_types = "cd") %>%
    left_join(read_tsv(fileRel, col_types = "cd"), by = "Language") %>%
    # positive integers
    transmute(
      language = Language,
      polAuth  = `Political Authority` + 1,
      relAuth  = `Religious Authority` + 1
    ) %>%
    # sort data to match phylogeny tips
    arrange(match(language, phylo[[1]]$tip.label))
}

# get trees to sample
getPostSamples <- function(numTrees) {
  set.seed(2113)
  out <- sample(1:1000, numTrees)
  return(out)
}

# fit OU stan model
fitOUModel <- function(modelStan, d, phylo, iter) {
  # choose tree
  tree <- phylo[[iter]]
  # cut up tree into segments
  times <- node.depth.edgelength(tree)
  # line up date of each node with the split points in the tree
  split_points <- sort(unique(times))
  node_time <- match(times, split_points)
  # create a sequence of nodes, respecting temporal order
  node_seq <- seq(from = 1, to = length(node_time))
  node_seq <- node_seq[order(node_time)]
  # find the "parent" node for each node and tip of the tree
  parent <- Ancestors(tree, node_seq, type = "parent")
  # parent time indicates amount of time since the parent node
  # scaled by the total depth of the tree
  parent_time <- rep(NA, length(node_seq))
  parent_time[1] <- -99 # placeholder for ancestral state
  for (i in 2:length(parent_time)) {
    parent_time[i] <- 
      (node.depth.edgelength(tree)[node_seq[i]] - node.depth.edgelength(tree)[parent[i]]) / 
      max(node.depth.edgelength(tree))
  }
  # get total num segments in the tree
  N_seg <- length(node_seq)
  # data list for stan
  data_list <- list(
    N_seg = N_seg,
    node_seq = node_seq,
    parent = parent,
    ts = parent_time,
    N = length(tree$tip.label),
    J = 2,
    y = cbind(d$polAuth, d$relAuth)
  )
  # fit model
  fit <- sampling(modelStan, data = data_list, chains = 1,
                  cores = 1, iter = 4000, init = "0",
                  control = list(adapt_delta = 0.999),
                  seed = 2113)
  return(fit)
}

# calculate delta theta for plot
delta_theta <- function(polAuth, relAuth, A, b, resp) {
  # political authority
  med_polAuth  <- median(polAuth)
  diff_polAuth <- mad(polAuth) 
  # religious authority
  med_relAuth  <- median(relAuth)
  diff_relAuth <- mad(relAuth)
  # calculate delta theta
  if (resp == "polAuth") {
    med_theta <- -(med_relAuth*A[1,2] + b[1])/A[1,1]
    diff_theta <- -((med_relAuth + diff_relAuth)*A[1,2] + b[1])/A[1,1]
    return( (diff_theta - med_theta) / diff_polAuth )
  } else if (resp == "relAuth") {
    med_theta <- -(med_polAuth*A[2,1] + b[2])/A[2,2]
    diff_theta <- -((med_polAuth + diff_polAuth)*A[2,1] + b[2])/A[2,2]
    return( (diff_theta - med_theta) / diff_relAuth )
  }
}

# get change in equilibrium trait value
getEquilibriumChange <- function(post) {
  # initialise vectors
  delta_theta_polAuth <- c()
  delta_theta_relAuth <- c()
  # calculate delta theta
  n_samps <- length(post$A[,1,1])
  for (i in 1:n_samps) {
    delta_theta_polAuth[i] <- delta_theta(polAuth = post$eta[i,1:97,1], 
                                          relAuth = post$eta[i,1:97,2], 
                                          A = post$A[i,,], b = post$b[i,], 
                                          resp = "polAuth")
    delta_theta_relAuth[i] <- delta_theta(polAuth = post$eta[i,1:97,1], 
                                          relAuth = post$eta[i,1:97,2], 
                                          A = post$A[i,,], b = post$b[i,], 
                                          resp = "relAuth")
  }
  # return data frame
  out <- data.frame(
    delta_theta = c(delta_theta_polAuth, delta_theta_relAuth),
    resp = rep(c("polAuth", "relAuth"), each = n_samps)
  )
  return(out)
}

# plot change in equilibrium trait value
plotEquilibriumChange <- function(df, geo = FALSE) {
  # get facet labels
  dat_text <- tibble(
    label1 = c(expression(Rel %->% Pol), 
               expression(Pol %->% Rel)),
    label2 = c(paste0("PP = ", round(mean(df$delta_theta[df$resp=="polAuth"] > 0), 2)),
               paste0("PP = ", round(mean(df$delta_theta[df$resp=="relAuth"] > 0), 2))),
    resp   = c("polAuth", "relAuth")
  )
  # plot
  out <-
    ggplot(df, aes(x = delta_theta)) + 
    geom_density(aes(fill = resp, colour = resp), alpha = 0.5) + 
    geom_vline(xintercept = 0, colour = "indianred", linetype = 'dashed', lwd = 1) +
    facet_wrap(resp ~ ., nrow = 2) +
    geom_text(data = dat_text, aes(x = 7.5, y = 0.15, label = label1, colour = resp), 
              size = 4.5, parse = TRUE) +
    geom_text(data = dat_text, aes(x = 7.5, y = 0.12, label = label2, colour = resp), 
              size = 4.5) +
    scale_fill_manual(values = c("#5387b6","#c55852")) +
    scale_colour_manual(values = c("#5387b6","#c55852")) +
    theme_bw(base_size = 14) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none") +
    scale_x_continuous(name = expression(paste(Delta, theta)["z"]), 
                       limits = c(-5, 15)) +
    ylab(NULL)
  # save plot
  ggsave(out, filename = paste0("figures/ouModel", ifelse(geo, "WithGeographicControl", ""), "/plotOU.pdf"), width = 5, height = 5)
  return(out)
}

# scatter plot of posterior median trait values for sample societies
plotMedianTraitValues <- function(d, post, phylo, geo = FALSE) {
  # extract median eta values for sample societies
  polAuth_median <- apply(post$eta[,1:97,1], 2, median)
  relAuth_median <- apply(post$eta[,1:97,2], 2, median)
  # plot
  out <-
    tibble(
      language = phylo[[1]]$tip.label,
      polAuth = scale(polAuth_median),
      relAuth = scale(relAuth_median)
      ) %>%
    ggplot(aes(x = polAuth, y = relAuth, label = language)) +
    geom_point() +
    geom_text_repel(size = 2) +
    theme_classic() +
    labs(y = "Religious authority (z-score)",
         x = "Political authority (z-score)")
  # save
  ggsave(out, filename = paste0("figures/ouModel", ifelse(geo, "WithGeographicControl", ""), "/plotMedTraitVals.pdf"), width = 5, height = 5)
  return(out)
}

# flow field diagram of evolution
plotFlowField <- function(post, geo = FALSE) {
  # get distribution of latent variables among study societes
  median_polAuth <- median(apply(post$eta[,1:97,1], 2, median))
  median_relAuth <- median(apply(post$eta[,1:97,2], 2, median))
  mad_polAuth    <- mad(apply(post$eta[,1:97,1], 2, median))
  mad_relAuth    <- mad(apply(post$eta[,1:97,2], 2, median))
  low_polAuth    <- median_polAuth - mad_polAuth*2.5 # -2.5 SD
  high_polAuth   <- median_polAuth + mad_polAuth*2.5 # +2.5 SD
  low_relAuth    <- median_relAuth - mad_relAuth*2.5 # -2.5 SD
  high_relAuth   <- median_relAuth + mad_relAuth*2.5 # +2.5 SD
  # get parameter values
  A <- apply(post$A, 2:3, median)
  b <- apply(post$b, 2, median)
  sigma <- apply(post$sigma, 2, median)
  # function for flow field diagram
  OU <- function(t, y, parameters) {
    dy <- numeric(2)
    dy[1] <- y[1]*A[1,1] + y[2]*A[1,2] + b[1]
    dy[2] <- y[2]*A[2,2] + y[1]*A[2,1] + b[2]
    list(dy)
  }
  # set up canvas
  pdf(paste0("figures/ouModel", ifelse(geo, "WithGeographicControl", ""), "/plotPhasePlane.pdf"), 
      width = 6, height = 6, pointsize = 12)
  par(pty = "s")
  # create flow field diagram
  OU.flowField <- 
    flowField(
      OU, 
      xlim = c(low_polAuth-0.2, high_polAuth+0.2), 
      ylim = c(low_relAuth-0.2, high_relAuth+0.2), 
      parameters = NA, add = FALSE, 
      xlab = "", ylab = "", 
      points = 12, col = "grey", 
      xaxt = 'n', yaxt = 'n', 
      arrow.type = "proportional", 
      frac = 1.5, xaxs = "i", yaxs = "i", 
      axes = FALSE, lwd = 2
    )
  mtext(side = 1, "Political authority (z-score)", at = median_polAuth, 
        line = 2.5, cex = 1.3)
  mtext(side = 2, "Religious authority (z-score)", at = median_relAuth, 
        line = 2.5, cex = 1.3)
  # add cutpoint lines
  abline(v = median(post$c1[,1]), lty = "dashed", lwd = 1.3, col = "#5387b6")
  abline(v = median(post$c1[,2]), lty = "dashed", lwd = 1.3, col = "#5387b6")
  abline(v = median(post$c1[,3]), lty = "dashed", lwd = 1.3, col = "#5387b6")
  abline(h = median(post$c2[,1]), lty = "dashed", lwd = 1.3, col = "#c55852")
  abline(h = median(post$c2[,2]), lty = "dashed", lwd = 1.3, col = "#c55852")
  abline(h = median(post$c2[,3]), lty = "dashed", lwd = 1.3, col = "#c55852")
  # add groupings
  mtext(side = 3, "Absent"    , at = -2.50, line = -0.1, cex = 0.75, col = "#5387b6")
  mtext(side = 3, "Sublocal"  , at = -0.95, line = -0.1, cex = 0.75, col = "#5387b6")
  mtext(side = 3, "Local"     , at =  0.50, line = -0.1, cex = 0.75, col = "#5387b6")
  mtext(side = 3, "Supralocal", at =  4.20, line = -0.1, cex = 0.75, col = "#5387b6")
  mtext(side = 4, "Absent"    , at = -2.40, line = -0.1, cex = 0.75, col = "#c55852", las = 1)
  mtext(side = 4, "Sublocal"  , at = -1.20, line = -0.1, cex = 0.75, col = "#c55852", las = 1)
  mtext(side = 4, "Local"     , at =  0.35, line = -0.1, cex = 0.75, col = "#c55852", las = 1)
  mtext(side = 4, "Supralocal", at =  4.20, line = -0.1, cex = 0.75, col = "#c55852", las = 1)
  # add nullclines to phase plane
  #nc <- 
  #  nullclines(
  #    OU, 
  #    xlim = c(low_polAuth, high_polAuth),
  #    ylim = c(low_relAuth, high_relAuth), 
  #    parameters = NA, points = 20, 
  #    axes = FALSE, col = c("#c55852","#5387b6"),
  #    add.legend = FALSE, lwd = 4)
  # add axes
  axis(1, at = c(low_polAuth, median_polAuth, high_polAuth),
       labels = (c(low_polAuth, median_polAuth, high_polAuth) - median_polAuth) / mad_polAuth)
  axis(2, at = c(low_relAuth, median_relAuth, high_relAuth),
       labels = (c(low_relAuth, median_relAuth, high_relAuth) - median_relAuth) / mad_relAuth)
  dev.off()
}

# make predictions for difference in trait values
plotSelectionGradient <- function(post, geo = FALSE) {
  # get distribution of latent variables among study societes
  median_polAuth <- median(apply(post$eta[,1:97,1], 2, median))
  median_relAuth <- median(apply(post$eta[,1:97,2], 2, median))
  mad_polAuth    <- mad(apply(post$eta[,1:97,1], 2, median))
  mad_relAuth    <- mad(apply(post$eta[,1:97,2], 2, median))
  low_polAuth    <- median_polAuth - mad_polAuth*2.5 # -2.5 SD
  high_polAuth   <- median_polAuth + mad_polAuth*2.5 # +2.5 SD
  low_relAuth    <- median_relAuth - mad_relAuth*2.5 # -2.5 SD
  high_relAuth   <- median_relAuth + mad_relAuth*2.5 # +2.5 SD
  # get parameter values
  A <- apply(post$A, 2:3, median)
  b <- apply(post$b, 2, median)
  sigma <- apply(post$sigma, 2, median)
  # ou sde function
  OU_sde <- function(polAuth, relAuth, resp) {
    if (resp == "polAuth") {
      delta_sigma <- ((A[1,1]*polAuth + A[1,2]*relAuth + b[1])/mad_polAuth) / 
        (sigma[1]/mad_polAuth)
    } else if (resp == "relAuth") {
      delta_sigma <- ((A[2,2]*relAuth + A[2,1]*polAuth + b[2])/mad_relAuth) / 
        (sigma[2]/mad_relAuth)
    }
    return(delta_sigma)
  }
  # get predictions for different levels of trait
  polAuth_seq <- seq(from = low_polAuth, to = high_polAuth, length.out = 20)
  relAuth_seq <- seq(from = low_relAuth, to = high_relAuth, length.out = 20)
  preds <- expand.grid(x = polAuth_seq, y = relAuth_seq)
  preds$polAuth <- NA
  preds$relAuth <- NA
  for (i in 1:nrow(preds)) {
    preds$polAuth[i] <- OU_sde(polAuth = preds$x[i], relAuth = preds$y[i], 
                               resp = "polAuth")
    preds$relAuth[i] <- OU_sde(polAuth = preds$x[i], relAuth = preds$y[i], 
                               resp = "relAuth")
  }
  out <- 
    preds %>% 
    pivot_longer(-c(x, y)) %>%
    mutate(name = factor(name, 
                         labels = c(
                           expression(paste(Delta,"Political authority") ), 
                           expression(paste(Delta,"Religious authority") )
                           ))) %>%
    ggplot(aes(x = (x - median_polAuth) / mad_polAuth, 
               y = (y - median_relAuth) / mad_relAuth,
               z = value, fill = value)) +
    facet_wrap(~ name, labeller = label_parsed) +
    geom_raster() +
    geom_contour(colour = "white", breaks = c(-1, 1)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_continuous_divergingx(palette = "Geyser", trans = "reverse", 
                                     guide = guide_legend(reverse = TRUE)) + 
    theme_bw(base_size = 14) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1.5, "lines"),
          strip.background = element_blank()) + 
    labs(fill = expression(frac(paste(Delta, alpha), sigma))) +
    xlab("Political authority (z-score)") +
    ylab("Religious authority (z-score)") + 
    coord_fixed()
  # save plot
  ggsave(filename = paste0("figures/ouModel", ifelse(geo, "WithGeographicControl", ""), "/plotSelectionGradient.pdf"), 
         width = 6, height = 3.5, dpi = 600)
  return(out)
}

# plot predictions for manifest variables
plotPredManifest <- function(post, geo = FALSE) {
  # get distribution of latent variables among study societes
  median_polAuth <- median(apply(post$eta[,1:97,1], 2, median))
  median_relAuth <- median(apply(post$eta[,1:97,2], 2, median))
  mad_polAuth    <- mad(apply(post$eta[,1:97,1], 2, median))
  mad_relAuth    <- mad(apply(post$eta[,1:97,2], 2, median))
  low_polAuth    <- median_polAuth - mad_polAuth*2.5 # -2.5 SD
  high_polAuth   <- median_polAuth + mad_polAuth*2.5 # +2.5 SD
  low_relAuth    <- median_relAuth - mad_relAuth*2.5 # -2.5 SD
  high_relAuth   <- median_relAuth + mad_relAuth*2.5 # +2.5 SD
  # sequences for prediction
  polAuth_seq <- seq(low_polAuth, high_polAuth, length.out = 30)
  relAuth_seq <- seq(low_relAuth, high_relAuth, length.out = 30)
  # plotting function
  plotFun <- function(cuts, seq, equilibrium = FALSE, self, other, ci, med, mad, xlab, ylab) {
    tibble(
      c1 = cuts[,1],
      c2 = cuts[,2],
      c3 = cuts[,3],
      Aself = post$A[,self,self],
      Across = post$A[,self,other],
      b = post$b[,self]
    ) %>%
      expand_grid(auth = seq) %>%
      # if retrieving equilbrium trait values for other variable, calculate here
      mutate(auth2 = if (equilibrium) -(auth*Across + b)/Aself else auth) %>%
      # on manifest scale
      mutate(Absent     = inv_logit(c1 - auth2),
             Sublocal   = inv_logit(c2 - auth2) - inv_logit(c1 - auth2),
             Local      = inv_logit(c3 - auth2) - inv_logit(c2 - auth2),
             Supralocal = 1 - inv_logit(c3 - auth2),
      ) %>%
      pivot_longer(cols = Absent:Supralocal, 
                   names_to = "Authority level", 
                   values_to = "pred") %>%
      group_by(auth, `Authority level`) %>%
      summarise(median = median(pred),
                lower  = quantile(pred, (1 - ci) / 2),
                upper  = quantile(pred, (1 + ci) / 2)) %>%
      mutate(`Authority level` = factor(`Authority level`, 
                                        levels = c("Absent", "Sublocal", 
                                                   "Local", "Supralocal"))) %>%
      ggplot(aes(x = auth, y = median, ymin = lower, ymax = upper)) +
      geom_ribbon(aes(fill = `Authority level`), alpha = 0.2) +
      geom_line(aes(colour = `Authority level`), size = 1) +
      scale_x_continuous(name = xlab, limits = (c(-2.5, 2.5)*mad) + med,
                         breaks = (c(-2.5, 0, 2.5)*mad) + med,
                         labels = function(x) (x - med) / mad) +
      scale_fill_manual(values = c("#999999", "#F0E442", "#E69F00", "#D55E00")) +
      scale_colour_manual(values = c("#999999", "#F0E442", "#E69F00", "#D55E00")) +
      ylab(ylab) +
      theme_classic()
  }
  # plotting
  # political authority (at different levels of pol z-score)
  pA <- plotFun(post$c1, polAuth_seq, equilibrium = FALSE,
                self = 1, other = 2, ci = 0.5, median_polAuth, mad_polAuth, 
                "Political authority (z-score)", 
                "Probability of political\nauthority level")
  # religious authority (at different levels of rel z-score)
  pB <- plotFun(post$c2, relAuth_seq, equilibrium = FALSE,
                self = 2, other = 1, ci = 0.5, median_relAuth, mad_relAuth, 
                "Religious authority (z-score)", 
                "Probability of religious\nauthority level")
  # equilibrium political authority (at different levels of rel z-score)
  pC <- plotFun(post$c1, relAuth_seq, equilibrium = TRUE,
                self = 1, other = 2, ci = 0.5, median_relAuth, mad_relAuth, 
                "Religious authority (z-score)", 
                "Probability of political authority\nlevel at equilibrium")
  # equilibrium religious authority (at different levels of pol z-score)
  pD <- plotFun(post$c2, polAuth_seq, equilibrium = TRUE,
                self = 2, other = 1, ci = 0.5, median_polAuth, mad_polAuth, 
                "Political authority (z-score)", 
                "Probability of religious authority\nlevel at equilibrium")
  # put together
  out <- plot_grid(pA + theme(legend.position = "none"), 
                   pC + theme(legend.position = "none"), 
                   pB + theme(legend.position = "none"), 
                   pD + theme(legend.position = "none"), 
                   nrow = 2)
  out <- plot_grid(out, get_legend(pA), nrow = 1, rel_widths = c(1, 0.25))
  # save
  ggsave(out, filename = paste0("figures/ouModel", ifelse(geo, "WithGeographicControl", ""), "/plotPred.pdf"), 
         height = 6, width = 7.5)
  return(out)
}

# create butterfly plot with z-scores
plotButterfly <- function(phylo, iter, post, socNames, d, geo = FALSE) {
  # to save memory, retain only parameters needed
  post <- post[c("c1","c2","eta")]
  # create ultrametric maximum clade credibility tree
  cons <- phylo %>% mcc() %>% force.ultrametric()
  # get posterior eta value from particular model
  getPostEta <- function(i, sampNode, var) {
    # length of posterior samples from particular model
    len <- dim(post$eta)[1] / length(iter)
    # if node isn't included in particular model, return NAs
    if (is.na(sampNode)) {
      return(rep(NA, len))
    # else return posterior samples
    } else {
      # get particular model number in sequence
      modelNum <- which(iter == i)
      # posterior samples to extract for particular model
      sampStart <- (len * (modelNum - 1)) + 1
      sampEnd   <- (len * modelNum)
      return(post$eta[sampStart:sampEnd,sampNode,var])
    }
  }
  # put together data
  p <-
    tibble(iter = iter) %>%
    mutate(
      # get consensus tree nodes (including tips)
      node = purrr::map(iter, function(x) c(1:97, matchNodes(cons, phylo[[x]])[,1])),
      # get all nodes in sampled tree that match consensus tree (including tips)
      sampNode = purrr::map(iter, function(x) c(1:97, matchNodes(cons, phylo[[x]])[,2]))
      ) %>%
    unnest(c(node, sampNode)) %>%
    arrange(node) %>%
    # get trait values at node from posteriors of different models
    mutate(polAuth = map2(iter, sampNode, getPostEta, var = 1),
           relAuth = map2(iter, sampNode, getPostEta, var = 2)) %>%
    unnest(c(polAuth, relAuth))
  # save memory again
  post <- post[c("c1","c2")]
  # calculate probs
  p$pol_probAbsent     <- inv_logit(post$c1[,1] - p$polAuth)
  p$pol_probSublocal   <- inv_logit(post$c1[,2] - p$polAuth) - p$pol_probAbsent
  p$pol_probLocal      <- inv_logit(post$c1[,3] - p$polAuth) - p$pol_probSublocal
  p$pol_probSupralocal <- 1 - p$pol_probLocal
  p$rel_probAbsent     <- inv_logit(post$c2[,1] - p$relAuth)
  p$rel_probSublocal   <- inv_logit(post$c2[,2] - p$relAuth) - p$rel_probAbsent
  p$rel_probLocal      <- inv_logit(post$c2[,3] - p$relAuth) - p$rel_probSublocal
  p$rel_probSupralocal <- 1 - p$rel_probLocal
  # density plot function
  plotDens <- function(p, Node, var) {
    # get fewer random samples to plot
    set.seed(1)
    samps <- sample(1:2e+05, 4000)
    # plot
    p %>%
      group_by(node) %>%
      slice(samps) %>%
      ungroup() %>%
      filter(node == Node) %>%
      pivot_longer(cols = starts_with(paste0(var, "_")),
                   names_to = "auth") %>%
      mutate(auth = ifelse(str_detect(auth, "Absent"), "Absent",
                           ifelse(str_detect(auth, "Sublocal"), "Sublocal",
                                  ifelse(str_detect(auth, "Local"), "Local", "Supralocal"))),
             auth = factor(auth, levels = c("Absent", "Sublocal", "Local", "Supralocal"))) %>%
      ggplot(aes(x = value, fill = auth)) +
      geom_density(alpha = 0.7) +
      scale_fill_manual(name = "Authority level",
                        values = c(c("#e9fafa", "#F0E442", "#E69F00", "#D55E00"))) +
      scale_x_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.5, 1)) +
      theme_classic() +
      theme(legend.position = "bottom",
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 8),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank())
  }
  # create density plots
  polDens1 <- plotDens(p, Node = 98 , var = "pol")
  polDens2 <- plotDens(p, Node = 102, var = "pol")
  polDens3 <- plotDens(p, Node = 126, var = "pol")
  polDens4 <- plotDens(p, Node = 130, var = "pol")
  polDens5 <- plotDens(p, Node = 132, var = "pol")
  relDens1 <- plotDens(p, Node = 98 , var = "rel")
  relDens2 <- plotDens(p, Node = 102, var = "rel")
  relDens3 <- plotDens(p, Node = 126, var = "rel")
  relDens4 <- plotDens(p, Node = 130, var = "rel")
  relDens5 <- plotDens(p, Node = 132, var = "rel")
  # get mean trait values for trees
  p <-
    p %>%
    group_by(node) %>%
    summarise(polAuth = median(polAuth, na.rm = TRUE),
              relAuth = median(relAuth, na.rm = TRUE)) %>%
    # match to data for plot
    mutate(language = d$language[node]) %>%
    left_join(d, by = "language") %>%
    mutate(polAuth.y = ifelse(polAuth.y == 1, "Absent",
                              ifelse(polAuth.y == 2, "Sublocal",
                                     ifelse(polAuth.y == 3, "Local", "Supralocal"))),
           polAuth.y = factor(polAuth.y, levels = c("Absent", "Sublocal", "Local", "Supralocal")),
           relAuth.y = ifelse(relAuth.y == 1, "Absent",
                              ifelse(relAuth.y == 2, "Sublocal",
                                     ifelse(relAuth.y == 3, "Local", "Supralocal"))),
           relAuth.y = factor(relAuth.y, levels = c("Absent", "Sublocal", "Local", "Supralocal")))
  # get new society names for plot
  cons$tip.label <- socNames$Society[match(cons$tip.label, socNames$Language)]
  # create plots
  pA <- 
    ggtree(cons, aes(colour = relAuth.x), right = TRUE, size = 0.4) %<+% p + 
    geom_tiplab(colour = "black", alpha = 1, size = 1.5, hjust = 0.5, offset = 1) +
    scale_colour_continuous(name = "Religious\nauthority",
                            low = "gray99", high = "black", 
                            breaks = -2:3, limits = c(-2.5, 4),
                            guide = guide_colourbar(ticks = FALSE,
                                                    label.position = "right")) +
    ggnewscale::new_scale_color() +
    geom_tippoint(aes(colour = relAuth.y), show.legend = FALSE) +
    scale_colour_manual(values = c("#e9fafa", "#F0E442", "#E69F00", "#D55E00")) +
    theme(legend.position = c(0.1, 0.08),
          legend.title = element_text(hjust = 0.35),
          plot.margin = margin(5, 0, 5, 5, unit = "mm")) +
    xlim(c(0, 6.5))
  pB <- 
    ggtree(cons, aes(colour = polAuth.x), right = TRUE, size = 0.4) %<+% p + 
    scale_x_reverse() +
    scale_colour_continuous(name = "Political\nauthority",
                            low = "gray99", high = "black", 
                            breaks = -2:3, limits = c(-2.5, 4),
                            guide = guide_colourbar(ticks = FALSE,
                                                    label.position = "left")) +
    ggnewscale::new_scale_color() +
    geom_tippoint(aes(colour = polAuth.y), show.legend = FALSE) +
    scale_colour_manual(values = c("#e9fafa", "#F0E442", "#E69F00", "#D55E00")) +
    theme(legend.position = c(0.9, 0.08),
          legend.title = element_text(hjust = 0.5),
          plot.margin = margin(5, 5, 5, 0, unit = "mm"))
  # put together
  leftCol <- plot_grid(relDens1 + theme(legend.position = "none"), 
                       relDens2 + theme(legend.position = "none"), 
                       relDens3 + theme(legend.position = "none"), 
                       relDens4 + theme(legend.position = "none"), 
                       relDens5 + theme(legend.position = "none"), 
                       ncol = 1, labels = letters[1:5], label_size = 8,
                       label_x = 0.1)
  rightCol <- plot_grid(polDens1 + theme(legend.position = "none"), 
                        polDens2 + theme(legend.position = "none"), 
                        polDens3 + theme(legend.position = "none"), 
                        polDens4 + theme(legend.position = "none"), 
                        polDens5 + theme(legend.position = "none"), 
                        ncol = 1, labels = letters[6:10], label_size = 8,
                        label_x = 0.8)
  out <- plot_grid(pA, pB, nrow = 1, rel_widths = c(1, 0.85))
  out <- 
    ggdraw(out) + 
    draw_plot(leftCol, -0.01, 0.27, 0.13, 0.55) +
    draw_plot(rightCol, 0.875, 0.27, 0.13, 0.55) +
    # label nodes
    draw_text("a", x = 0.05, y = 0.93,  size = 7) +
    draw_text("b", x = 0.14, y = 0.85,  size = 7) +
    draw_text("c", x = 0.17, y = 0.51,  size = 7) +
    draw_text("d", x = 0.24, y = 0.237, size = 7) +
    draw_text("e", x = 0.33, y = 0.20,  size = 7) +
    draw_text("f", x = 0.96, y = 0.93,  size = 7) +
    draw_text("g", x = 0.87, y = 0.85,  size = 7) +
    draw_text("h", x = 0.83, y = 0.51,  size = 7) +
    draw_text("i", x = 0.76, y = 0.237, size = 7) +
    draw_text("j", x = 0.665, y = 0.20, size = 7)
  out <- plot_grid(out, get_legend(polDens1), nrow = 2, rel_heights = c(1, 0.04))
  ggsave(out, filename = paste0("figures/ouModel", ifelse(geo, "WithGeographicControl", ""), "/plotButterfly.pdf"), width = 6.5, height = 7)
  return(out)
}

# create traceplot for main parameters
plotTrace <- function(post, numTrees, numChains, geo = FALSE) {
  # number of iterations
  numIter <- dim(post$alpha_auto)[1] / numTrees
  # plot
  out <-
    tibble(
      chain      = rep(1:numChains, each = numIter),
      iter       = rep(1:numIter, times = numChains),
      `A[1,1]`   = post$A[1:(numChains*numIter),1,1],
      `A[1,2]`   = post$A[1:(numChains*numIter),1,2],
      `A[2,1]`   = post$A[1:(numChains*numIter),2,1],
      `A[2,2]`   = post$A[1:(numChains*numIter),2,2],
      `Q[1,1]`   = post$Q[1:(numChains*numIter),1,1],
      `Q[2,2]`   = post$Q[1:(numChains*numIter),2,2],
      `sigma[1]` = post$sigma[1:(numChains*numIter),1],
      `sigma[2]` = post$sigma[1:(numChains*numIter),2],
      `b[1]`     = post$b[1:(numChains*numIter),1],
      `b[2]`     = post$b[1:(numChains*numIter),2],
      `c1[1]`    = post$c1[1:(numChains*numIter),1],
      `c1[2]`    = post$c1[1:(numChains*numIter),2],
      `c1[3]`    = post$c1[1:(numChains*numIter),3],
      `c2[1]`    = post$c2[1:(numChains*numIter),1],
      `c2[2]`    = post$c2[1:(numChains*numIter),2],
      `c2[3]`    = post$c2[1:(numChains*numIter),3]
      ) %>%
    pivot_longer(cols = !c(chain, iter),
                 names_to = "parameter",
                 values_to = "value") %>%
    mutate(chain = factor(chain)) %>%
    ggplot(aes(x = iter, y = value, colour = chain)) +
    geom_line(alpha = 0.3) +
    facet_wrap(~ parameter, scales = "free") +
    labs(x = "Iteration", y = NULL) +
    theme_classic() +
    theme(strip.background = element_blank(),
          legend.position = "none")
  # save plot
  ggsave(out, filename = paste0("figures/ouModel", ifelse(geo, "WithGeographicControl", ""), "/plotTrace.pdf"), 
         width = 8, height = 5)
  return(out)
}

# simulate data from generative OU model
simOUData <- function(modelSimStan, phylo, iter, output="trait_values") {
  # tree for simulation
  tree <- phylo[[iter]]
  # cut up tree into segments
  times <- node.depth.edgelength(tree)
  # line up date of each node with the split points in the tree
  split_points <- sort(unique(times))
  node_time <- match(times, split_points)
  # create a sequence of nodes, respecting temporal order
  node_seq <- seq(from = 1, to = length(node_time))
  node_seq <- node_seq[order(node_time)]
  # find the "parent" node for each node and tip of the tree
  parent <- Ancestors(tree, node_seq, type = "parent")
  # parent time indicates amount of time since the parent node
  # scaled by the total depth of the tree
  parent_time <- rep(NA, length(node_seq))
  parent_time[1] <- -99 # placeholder for ancestral state
  for (i in 2:length(parent_time)) {
    parent_time[i] <- 
      (node.depth.edgelength(tree)[node_seq[i]] - node.depth.edgelength(tree)[parent[i]]) / 
      max(node.depth.edgelength(tree))
  }
  # get total num segments in the tree
  N_seg <- length(node_seq)
  # data list for stan
  data_list <- list(
    N_seg = N_seg,
    node_seq = node_seq,
    parent = parent,
    ts = parent_time,
    N = length(tree$tip.label),
    J = 2
  )
  # simulate 100 data frames from model
  fit_sim <- sampling(modelSimStan, data = data_list, 
                      algorithm = "Fixed_param", 
                      chains = 100, iter = 1, seed = iter)
  
  if (output == "trait_values") {
  out <- extract.samples(fit_sim)$y[,,]
  # randomly select one data frame with all four authority levels for both vars
  set.seed(iter)
  polAllFour <- apply(out[,,1], 1, function(x) length(table(x)))
  relAllFour <- apply(out[,,2], 1, function(x) length(table(x)))
  allFour <- polAllFour == 4 & relAllFour == 4
  out <- out[sample(which(allFour), 1),,]
  # as data frame
  out <- data.frame(language = tree$tip.label,
                    polAuth = out[,1],
                    relAuth = out[,2])
  }
  
  if (output == "parameters") {
    out <- extract.samples(fit_sim)
  }
  
  return(out)
}

# plot simulation results
plotSim <- function(simPost) {
  # fixed parameters
  f <- tibble(
    par = c("A[1,1]", "A[1,2]", "A[2,1]",
            "A[2,2]", "Q[1,1]", "Q[2,2]"),
    value = c(-1, 0, 4, -0.5, 6, 3)
  )
  # loop over models and extract posterior
  d <- tibble()
  for (i in 1:length(simPost)) {
    post <- tibble(
      model = i,
      `A[1,1]` = simPost[[i]]$A[,1,1],
      `A[1,2]` = simPost[[i]]$A[,1,2],
      `A[2,1]` = simPost[[i]]$A[,2,1],
      `A[2,2]` = simPost[[i]]$A[,2,2],
      `Q[1,1]` = simPost[[i]]$Q[,1,1],
      `Q[2,2]` = simPost[[i]]$Q[,2,2]
    )
    d <- rbind(d, post)
  }
  # plot
  out <-
    d %>%
    pivot_longer(cols = -model, names_to = "par") %>%
    group_by(model, par) %>%
    summarise(
      med = median(value),
      low = quantile(value, 0.025),
      upp = quantile(value, 0.975)
      ) %>%
    ggplot(aes(x = med, y = model, xmin = low, xmax = upp)) +
    geom_point(size = 0.005) +
    geom_errorbarh(size = 0.01, height = 0) +
    facet_wrap(. ~ par, nrow = 2, scales = "free_x") +
    geom_vline(data = f, aes(xintercept = value), 
               size = 0.5, colour = "red", linetype = "dashed") +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank())
  # save
  ggsave(out, filename = "figures/plotSim.pdf", width = 6, height = 4)
  return(out)
}

# generate prior
getPrior <- function(seed) {
  set.seed(seed)
  rnorm(200000, 0, 2)
}

# plot linguistic vs. geographic distance
plotDistance <- function(d, phylo, lonLat) {
  # get linguistic distances from maximum clade credibility tree
  cons <- mcc(phylo)
  lingDist <- cophenetic.phylo(cons)
  # wrangle lonLat
  lonLat <- 
    lonLat %>%
    # only languages in tree
    filter(Language %in% cons$tip.label) %>%
    # sort data to match phylogeny tips
    arrange(match(Language, cons$tip.label))
  # get geographic distances between longitude and latitude points
  geoDist <- distm(cbind(lonLat$Longitude, lonLat$Latitude))
  rownames(geoDist) <- colnames(geoDist) <- lonLat$Language
  # plot
  dist <-
    tibble(
      lingDist = lingDist[lower.tri(lingDist)],
      geoDist = log(geoDist[lower.tri(geoDist)])
    )
  cor <- round(cor(dist$lingDist, dist$geoDist), 2)
  out <-
    ggplot(dist, aes(x = geoDist, y = lingDist)) +
    geom_point(alpha = 0.2, size = 0.5) +
    geom_smooth(method = "lm") +
    labs(x = "Geographic distance (log)",
         y = "Linguistic distance",
         title = "Relationship between linguistic and geographic\ndistance across Austronesian societies",
         subtitle = paste0("Correlation = ", cor)) +
    theme_classic()
  # save plot
  ggsave(out, filename = paste0("figures/plotDist.pdf"), height = 5, width = 5)
  return(out)
}

# initialise phylogenetic GLMM in brms
getPhyloGLMMinitial <- function(d, phylo) {
  # use first tree for initialisation
  phylo <- phylo[[1]]
  # get covariance matrix
  A <- vcv.phylo(phylo, corr = TRUE)
  # add variable for modelling
  d$language2 <- d$language
  # initialise model
  bf1 <- bf(polAuth ~ 1 + (1 |a| gr(language, cov = A)) + (1 |b| language2)) + cumulative()
  bf2 <- bf(relAuth ~ 1 + (1 |a| gr(language, cov = A)) + (1 |b| language2)) + cumulative()
  brm(bf1 + bf2 + set_rescor(FALSE), data = d, data2 = list(A = A),
      prior = c(prior(normal(0, 2), class = Intercept, resp = polAuth),
                prior(normal(0, 2), class = Intercept, resp = relAuth),
                prior(exponential(2), class = sd, resp = polAuth),
                prior(exponential(2), class = sd, resp = relAuth),
                prior(lkj(3), class = cor)),
      chains = 0, seed = 2113, backend = "cmdstanr")
}

# fit phylogenetic GLMM in brms
fitPhyloGLMM <- function(modelPhyloGLMM, phylo, iter) {
  # loop over trees
  phylo <- phylo[[iter]]
  # get covariance matrix
  A <- vcv.phylo(phylo, corr = TRUE)
  # fit model
  update(modelPhyloGLMM, data2 = list(A = A),
         iter = 4000, chains = 1, cores = 1, seed = 2113,
         backend = "cmdstanr")
}

# plot results of phylogenetic GLMM
plotPhyloGLMM <- function(postPhyloGLMM) {
  # set theme for both plots
  plotTheme <-
    theme_classic() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank())
  # phylogenetic correlation
  pA <-
    ggplot() +
    geom_density(data = tibble(phyloCor = as.numeric(postPhyloGLMM[,,1])),
                 aes(x = phyloCor), fill = "indianred", colour = "indianred") +
    geom_pointinterval(data = tibble(median = median(postPhyloGLMM[,,1]),
                                     lower = hdi(postPhyloGLMM[,,1])[,1],
                                     upper = hdi(postPhyloGLMM[,,1])[,2],
                                     y = 0), 
                       aes(y = y, x = median, xmin = lower, xmax = upper),
                       fatten_point = 3, size = 6) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ylim(c(-0.5, 3.5)) +
    xlim(c(-1, 1)) +
    xlab("Phylogenetic correlation") +
    plotTheme
  # residual correlation
  pB <-
    ggplot() +
    geom_density(data = tibble(residCor = as.numeric(postPhyloGLMM[,,2])),
                 aes(x = residCor), fill = "grey", colour = "grey") +
    geom_pointinterval(data = tibble(median = median(postPhyloGLMM[,,2]),
                                     lower = hdi(postPhyloGLMM[,,2])[,1],
                                     upper = hdi(postPhyloGLMM[,,2])[,2],
                                     y = 0), 
                       aes(y = y, x = median, xmin = lower, xmax = upper),
                       fatten_point = 3, size = 6) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ylim(c(-0.5, 3.5)) +
    xlim(c(-1, 1)) +
    xlab("Residual correlation") +
    plotTheme
  # put together
  out <- plot_grid(pA, pB, nrow = 2)
  # save
  ggsave(out, filename = "figures/plotPhyloCor.pdf", width = 6, height = 4.5)
  return(out)
}

# fit OU stan model with geographic control
fitOUModelGeo <- function(modelStan, d, phylo, iter, lonLat) {
  # choose tree
  tree <- phylo[[iter]]
  # cut up tree into segments
  times <- node.depth.edgelength(tree)
  # line up date of each node with the split points in the tree
  split_points <- sort(unique(times))
  node_time <- match(times, split_points)
  # create a sequence of nodes, respecting temporal order
  node_seq <- seq(from = 1, to = length(node_time))
  node_seq <- node_seq[order(node_time)]
  # find the "parent" node for each node and tip of the tree
  parent <- Ancestors(tree, node_seq, type = "parent")
  # parent time indicates amount of time since the parent node
  # scaled by the total depth of the tree
  parent_time <- rep(NA, length(node_seq))
  parent_time[1] <- -99 # placeholder for ancestral state
  for (i in 2:length(parent_time)) {
    parent_time[i] <- 
      (node.depth.edgelength(tree)[node_seq[i]] - node.depth.edgelength(tree)[parent[i]]) / 
      max(node.depth.edgelength(tree))
  }
  # get total num segments in the tree
  N_seg <- length(node_seq)
  # construct geographic distance matrix
  d <- left_join(d, lonLat, by = c("language" = "Language"))
  geo_dist <- distm(cbind(d$Longitude, d$Latitude))
  geo_dist <- geo_dist / max(geo_dist)
  # data list for stan
  data_list <- list(
    N_seg = N_seg,
    node_seq = node_seq,
    parent = parent,
    ts = parent_time,
    N = length(tree$tip.label),
    J = 2,
    y = cbind(d$polAuth, d$relAuth),
    geo_dist = geo_dist
  )
  # fit model
  fit <- sampling(modelStan, data = data_list, chains = 1,
                  cores = 1, iter = 4000, init = "0",
                  control = list(adapt_delta = 0.999),
                  seed = 2113)
  return(fit)
}


# Simulate co-evolution of traits using Euler–Maruyama approximation for the system of stochastic differential equations
SDE_sim <- function(init_eta=c(0,0), parameters, time_depth, time=500, accuracy=10) {
  
  # Accuracy refers to the accuracy of our simulation. The system is co-evolving in continuous time, but we'll approximate it with small discrete time steps. I set the default as 10 year time steps
  n_times = time/accuracy
  
  # A matrix to hold our time series
  eta = matrix(NA, nrow=n_times, ncol=2)
  eta[1,] = init_eta
  
  A = parameters$A
  b = parameters$b
  G = parameters$G
  
  dt <- accuracy/time_depth
  
  for (t in 2:n_times) {
    eta[t,] = eta[t-1,] + (A %*% eta[t-1,] + b)*dt + (G %*% rnorm(2, 0, 1))*dt
  }
  
  return(eta)
}



