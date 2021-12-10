library(targets)
library(tidyverse)

set.seed(1)

N <- 200000 # num sims
tree_depth <- 5400 # tree depth
dt_500 <- 500/tree_depth # 500 year change

# prior change due to drift

sigma_prior <- abs(rnorm(N, 0, 1)) # current prior
prior_change_z <- rnorm(N, 0, 1) # amount of change due to drift, unscaled
prior_change_drift <- prior_change_z * sigma_prior * dt_500

# posterior change due to drift

tar_load(post)
# pol authority
sigma_post1 <- post$sigma[,1] # posterior 1
post_change_z1 <- rnorm(N, 0, 1) # should this still be std_normal?
post_change_drift1 <- post_change_z1 * sigma_post1 * dt_500
# rel authority
sigma_post2 <- post$sigma[,2] # posterior 2
post_change_z2 <- rnorm(N, 0, 1) # should this still be std_normal?
post_change_drift2 <- post_change_z2 * sigma_post2 * dt_500

# put together
tibble(
  drift = c(prior_change_drift, prior_change_drift,
            post_change_drift1, post_change_drift2),
  priorPost = rep(c("Prior", "Posterior"), each = N*2),
  variable = rep(rep(c("Political authority", "Religious authority"), each = N), times = 2)
) %>%
  mutate(priorPost = factor(priorPost, levels = c("Prior", "Posterior"))) %>%
  ggplot(aes(x = drift)) +
  geom_histogram() +
  facet_grid(priorPost ~ variable) +
  scale_x_continuous(name = "Expected drift in trait (500 years)", 
                     limits = c(-0.8, 0.8)) +
  theme_classic()
