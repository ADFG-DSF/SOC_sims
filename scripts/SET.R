# SET threshold

# Author: Adam Reimer
# Version: 2025-01-17

# Packages
packs <- c("tidyverse", "jagsUI", "mgcv")
lapply(packs, require, character.only = TRUE)

# source functions
function_files <- list.files(path=".\\functions")
lapply(function_files, function(x) source(paste0(".\\functions\\", x)))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Experimental script

# Assumed population ------------------------------------------------------
lnalpha <- 1.5
beta <- 0.0001
sigma <- 0.5
phi <- 0.0

# break up the simulations into chunks
years <- 40

# 70%/90% MSY goal bounds
ub_pctMSY = optimise(get_bounds, #'true' OYP bounds
                     get_Smsy(lnalpha, beta):(get_Smsy(lnalpha, beta)*5), 
                     lnalpha = lnalpha, 
                     beta = beta,
                     pct_MSY = 0.7,
                     correct = TRUE,
                     sigma = sigma,
                     phi = phi)$minimum
lb_pctMSY = optimise(get_bounds, #'true' OYP bounds
                     1:get_Smsy(lnalpha, beta), 
                     lnalpha = lnalpha, 
                     beta = beta,
                     pct_MSY = 0.9,
                     correct = TRUE,
                     sigma = sigma,
                     phi = phi)$minimum


# Gamma SET and eg range examples -----------------------------------------
# Dataset to demonstrate relationship between gamma set and gamma eg ranges
example_gamma <- 
  expand.grid(lnalpha = seq(0.5, 2.5, length.out = 5),
              beta = 1 / c(10000, 60000, 200000, 1000000),
              gamma = seq(1, 1.6, length.out = 4)) %>%
  mutate(Smax = 1/beta,
         Rmax = Ricker(lnalpha, beta, Smax),
         a = gamma_par(Smax, Rmax, gamma)[[1]],
         b = gamma_par(Smax, Rmax, gamma)[[2]])%>%
  rowwise() %>%
  mutate(Smsy =   
           optimize(
             function(a, b, g, x){
               SRgamma(alpha = a, beta = b, gamma = g, S = x) - x
             },
             interval = c(0, 4 * 1 / beta),
             maximum = TRUE,
             a = a,
             b = b,
             g = gamma)$maximum,
         Rmsy = SRgamma(alpha = a, beta = b, gamma = gamma, S = Smsy),
         MSY = Rmsy - Smsy) %>%
  mutate(lb = 
           optimize(
             function(a, b, g, pct_MSY, MSY, x){
               ((SRgamma(alpha = a, beta = b, gamma = g, S = x) - x) - pct_MSY * MSY)^2
             },
             interval = c(Smsy / 10, Smsy),
             a = a,
             b = b,
             g = gamma,
             pct_MSY = 0.9,
             MSY = MSY)$minimum,
         ub = 
           optimize(
             function(a, b, g, pct_MSY, MSY, x){
               ((SRgamma(alpha = a, beta = b, gamma = g, S = x) - x) - pct_MSY * MSY)^2
             },
             interval = c(Smsy, Smsy * 3),
             a = a,
             b = b,
             g = gamma,
             pct_MSY = 0.7,
             MSY = MSY)$minimum,
         set = ifelse(gamma != 1, (gamma - 1) / b, 10)) %>% # set a genetics based estimate as the lower
  ungroup()

#Plot one beta at a time
plot_gamma <- 
  example_gamma %>%
  filter(Smax == 1000000)

# R vrs. S with goal ranges and set
# Note SET is invariant to alpha, (gamma - 1) / b
plot_gamma %>%
  slice(rep(1:dim(plot_gamma)[1], each = 101)) %>%
  mutate(S0 = 
           rep(
             seq(0, round(max(plot_gamma$ub) * 1.2, -3), by = round(max(plot_gamma$ub) * 1.2, -3) / 100), 
             times = dim(plot_gamma)[1]),
         S = ifelse(S0 == 0, 10, S0),
         R_fit = SRgamma(alpha = a, beta = b, gamma = gamma, S = S),
         Y_fit = R_fit - S,
         lnRS_fit = log(R_fit / S)) %>%
  ggplot(aes(x = S, y = R_fit)) +
  geom_line() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_vline(aes(xintercept = set, color = as.character(gamma))) +
  geom_rect(aes(xmin = lb, xmax = ub, ymin = 0, ymax = Inf, fill = as.character(gamma)), alpha= 0.01) +
  scale_x_continuous(limits = c(0, max(plot_gamma$ub) * 1.2)) +
  facet_grid(lnalpha ~ gamma, scales = "free")

# Y vrs. S with goal ranges and set
plot_gamma %>%
  slice(rep(1:dim(plot_gamma)[1], each = 101)) %>%
  mutate(S0 = 
           rep(
             seq(0, round(max(plot_gamma$ub) * 1.2, -3), by = round(max(plot_gamma$ub) * 1.2, -3) / 100),  
             times = dim(plot_gamma)[1]),
         S = ifelse(S0 == 0, 10, S0),
         R_fit = SRgamma(alpha = a, beta = b, gamma = gamma, S = S),
         Y_fit = R_fit - S,
         lnRS_fit = log(R_fit / S)) %>%
  filter(Y_fit > 0) %>%
  ggplot(aes(x = S, y = Y_fit)) +
  geom_line() +
  geom_vline(aes(xintercept = set, color = as.character(gamma))) +
  geom_rect(aes(xmin = lb, xmax = ub, ymin = 0, ymax = Inf, fill = as.character(gamma)), alpha= 0.01) +
  scale_x_continuous(limits = c(0, max(plot_gamma$ub * 1.2))) +
  facet_grid(lnalpha ~ gamma, scales = "free")

# ln(R/S) vrs. S with goal ranges and set
plot_gamma %>%
  slice(rep(1:dim(plot_gamma)[1], each = 101)) %>%
  mutate(S0 = 
           rep(
             seq(0, round(max(plot_gamma$ub) * 1.2, -3), by = round(max(plot_gamma$ub) * 1.2, -3) / 100), 
             times = dim(plot_gamma)[1]),
         S = ifelse(S0 == 0, 10, S0),
         R_fit = SRgamma(alpha = a, beta = b, gamma = gamma, S = S),
         Y_fit = R_fit - S,
         lnRS_fit = log(R_fit / S)) %>%
  ggplot(aes(x = S, y = lnRS_fit)) +
  geom_line() +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = set, color = as.character(gamma))) +
  geom_rect(aes(xmin = lb, xmax = ub, ymin = -Inf, ymax = Inf, fill = as.character(gamma)), alpha= 0.01) +
  scale_x_continuous(limits = c(0, max(plot_gamma$ub * 1.6))) +
  coord_cartesian(ylim = c(-1, 3)) +
  facet_grid(lnalpha ~ gamma, scales = "free")

# set / lb 
# set a small fraction of eg lb except when depensation is strong and productivity is weak.
example_gamma %>% 
  mutate(pct = set / lb) %>%
  ggplot(aes(x = gamma, y = pct, color = as.character(lnalpha))) +
  geom_line() +
  geom_jitter(width = 0.015, height = 0.01)





# Ricker SR, Gamma estimates -------------------------------------------
# Can we estimate the gamma parameter when the data comes from a Ricker?
# Simulate a dataset
dat_Ricker <-
  sim_Ricker(lnalpha, 
             beta, 
             sigW = sigma, 
             phi = phi, 
             age0 = c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02), 
             Sims0 = 1520, 
             Hfun = H_goal,
             lb_goal = lb_pctMSY,
             ub_goal = ub_pctMSY,
             sigF = 0.4,
             sigN = 0.0) %>% 
  select(sim, R, S) %>%
  mutate(group = (sim - 1) %/% years + 1,
         lnRS = ifelse(is.finite(log(R/S)) & !is.na(log(R/S)), log(R/S), NA),
         Y = R - S) %>%
  group_by(group) %>%
  mutate(year = row_number())

# Look the the simulated data
# The SR relationship that generated it
# The similar Gamma
# a smooth through the data
# just trying to verify we have some reasonable expectation to estimate depensation
dat_Ricker %>%
  filter(group < 7) %>%
  ggplot(aes(x = S, y = lnRS, color = as.character(group))) +
  geom_point() +
  geom_smooth(se = FALSE) + 
  geom_function(fun = function(x) log((0.0739 * x^1.5 * exp(- 0.00015 * x)) / x)) +
  geom_function(fun = function(x) log((exp(1.5) * x * exp(- 0.0001 * x)) / x), linetype = 2) +
  scale_x_continuous(limits = c(0, round(max(dat_Ricker$S), -3)), n.breaks = 3) +
  annotate(geom = "rect", xmin = lb_pctMSY, xmax = ub_pctMSY, ymin = 0, ymax = Inf, alpha = 0.2, color = NA) +
  geom_hline(aes(yintercept = 0)) +
  facet_wrap(. ~ group)

# * JAGS ----
# data
jagsdat_Ricker <- 
  list(nyrs = max(dat_Ricker$year),
       G = max(dat_Ricker$group),
       S = matrix(dat_Ricker$S, nrow = max(dat_Ricker$year), ncol = max(dat_Ricker$group), byrow = FALSE),
       R = matrix(dat_Ricker$R, nrow = max(dat_Ricker$year), ncol = max(dat_Ricker$group), byrow = FALSE),
       g = dat_Ricker$group,
       ar1 = 0,
       kf = 0)

parameters <- c('beta', 'sigma', 'lnalpha', 'gamma')

# Run
post_Ricker <- 
  jags(data = jagsdat_Ricker,
       parameters.to.save = parameters,
       model.file = ".\\scripts\\gamma_groups.txt",
       n.chains = 3,
       n.iter = 1e4,
       n.burnin = 5e3,
       n.thin = 10,
       parallel = TRUE,
       store.data = TRUE
)

est_Ricker <- 
  data.frame(group = 1:post_Ricker$data$G,
             a = exp(post_Ricker$q50$lnalpha),
             beta = post_Ricker$q50$beta,
             gamma = post_Ricker$q50$gamma)

# * Plot estimates and data -------------------------------------------------
# The assumed relationship for each group
assumed_Ricker <- 
  data.frame(group = 1:post_Ricker$data$G,
             lnalpha = lnalpha,
             beta = beta) %>%
  slice(rep(1:post_Ricker$data$G, each = round(max(dat_Ricker$S), -3) %/% 250 + 1)) %>%
  mutate(S0 = rep(seq(0, round(max(dat_Ricker$S), -3), by = 250), times = post_Ricker$data$G),
         S = ifelse(S0 == 0, 10, S0),
         R_fit = Ricker(lnalpha = lnalpha, beta = beta, S = S),
         Y_fit = R_fit - S,
         lnRS_fit = log(R_fit / S),
         source = "assumed") %>%
  pivot_longer(dplyr::ends_with("fit"), 
               names_to = "stat", 
               values_to = "value", 
               names_pattern = "(.*)_fit") %>%
  select(S, group, source, stat, value)

# Data for each group
data_Ricker <- 
  dat_Ricker %>% 
  pivot_longer(c(R, Y, lnRS), names_to = "stat", values_to = "value") %>%
  mutate(source = "assumed")

# The fitted relationship for each group
estimated_Ricker <- 
  est_Ricker %>%
  slice(rep(1:post_Ricker$data$G, each = round(max(dat_Ricker$S), -3) %/% 250 + 1)) %>%
  mutate(S0 = rep(seq(0, round(max(dat_Ricker$S), -3), by = 250), times = post_Ricker$data$G),
         S = ifelse(S0 == 0, 10, S0),
         R_fit = SRgamma(alpha = a, beta = beta, gamma = gamma, S = S),
         Y_fit = R_fit - S,
         lnRS_fit = log(R_fit / S),
         source = "estimated") %>%
  pivot_longer(dplyr::ends_with("fit"), 
               names_to = "stat", 
               values_to = "value", 
               names_pattern = "(.*)_fit")%>%
  select(S, group, source, stat, value)

ref_Ricker <- 
  est_Ricker %>%
  rowwise() %>%
  mutate(Smsy =   
           optimize(
             function(a, b, g, x){
               SRgamma(alpha = a, beta = b, gamma = g, S = x) - x
             },
             interval = c(0, 4 * 1 / beta),
             maximum = TRUE,
             a = a,
             b = beta,
             g = gamma)$maximum,
         Rmsy = SRgamma(alpha = a, beta = beta, gamma = gamma, S = Smsy),
         MSY = Rmsy - Smsy) %>%
  mutate(lb = 
           optimize(
             function(a, b, g, pct_MSY, MSY, x){
               ((SRgamma(alpha = a, beta = b, gamma = g, S = x) - x) - pct_MSY * MSY)^2
             },
             interval = c(50, Smsy),
             a = a,
             b = beta,
             g = gamma,
             pct_MSY = 0.9,
             MSY = MSY)$minimum,
         ub = 
           optimize(
             function(a, b, g, pct_MSY, MSY, x){
               ((SRgamma(alpha = a, beta = b, gamma = g, S = x) - x) - pct_MSY * MSY)^2
             },
             interval = c(Smsy, Smsy * 3),
             a = a,
             b = beta,
             g = gamma,
             pct_MSY = 0.7,
             MSY = MSY)$minimum,
         set = 
           optimize(
             function(a, b, g, x){
               -log(SRgamma(alpha = a, beta = b, gamma = g, S = x) / x)
             },
             interval = c(10, Smsy),
             #maximum = TRUE,
             a = a,
             b = beta,
             g = gamma)$minimum)

# ** plot -------------
# Can we estimate the gamma parameter when the data comes from a Ricker?
# mixed bag. Estimate both depensation, compensation and hyper compensation depending on the dataset.
pull_groups <- sample(1:max(data_Ricker$group), 3, replace = FALSE)
rbind(estimated_Ricker, assumed_Ricker) %>%
  filter(group %in% pull_groups) %>%
  ggplot(aes(x = S, y = value, color = source)) +
  geom_line() +
  geom_point(data = data_Ricker[data_Ricker$group %in% pull_groups, ]) +
  scale_x_continuous(limits = c(0, 20000), n.breaks = 5) +
  geom_rect(data = ref_Ricker[ref_Ricker$group %in% pull_groups, ], 
            aes(xmin = lb, xmax = ub, ymin = -Inf, ymax = Inf), 
            inherit.aes = FALSE, alpha= 0.2) +
  annotate(geom = "rect", 
           xmin = lb_pctMSY, xmax = ub_pctMSY, ymin = -Inf, ymax = Inf, 
           alpha = 0.2, fill = "orange") +
  scale_color_manual(values = c("orange", "grey")) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(data = ref_Ricker[ref_Ricker$group %in% pull_groups, ], 
             aes(xintercept = set),
             linetype = 2) +
  geom_vline(data = ref_Ricker[ref_Ricker$group %in% pull_groups, ], 
             aes(xintercept = Smsy), 
             linetype = 3) +
  facet_grid(stat ~ group, scales = "free_y")






# comparable Gamma SR -------------------------------------------
# Can we estimate the gamma parameter when the data comes from a Gamma?
# arbitary gamma but contains obvious depensatory dynamics
gamma <- 1.5

# gamma params with the same Smsy and Rmsy
gamma_params <- gamma_par(1/beta, Ricker(lnalpha, beta, 1/beta), gamma)

# Simulate a dataset
dat_gamma <-
  sim_SRgamma(gamma_params[[1]], 
              beta = gamma_params[[2]], 
              gamma = gamma_params[[3]], 
              sigW = sigma, 
              phi = phi, 
              age0 = c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02), 
              Sims0 = 1520, 
              Hfun = H_goal,
              lb_goal = lb_pctMSY,
              ub_goal = ub_pctMSY,
              sigF = 0.4,
              sigN = 0.0) %>% 
  select(sim, R, S) %>%
  mutate(group = (sim - 1) %/% years + 1,
         lnRS = ifelse(is.finite(log(R/S)) & !is.na(log(R/S)), log(R/S), NA),
         Y = R - S) %>%
  group_by(group) %>%
  mutate(year = row_number())

# Look the the simulated data
# The SR relationship that generated it
# The similar Ricker
# a smooth through the data
# just trying to verify we have some reasonable expectation to estimate depensation
dat_gamma %>%
  filter(group < 7) %>%
  ggplot(aes(x = S, y = lnRS, color = as.character(group))) +
  geom_point() +
  geom_smooth(se = FALSE) + 
  geom_function(fun = function(x) log((0.0739 * x^1.5 * exp(- 0.00015 * x)) / x)) +
  geom_function(fun = function(x) log((exp(1.5) * x * exp(- 0.0001 * x)) / x), linetype = 2) +
  scale_x_continuous(limits = c(0, round(max(dat_gamma$S), -3)), n.breaks = 3) +
  annotate(geom = "rect", xmin = lb_pctMSY, xmax = ub_pctMSY, ymin = 0, ymax = Inf, alpha = 0.2, color = NA) +
  geom_hline(aes(yintercept = 0)) +
  facet_wrap(. ~ group)

# * JAGS ----
# data
jagsdat_gamma <- 
  list(nyrs = max(dat_gamma$year),
       G = max(dat_gamma$group),
       S = matrix(dat_gamma$S, nrow = max(dat_gamma$year), ncol = max(dat_gamma$group), byrow = FALSE),
       R = matrix(dat_gamma$R, nrow = max(dat_gamma$year), ncol = max(dat_gamma$group), byrow = FALSE),
       g = dat_gamma$group,
       ar1 = 0,
       kf = 0)

# Run
post_gamma <- jags(data = jagsdat_gamma,
                   parameters.to.save = parameters,
                   model.file = ".\\scripts\\gamma_groups.txt",
                   n.chains = 3,
                   n.iter = 1e4,
                   n.burnin = 5e3,
                   n.thin = 10,
                   parallel = TRUE,
                   store.data = TRUE
)

est_gamma <- 
  data.frame(group = 1:post_gamma$data$G,
             a = exp(post_gamma$q50$lnalpha),
             beta = post_gamma$q50$beta,
             gamma = post_gamma$q50$gamma)

# * Plot estimates and data -------------------------------------------------
# The assumed relationship for each group
assumed_gamma <- 
  data.frame(group = 1:post_gamma$data$G,
             a = gamma_params[[1]],
             beta = gamma_params[[2]],
             gamma = gamma_params[[3]]) %>%
  slice(rep(1:post_gamma$data$G, each = round(max(dat_gamma$S), -3) %/% 250 + 1)) %>%
  mutate(S0 = rep(seq(0, round(max(dat_gamma$S), -3), by = 250), times = post_gamma$data$G),
         S = ifelse(S0 == 0, 10, S0),
         R_fit = SRgamma(alpha = a, beta = beta, gamma = gamma, S = S),
         Y_fit = R_fit - S,
         lnRS_fit = log(R_fit / S),
         source = "assumed") %>%
  pivot_longer(dplyr::ends_with("fit"), 
               names_to = "stat", 
               values_to = "value", 
               names_pattern = "(.*)_fit")

# Data for each group
data_gamma <- 
  dat_gamma %>% 
  pivot_longer(c(R, Y, lnRS), names_to = "stat", values_to = "value") %>%
  mutate(source = "assumed")

# The fitted relationship for each group
estimated_gamma <- 
  est_gamma %>%
  slice(rep(1:post_gamma$data$G, each = round(max(dat_gamma$S), -3) %/% 250 + 1)) %>%
  mutate(S0 = rep(seq(0, round(max(dat_gamma$S), -3), by = 250), times = post_gamma$data$G),
         S = ifelse(S0 == 0, 10, S0),
         R_fit = SRgamma(alpha = a, beta = beta, gamma = gamma, S = S),
         Y_fit = R_fit - S,
         lnRS_fit = log(R_fit / S),
         source = "estimated") %>%
  pivot_longer(dplyr::ends_with("fit"), 
               names_to = "stat", 
               values_to = "value", 
               names_pattern = "(.*)_fit")

ref_gamma <- 
  est_gamma %>%
  rowwise() %>%
  mutate(Smsy =   
           optimize(
             function(a, b, g, x){
               SRgamma(alpha = a, beta = b, gamma = g, S = x) - x
             },
             interval = c(0, 4 * 1 / beta),
             maximum = TRUE,
             a = a,
             b = beta,
             g = gamma)$maximum,
         Rmsy = SRgamma(alpha = a, beta = beta, gamma = gamma, S = Smsy),
         MSY = Rmsy - Smsy) %>%
  mutate(lb = 
           optimize(
             function(a, b, g, pct_MSY, MSY, x){
               ((SRgamma(alpha = a, beta = b, gamma = g, S = x) - x) - pct_MSY * MSY)^2
             },
             interval = c(50, Smsy),
             a = a,
             b = beta,
             g = gamma,
             pct_MSY = 0.9,
             MSY = MSY)$minimum,
         ub = 
           optimize(
             function(a, b, g, pct_MSY, MSY, x){
               ((SRgamma(alpha = a, beta = b, gamma = g, S = x) - x) - pct_MSY * MSY)^2
             },
             interval = c(Smsy, Smsy * 3),
             a = a,
             b = beta,
             g = gamma,
             pct_MSY = 0.7,
             MSY = MSY)$minimum,
         set = 
           optimize(
             function(a, b, g, x){
               -log(SRgamma(alpha = a, beta = b, gamma = g, S = x) / x)
             },
             interval = c(10, Smsy),
             #maximum = TRUE,
             a = a,
             b = beta,
             g = gamma)$minimum)

# ** plot -------------
# Can we estimate the gamma parameter when the data comes from a Gamma?
# Mostly yes.
pull_groups <- sample(1:max(dat_gamma$group), 3, replace = FALSE)
rbind(estimated_gamma, assumed_gamma) %>%
  filter(group %in% pull_groups) %>%
  ggplot(aes(x = S, y = value, color = source)) +
  geom_line() +
  geom_point(data = data_gamma[data_gamma$group %in% pull_groups, ]) +
  scale_x_continuous(limits = c(0, 20000), n.breaks = 5) +
  geom_rect(data = ref_gamma[ref_gamma$group %in% pull_groups, ], 
            aes(xmin = lb, xmax = ub, ymin = -Inf, ymax = Inf), 
            inherit.aes = FALSE, alpha= 0.2) +
  annotate(geom = "rect", 
           xmin = lb_pctMSY, xmax = ub_pctMSY, ymin = -Inf, ymax = Inf, 
           alpha = 0.2, fill = "orange") +
  scale_color_manual(values = c("orange", "grey")) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(data = ref_gamma[ref_gamma$group %in% pull_groups, ], 
             aes(xintercept = set),
             linetype = 2) +
  geom_vline(data = ref_gamma[ref_gamma$group %in% pull_groups, ], 
             aes(xintercept = Smsy), 
             linetype = 3) +
  facet_grid(stat ~ group, scales = "free_y")

# Compare parameter estimates
data.frame(gamma_Ricker = post_Ricker$q50$gamma,
           gamma_gamma = post_gamma$q50$gamma) %>%
  pivot_longer(cols = starts_with("gamma"), 
               names_to = "SR", 
               names_pattern = "gamma_(.*)",
               values_to = "gamma") %>%
  ggplot(aes(x = SR, y = gamma)) +
    geom_boxplot()






# time varying Ricker SR, Gamma estimates -------------------------------------------
# Can we estimate the gamma parameter when the data comes from a time-varying Ricker?
# First 20 yrs from "normal" productivity, next 20 year Seq = the lower bound
# Simulate a dataset
dat_tvRicker <-
  sim_Ricker(rep(c(rep(lnalpha, 20), rep(lnalpha * 0.25, 20)), times = 38), 
             beta, 
             sigW = sigma, 
             phi = phi, 
             age0 = c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02), 
             Sims0 = 1520, 
             Hfun = H_goal,
             lb_goal = lb_pctMSY,
             ub_goal = ub_pctMSY,
             power = 0.5,
             sigF = 0.2,
             sigN = 0.2) %>% 
  select(sim, R, S) %>%
  mutate(group = (sim - 1) %/% years + 1,
         lnRS = ifelse(is.finite(log(R/S)) & !is.na(log(R/S)), log(R/S), NA),
         Y = R - S) %>%
  group_by(group) %>%
  mutate(year = row_number())

# Look the the simulated data
# The SR relationship that generated it
# The similar Gamma
# a smooth through the data
# just trying to verify we have some reasonable expectation to estimate depensation

dat_tvRicker %>%
  ungroup() %>%
  filter(group %in% sample(1:38, 6)) %>%
  ggplot(aes(x = S, y = lnRS, color = as.character(group))) +
  geom_point() +
  geom_smooth(se = FALSE) + 
  geom_function(fun = function(x) log((0.0739 * x^1.5 * exp(- 0.00015 * x)) / x)) +
  geom_function(fun = function(x) log((exp(lnalpha) * x * exp(- 0.0001 * x)) / x), linetype = 2) +
  geom_function(fun = function(x) log((exp(lnalpha * 0.25) * x * exp(- 0.0001 * x)) / x), linetype = 3) +
  scale_x_continuous(limits = c(0, round(max(dat_tvRicker$S), -3)), n.breaks = 3) +
  annotate(geom = "rect", xmin = lb_pctMSY, xmax = ub_pctMSY, ymin = 0, ymax = Inf, alpha = 0.2, color = NA) +
  geom_hline(aes(yintercept = 0)) +
  facet_wrap(. ~ group)

# * JAGS ----
# data
jagsdat_tvRicker <- 
  list(nyrs = max(dat_tvRicker$year),
       G = max(dat_tvRicker$group),
       S = matrix(dat_tvRicker$S, nrow = max(dat_tvRicker$year), ncol = max(dat_tvRicker$group), byrow = FALSE),
       R = matrix(dat_tvRicker$R, nrow = max(dat_tvRicker$year), ncol = max(dat_tvRicker$group), byrow = FALSE),
       g = dat_tvRicker$group,
       ar1 = 0,
       kf = 0)

parameters <- c('beta', 'sigma', 'lnalpha', 'gamma')

# Run
post_tvRicker <- 
  jags(data = jagsdat_tvRicker,
       parameters.to.save = parameters,
       model.file = ".\\scripts\\gamma_groups.txt",
       n.chains = 3,
       n.iter = 1e4,
       n.burnin = 5e3,
       n.thin = 10,
       parallel = TRUE,
       store.data = TRUE
  )

est_tvRicker <- 
  data.frame(group = 1:post_tvRicker$data$G,
             a = exp(post_tvRicker$q50$lnalpha),
             beta = post_tvRicker$q50$beta,
             gamma = post_tvRicker$q50$gamma)

# * Plot estimates and data -------------------------------------------------
# The assumed relationship for each group
assumed_tvRicker <- 
  data.frame(group = 1:post_tvRicker$data$G,
             lnalpha = lnalpha,
             beta = beta) %>%
  slice(rep(1:post_tvRicker$data$G, each = round(max(dat_tvRicker$S), -3) %/% 250 + 1)) %>%
  mutate(S0 = rep(seq(0, round(max(dat_tvRicker$S), -3), by = 250), times = post_tvRicker$data$G),
         S = ifelse(S0 == 0, 10, S0),
         R1_fit = Ricker(lnalpha = lnalpha, beta = beta, S = S),
         Y1_fit = R1_fit - S,
         lnRS1_fit = log(R1_fit / S),
         R2_fit = Ricker(lnalpha = lnalpha * 0.25, beta = beta, S = S),
         Y2_fit = R2_fit - S,
         lnRS2_fit = log(R2_fit / S)) %>%
  pivot_longer(dplyr::ends_with("fit"), 
               names_to = "stat0", 
               values_to = "value", 
               names_pattern = "(.*)_fit") %>%
  mutate(stat = gsub("(.*)\\d", "\\1", stat0),
         regime = gsub("(.*)(\\d)", "\\2", stat0),
         source = paste0("assumed", "_", regime)) %>%
  select(S, group, source, stat, value)

# Data for each group
data_tvRicker <- 
  dat_tvRicker %>% 
  pivot_longer(c(R, Y, lnRS), names_to = "stat", values_to = "value") %>%
  mutate(source = ifelse(year <= 20, "assumed_1", "assumed_2")) #not linked to simulation dataset!!!

# The fitted relationship for each group
estimated_tvRicker <- 
  est_tvRicker %>%
  slice(rep(1:post_tvRicker$data$G, each = round(max(dat_tvRicker$S), -3) %/% 250 + 1)) %>%
  mutate(S0 = rep(seq(0, round(max(dat_tvRicker$S), -3), by = 250), times = post_tvRicker$data$G),
         S = ifelse(S0 == 0, 10, S0),
         R_fit = SRgamma(alpha = a, beta = beta, gamma = gamma, S = S),
         Y_fit = R_fit - S,
         lnRS_fit = log(R_fit / S),
         source = "estimated") %>%
  pivot_longer(dplyr::ends_with("fit"), 
               names_to = "stat", 
               values_to = "value", 
               names_pattern = "(.*)_fit")%>%
  select(S, group, source, stat, value)

ref_tvRicker <- 
  est_tvRicker %>%
  rowwise() %>%
  mutate(Smsy =   
           optimize(
             function(a, b, g, x){
               SRgamma(alpha = a, beta = b, gamma = g, S = x) - x
             },
             interval = c(0, 4 * 1 / beta),
             maximum = TRUE,
             a = a,
             b = beta,
             g = gamma)$maximum,
         Rmsy = SRgamma(alpha = a, beta = beta, gamma = gamma, S = Smsy),
         MSY = Rmsy - Smsy) %>%
  mutate(lb = 
           optimize(
             function(a, b, g, pct_MSY, MSY, x){
               ((SRgamma(alpha = a, beta = b, gamma = g, S = x) - x) - pct_MSY * MSY)^2
             },
             interval = c(50, Smsy),
             a = a,
             b = beta,
             g = gamma,
             pct_MSY = 0.9,
             MSY = MSY)$minimum,
         ub = 
           optimize(
             function(a, b, g, pct_MSY, MSY, x){
               ((SRgamma(alpha = a, beta = b, gamma = g, S = x) - x) - pct_MSY * MSY)^2
             },
             interval = c(Smsy, Smsy * 3),
             a = a,
             b = beta,
             g = gamma,
             pct_MSY = 0.7,
             MSY = MSY)$minimum,
         set = 
           optimize(
             function(a, b, g, x){
               -log(SRgamma(alpha = a, beta = b, gamma = g, S = x) / x)
             },
             interval = c(10, Smsy),
             #maximum = TRUE,
             a = a,
             b = beta,
             g = gamma)$minimum)

# ** plot -------------
# Can we estimate the gamma parameter when the data comes from a time-varying Ricker?
# poorly. Often beta is underestimated gamma all over the place. I think this is expected as
# most of our data from high escapements comes from the "normal" Ricker while 
# most of our small escapements come from the"low" Ricker. 
# To me this suggested we need to account time-varying productivity during estimation 
pull_groups <- sample(1:max(data_tvRicker$group), 3, replace = FALSE)
rbind(estimated_tvRicker, assumed_tvRicker) %>%
  filter(group %in% pull_groups) %>%
  ggplot(aes(x = S, y = value, color = source)) +
  geom_line() +
  geom_point(data = data_tvRicker[data_tvRicker$group %in% pull_groups, ]) +
  scale_x_continuous(limits = c(0, 20000), n.breaks = 5) +
  geom_rect(data = ref_tvRicker[ref_tvRicker$group %in% pull_groups, ], 
            aes(xmin = lb, xmax = ub, ymin = -Inf, ymax = Inf), 
            inherit.aes = FALSE, alpha= 0.2) +
  annotate(geom = "rect", 
           xmin = lb_pctMSY, xmax = ub_pctMSY, ymin = -Inf, ymax = Inf, 
           alpha = 0.2, fill = "orange") +
  scale_color_manual(values = c("orange", "red","grey")) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(data = ref_tvRicker[ref_tvRicker$group %in% pull_groups, ], 
             aes(xintercept = set),
             linetype = 2) +
  geom_vline(data = ref_tvRicker[ref_tvRicker$group %in% pull_groups, ], 
             aes(xintercept = Smsy), 
             linetype = 3) +
  facet_grid(stat ~ group, scales = "free_y")

# Compare parameter estimates
data.frame(gamma_Ricker = post_Ricker$q50$gamma,
           gamma_tvRicker = post_tvRicker$q50$gamma,
           gamma_gamma = post_gamma$q50$gamma) %>%
  pivot_longer(cols = starts_with("gamma"), 
               names_to = "SR", 
               names_pattern = "gamma_(.*)",
               values_to = "gamma") %>%
  ggplot(aes(x = SR, y = gamma)) +
  geom_boxplot()

data.frame(beta_Ricker = post_Ricker$q50$beta,
           beta_tvRicker = post_tvRicker$q50$beta,
           beta_gamma = post_gamma$q50$beta) %>%
  pivot_longer(cols = starts_with("beta"), 
               names_to = "SR", 
               names_pattern = "beta_(.*)",
               values_to = "beta") %>%
  ggplot(aes(x = SR, y = beta)) +
  geom_boxplot()

data.frame(lnalpha_Ricker = post_Ricker$q50$lnalpha,
           lnalpha_tvRicker = post_tvRicker$q50$lnalpha,
           lnalpha_gamma = post_gamma$q50$lnalpha) %>%
  pivot_longer(cols = starts_with("lnalpha"), 
               names_to = "SR", 
               names_pattern = "lnalpha_(.*)",
               values_to = "lnalpha") %>%
  ggplot(aes(x = SR, y = lnalpha)) +
  geom_boxplot()



# tvRicker-rt --------------------------------------------------------
# Can we estimate the gamma parameter when the data comes from a time varying Ricker 
# with a regime transition model?
# * JAGS ----
# Trying to take beta mostly from the "normal" regime.
# gamma coming only from the low regime.
# convergence is not good enough
# data
post_tvRicker_rt <- 
  jags(data = jagsdat_tvRicker,
       parameters.to.save = c("beta", "lnalpha", "gamma", "pi", "lambda", "sigma"),
       model.file = ".\\scripts\\rt_groups.txt",
       n.chains = 3,
       n.iter = 1e5,
       n.burnin = 5e4,
       n.thin = 10,
       parallel = TRUE,
       store.data = TRUE
  )

est_tvRicker_rt <- 
  data.frame(group = 1:post_tvRicker_rt$data$G,
             a1 = exp(post_tvRicker_rt$q50$lnalpha[1, ]),
             a2 = exp(post_tvRicker_rt$q50$lnalpha[2, ]),
             b1 = post_tvRicker_rt$q50$beta[1, ],
             b2 = post_tvRicker_rt$q50$beta[2, ],
             g1 = post_tvRicker_rt$q50$gamma[1, ],
             g2 = post_tvRicker_rt$q50$gamma[2, ]) %>%
  pivot_longer(cols = c("a1", "a2", "b1", "b2", "g1", "g2"), 
               names_to = c("param", "regime"),
               names_pattern = "(.)(\\d)",
               values_to = "value") %>%
  pivot_wider(id_cols = c(group, regime), 
               names_from = "param", 
               values_from = "value")

# * Plot estimates and data -------------------------------------------------
assumed_tvRicker_rt <-
  assumed_tvRicker %>%
  mutate(regime = gsub("assumed_(\\d)", "\\1", source),
         source = gsub("(assumed)_(\\d)", "\\1", source))

data_tvRicker_rt <-
  data_tvRicker %>%
  mutate(regime = gsub("assumed_(\\d)", "\\1", source),
         source = gsub("(assumed)_(\\d)", "\\1", source))

# The fitted relationship for each group
estimated_tvRicker_rt <- 
  est_tvRicker_rt %>%
  slice(rep(1:(2 * post_tvRicker$data$G), each = round(max(dat_tvRicker$S), -3) %/% 250 + 1)) %>%
  mutate(S0 = rep(seq(0, round(max(dat_tvRicker$S), -3), by = 250), times = 2 * post_tvRicker$data$G),
         S = ifelse(S0 == 0, 10, S0),
         R_fit = SRgamma(alpha = a, beta = b, gamma = g, S = S),
         Y_fit = R_fit - S,
         lnRS_fit = log(R_fit / S),
         source = "estimated") %>%
  pivot_longer(dplyr::ends_with("fit"), 
               names_to = "stat", 
               values_to = "value", 
               names_pattern = "(.*)_fit")%>%
  select(S, group, regime, source, stat, value)

# * Simulate yield across regimes ------------------------------------------------
# this allows for Smsy and eg range estimates to account for both regimes.
# Not 100% sure we want to do this... Perhaps we should not be targeting yield 
# the the stock is barely replacing itself.

#indices to use to draw 1000 samples from posterior
draw <- sample(1:dim(post_tvRicker_rt$sims.list$lnalpha)[1],
               size = 1000,
               replace = FALSE)

#draw pi
sims_pi <- list(pi = post_tvRicker_rt$sims.list$pi[draw, , ,])
stationary_pi <- list()
for(g in 1:dim(sims_pi$pi)[4]){
  stationary_pi[[g]] <- array(NA, dim(sims_pi$pi)[1:3])
  for(i in 1:dim(sims_pi$pi)[1]){
    stationary_pi[[g]][i,,] <- 
      sims_pi$pi[i,,,g] %*%
      sims_pi$pi[i,,,g] %*%
      sims_pi$pi[i,,,g] %*%
      sims_pi$pi[i,,,g] %*%
      sims_pi$pi[i,,,g] %*%
      sims_pi$pi[i,,,g] %*%
      sims_pi$pi[i,,,g] %*%
      sims_pi$pi[i,,,g] %*%
      sims_pi$pi[i,,,g] %*%
      sims_pi$pi[i,,,g] %*%
      sims_pi$pi[i,,,g]
  }
}
head(stationary_pi[[1]])
lapply(stationary_pi, function(x) apply(x, MARGIN = c(2:3), median))
head(sims_pi$pi)

#function to simulate future regimes
get_regime <- function(sample, sims){
  regime <- 
    data.frame(
      year = 1:100,
      rand = runif(100, 0, 1),
      phi_11 = rep(diag(sims[sample,,])[1], 100),
      phi_22 = rep(diag(sims[sample,,])[2], 100),
      regime = c(sample(1:2, 1, prob = c(.5, .5)), rep(NA, 99)))
  
  for(i in 2:100){
    if(regime$regime[i - 1] == 1){
      if(regime$rand[i] <= regime$phi_11[i]){regime$regime[i] = 1}
      else{regime$regime[i] = 2}
    }
    else{
      if(regime$rand[i] <= regime$phi_22[i]){regime$regime[i] = 2}
      else{regime$regime[i] = 1}
    }
  }
  as.integer(regime$regime)
}

#draw lnalpha, beta and sigma and format to join regime dataframe
sims_lnalpha <-
  lapply(1:post_tvRicker_rt$data$G, function(g){
    data.frame(sample = 1:1000,
               lnalpha = post_tvRicker_rt$sims.list$lnalpha[draw, , g]) %>%
    pivot_longer(cols = starts_with("lnalpha"), 
                 names_to = "regime", 
                 names_prefix = "lnalpha.", 
                 values_to = "lnalpha") %>%
    mutate(regime = as.integer(regime))
  })

sims_beta <-
  lapply(1:post_tvRicker_rt$data$G, function(g){
    data.frame(sample = 1:1000,
               beta = post_tvRicker_rt$sims.list$beta[draw, , g]) %>%
    pivot_longer(cols = starts_with("beta"), 
                 names_to = "regime", 
                 names_prefix = "beta.", 
                 values_to = "beta") %>%
    mutate(regime = as.integer(regime))
  })

sims_gamma <-
  lapply(1:post_tvRicker_rt$data$G, function(g){
    data.frame(sample = 1:1000,
               gamma = post_tvRicker_rt$sims.list$gamma[draw, , g]) %>%
    pivot_longer(cols = starts_with("gamma"), 
                 names_to = "regime", 
                 names_prefix = "gamma.", 
                 values_to = "gamma") %>%
    mutate(regime = as.integer(regime))
  })

sims_sigma <-
  lapply(1:post_tvRicker_rt$data$G, function(g){
    data.frame(sample = 1:1000,
               sigma = post_tvRicker_rt$sims.list$sigma[draw, , g]) %>%
      pivot_longer(cols = starts_with("sigma"), 
                   names_to = "regime", 
                   names_prefix = "sigma.", 
                   values_to = "sigma") %>%
      mutate(regime = as.integer(regime))
  })


#create a dataframe with 100 years of future regimes for each 1000 samples
regimes <-
  lapply(1:post_tvRicker_rt$data$G, function(g){
    lapply(1:1000, function(x){
      data.frame(sample = x, year = 1:100, regime = get_regime(x, stationary_pi[[g]]))
    }) %>%
    do.call(rbind, .) 
  })

#join with appropriate SR parameter estimates to give SR parameter estimates 
# by sample and year after accounting for regime membership 
SRdata <-
  lapply(1:post_tvRicker_rt$data$G, function(g){
    regimes[[g]] %>%
      left_join(sims_lnalpha[[g]], by = c("sample", "regime")) %>%
      left_join(sims_beta[[g]], by = c("sample", "regime")) %>%
      left_join(sims_gamma[[g]], by = c("sample", "regime")) %>%
      left_join(sims_sigma[[g]], by = c("sample", "regime")) %>%
      mutate(epsilon = exp(rnorm(n(), 0, sigma)))
  })

#Calculate yield and recruitment for a wide range of S across samples and years
#R_epsilon and Y_epsilon used to estimate EG range using simulation and GLMM fit.
#Rp, Yp, Smsy_p, MSY_p, OYP_90_p for mean based OYP
#R, Y, Smsy, MSY, OYP_90 for median based OYP
est_Y <-
  lapply(SRdata, function(x){
    x %>%
    slice(rep(1:n(), each = 100)) %>%
    mutate(lnalpha_p = lnalpha + sigma * sigma / 2,
           S = rep(seq(400, 40000, length.out = 100), times = 100 * 1000), #2),
           R = SRgamma(exp(lnalpha), beta, gamma, S),
           #Rp = S * exp(lnalpha_p - beta  * S),
           R_epsilon = R * epsilon,
           Y = R - S,
           #Yp = Rp - S,
           Y_epsilon = R_epsilon - S#,
           #Smsy = lnalpha / (beta) * (0.5 - 0.07 * lnalpha),
           #Smsy_p = lnalpha_p / (beta) * (0.5 - 0.07 * lnalpha_p),
           #MSY = Smsy * exp(lnalpha - beta * Smsy) - Smsy,
           #MSY_p = Smsy_p * exp(lnalpha - beta * Smsy_p) - Smsy_p,
           #OYP_90 = ifelse(Y > 0.9 * MSY, 1, 0),
           #OYP_90_p = ifelse(Yp > 0.9 * MSY_p, 1, 0),
           #Smax = 1 / (beta),
           #MSR = Smax * exp(lnalpha - beta * Smax),
           #ORP_90 = ifelse(R > 0.9 * MSR, 1, 0)
           )
  })


# Simulation / GLMM eg ranges
#build dataset
# Seems like this should be a weighted average to account for the very skewed range of 
# escapements we are likely to observed under each regime.
# I'm not sure how to do that but here is a very rough approach
sim_fullRicker <- 
  sim_Ricker(1.5, 
             0.0001, 
             sigW = 0.5, 
             phi = 0, 
             age0 = c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02), 
             Sims0 = 100000, 
             Hfun = H_goal,
             lb_goal = lb_pctMSY,
             ub_goal = ub_pctMSY,
             sigF = 0.4,
             sigN = 0.0) %>%
  mutate(S_group = cut(S, 
                       breaks = seq(0, 40000, length.out = 101), 
                       labels = seq(400, 40000, length.out = 100))) %>%
  group_by(S_group) %>%
  summarise(n_full = n())
sim_redRicker <- 
  sim_Ricker(lnalpha * 0.25, 
             0.0001, 
             sigW = 0.5, 
             phi = 0, 
             age0 = c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02), 
             Sims0 = 100000, 
             Hfun = H_goal,
             lb_goal = lb_pctMSY,
             ub_goal = ub_pctMSY,
             sigF = 0.4,
             sigN = 0.0) %>%
  mutate(S_group = cut(S, 
                       breaks = seq(0, 40000, length.out = 101), 
                       labels = seq(400, 40000, length.out = 100))) %>%
  group_by(S_group) %>%
  summarise(n_reduced = n())
weights <-
  sim_fullRicker %>% 
  full_join(sim_redRicker, by = "S_group") %>%
  mutate(n_full = ifelse(is.na(n_full), 0, n_full),
         n_reduced = ifelse(is.na(n_reduced), 0, n_reduced),
         weight_2 = n_reduced / (n_full + n_reduced),
         weight_1 = 1 - weight_2,
         S = as.numeric(as.character(S_group))) %>%
  pivot_longer(cols = starts_with("weight"), 
               names_to = "regime", 
               names_pattern = ".*_(\\d)", 
               values_to = "weight") %>%
  mutate(regime = as.numeric(regime))

sim_Y <-
  lapply(est_Y, function(x){
    x %>%
    mutate(Y_epsilon_clean = 
             ifelse(Y_epsilon > quantile(x$Y_epsilon, probs = 0.99), 
                    quantile(x$Y_epsilon, probs = 0.99),
                    Y_epsilon)) %>%
    group_by(regime, S) %>%
    summarise(Y_epsilon_mean = mean(Y_epsilon_clean)) %>%
    left_join(weights, by = c("regime", "S")) %>%
    mutate(weighted_summand = Y_epsilon_mean * weight) %>%
    group_by(S) %>%
    summarise(Y_epsilon_wmean = sum(weighted_summand) / sum(weight))
  })

#  * Smsy and eg range (mean)
mod_Ymean <- lapply(sim_Y, function(x) gam(Y_epsilon_wmean ~ s(S), data = x))
fit_Ymean <- lapply(mod_Ymean, 
                    function(x){
                      predict(x, newdata = data.frame(S = 1:40000)) %>% 
                        round()
                    })
#Smsy
Smsy_tvRicker_rt <- sapply(fit_Ymean, function(x) as.numeric(names(x[which.max(x)])))
#eg bounds based on 90% and 70% of MSY
eg_tvRicker_rt <- 
    sapply(fit_Ymean, function(x){
      lb <- x[x > round(max(x) * 0.9)] %>% 
      names() %>%
      as.numeric() %>%
      min()
      
      ub <- x[x > round(max(x) * 0.7)] %>% 
        names() %>%
        as.numeric() %>%
        max()
      c(lb, ub)
    }) %>%
  t() %>%
  as.data.frame() %>%
  setNames(c("lb", "ub")) %>%
  mutate(group = 1:length(fit_Ymean),
         Smsy = Smsy_tvRicker_rt)

ref_tvRicker_rt <- 
  est_tvRicker_rt %>%
  left_join(eg_tvRicker_rt, by = "group") %>%
  rowwise() %>%
  mutate(set = #(g-1)/d will also do it
           optimize(
             function(a, b, g, x){
               -log(SRgamma(alpha = a, beta = b, gamma = g, S = x) / x)
             },
             interval = c(10, Smsy),
             #maximum = TRUE,
             a = a,
             b = b,
             g = g)$minimum)

# ** plot -------------
# Can we estimate the gamma parameter when the data comes from a time-varying Ricker 
# using the regime transition model?
# Mostly works.
# The estimated goal accounts for both regimes and I'm not sure that is appropriate. Maybe 
# a better approach is to set goals based on the "normal" years and a recovery goal of some type 
# during the conservation concern years.
# Trying to take beta mostly from the "normal" regime.
# gamma coming only from the low regime.
#pull_groups <- sample(1:max(data_tvRicker$group), 3, replace = FALSE)
rbind(estimated_tvRicker_rt, assumed_tvRicker_rt) %>%
  filter(group %in% pull_groups) %>%
  ggplot(aes(x = S, y = value, color = source, shape = regime)) +
  geom_line() +
  geom_point(data = data_tvRicker_rt[data_tvRicker_rt$group %in% pull_groups, ]) +
  scale_x_continuous(limits = c(0, 20000), n.breaks = 5) +
  geom_rect(data = ref_tvRicker_rt[ref_tvRicker_rt$group %in% pull_groups, ], 
            aes(xmin = lb, xmax = ub, ymin = -Inf, ymax = Inf), 
            inherit.aes = FALSE, alpha= 0.2) +
  annotate(geom = "rect", 
           xmin = lb_pctMSY, xmax = ub_pctMSY, ymin = -Inf, ymax = Inf, 
           alpha = 0.2, fill = "orange") +
  scale_color_manual(values = c("orange", "grey")) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(data = ref_tvRicker_rt[ref_tvRicker_rt$group %in% pull_groups, ], 
             aes(xintercept = set),
             linetype = 2) +
  geom_vline(data = ref_tvRicker_rt[ref_tvRicker_rt$group %in% pull_groups, ], 
             aes(xintercept = Smsy), 
             linetype = 3) +
  facet_grid(stat ~ group, scales = "free_y")

#Demonstration of how the weighting is poor
#In essesence there is not yield to weight under regime 2
Ymean <- est_Y[[10]] %>%
  mutate(Y_epsilon_clean = 
           ifelse(Y_epsilon > quantile(est_Y[[1]]$Y_epsilon, probs = 0.99), 
                  quantile(est_Y[[1]]$Y_epsilon, probs = 0.99),
                  Y_epsilon)) %>%
  group_by(regime, S) %>%
  summarise(Y_epsilon_mean = mean(Y_epsilon_clean))

Yw <- est_Y[[10]] %>%
  mutate(Y_epsilon_clean = 
           ifelse(Y_epsilon > quantile(est_Y[[1]]$Y_epsilon, probs = 0.99), 
                  quantile(est_Y[[1]]$Y_epsilon, probs = 0.99),
                  Y_epsilon)) %>%
  group_by(regime, S) %>%
  summarise(Y_epsilon_mean = mean(Y_epsilon_clean)) %>%
  left_join(weights, by = c("regime", "S")) %>%
  mutate(weighted_summand = Y_epsilon_mean * weight)

Ywmean <- est_Y[[10]] %>%
  mutate(Y_epsilon_clean = 
           ifelse(Y_epsilon > quantile(est_Y[[1]]$Y_epsilon, probs = 0.99), 
                  quantile(est_Y[[1]]$Y_epsilon, probs = 0.99),
                  Y_epsilon)) %>%
  group_by(regime, S) %>%
  summarise(Y_epsilon_mean = mean(Y_epsilon_clean)) %>%
  left_join(weights, by = c("regime", "S")) %>%
  mutate(weighted_summand = Y_epsilon_mean * weight) %>%
  group_by(S) %>%
  summarise(Y_epsilon_wmean = sum(weighted_summand) / sum(weight))

plot(Ymean$S, Ymean$Y_epsilon_mean, col = Ymean$regime)
points(Yw$S, Yw$weighted_summand, col = Yw$regime, pch = 2)
lines(Ywmean$S, Ywmean$Y_epsilon_wmean)
#end of demonstration

# compare parameter estimates
data.frame(gamma_Ricker = post_Ricker$q50$gamma,
           gamma_tvRicker = post_tvRicker$q50$gamma,
           gamma_tvRicker_rt = post_tvRicker_rt$q50$gamma[2, ],
           gamma_gamma = post_gamma$q50$gamma) %>%
  pivot_longer(cols = starts_with("gamma"), 
               names_to = "SR", 
               names_pattern = "gamma_(.*)",
               values_to = "gamma") %>%
  ggplot(aes(x = SR, y = gamma)) +
  geom_boxplot()

data.frame(beta_Ricker = post_Ricker$q50$beta,
           beta_tvRicker = post_tvRicker$q50$beta,
           beta_tvRicker_rt1 = post_tvRicker_rt$q50$beta[1, ],
           beta_tvRicker_rt2 = post_tvRicker_rt$q50$beta[2, ],
           beta_gamma = post_gamma$q50$beta) %>%
  pivot_longer(cols = starts_with("beta"), 
               names_to = "SR", 
               names_pattern = "beta_(.*)",
               values_to = "beta") %>%
  ggplot(aes(x = SR, y = beta)) +
  geom_boxplot()

data.frame(lnalpha_Ricker = post_Ricker$q50$lnalpha,
           lnalpha_tvRicker = post_tvRicker$q50$lnalpha,
           lnalpha_tvRicker_rt1 = post_tvRicker_rt$q50$lnalpha[1, ],
           lnalpha_tvRicker_rt2 = post_tvRicker_rt$q50$lnalpha[2, ],
           lnalpha_gamma = post_gamma$q50$lnalpha) %>%
  pivot_longer(cols = starts_with("lnalpha"), 
               names_to = "SR", 
               names_pattern = "lnalpha_(.*)",
               values_to = "lnalpha") %>%
  ggplot(aes(x = SR, y = lnalpha)) +
  geom_boxplot()

# Set compared to lb and Smsy
ref_tvRicker_rt %>% 
  mutate(pct_Smsy = set / Smsy, 
         pct_lb = set / lb) %>% 
  ggplot(aes(x = regime, y = pct_Smsy)) + 
    geom_boxplot()

ref_tvRicker_rt %>% 
  mutate(pct_Smsy = set / Smsy, 
         pct_lb = set / lb) %>% 
  ggplot(aes(x = regime, y = pct_lb)) + 
  geom_boxplot()

