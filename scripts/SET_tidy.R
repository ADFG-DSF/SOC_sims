# SET threshold

# Author: Adam Reimer
# Version: 2025-02-13

# Packages
packs <- c("tidyverse", "jagsUI", "mgcv", "ggpubr")
lapply(packs, require, character.only = TRUE)

# source functions
function_files <- list.files(path=".\\functions")
lapply(function_files, function(x) source(paste0(".\\functions\\", x)))

#The concept of using the peak of the ln(R/S) curve from a estimated gamma SR explored in
# the SET.R script seems like it might work. This script is using a tidy workflow so its 
# different parameter combinations.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Ricker SR -------
scenarios_R <-
  expand.grid(lnalpha = c(1, 1.5, 2),
              beta = 0.0001, #c(0.001, 0.0001, 0.000001),
              sigma = c(0.25, 0.5, 0.75),
              phi = 0,
              pct_lb = c(0.5, 0.7, 0.9)) %>%
  mutate(scenario = 1:n(),
         Smsy = get_Smsy(lnalpha, beta, correct = FALSE), #notice no log alpha correction for right now.
         power = lnalpha / 2 - .2) %>% # very crude. 
  rowwise() %>%
  mutate(lb_pctMSY = optimise(f = get_bounds, #'true' OYP bounds
                              interval = 1:get_Smsy(lnalpha, beta),
                              lnalpha = lnalpha,
                              beta = beta,
                              pct_MSY = pct_lb,
                              correct = TRUE,
                              sigma = sigma,
                              phi = phi)$minimum,
         ub_pctMSY = optimise(f = get_bounds, #'true' OYP bounds
                              interval = get_Smsy(lnalpha, beta):(get_Smsy(lnalpha, beta)*5),
                              lnalpha = lnalpha,
                              beta = beta,
                              pct_MSY = 0.7,
                              correct = TRUE,
                              sigma = sigma,
                              phi = phi)$minimum) %>%
  mutate(lnalpha_red = lb_pctMSY * beta) %>% 
  ungroup ()

# Replicate each scenario and simulate SR data
rep_scenarios_Ricker <-
  scenarios_R %>%
  slice(rep(1:nrow(scenarios), each = 20)) %>%
  mutate(rep = rep(1:20, times = nrow(scenarios))) %>%
  rowwise() %>%
  mutate(data = list(
           sim_Ricker(c(rep(lnalpha, 40 / 2), rep(lnalpha_red, 40 / 2)),
                      beta,
                      sigW = sigma,
                      phi = 0,
                      age0 = c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02),
                      Sims0 = 40,
                      Hfun = H_goal,
                      lb_goal = lb_pctMSY,
                      ub_goal = ub_pctMSY,
                      power = power,
                      sigF = 0.2,
                      sigN = 0.2))) %>%
  ungroup() %>%
  mutate(data_jags =
           map(data, ~.x %>%
                 {list(nyrs = max(.$sim),
                       S = .$S,
                       R = .$R,
                       ar1 = 0)}),
         mod_gamma = map(data_jags, ~ jags(data = .x,
                                           parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
                                           model.file = ".\\scripts\\SRgamma_changepoint.txt",
                                           n.chains = 3,
                                           n.iter = 5e4,
                                           n.burnin = 1e4,
                                           n.thin = 480, 
                                           parallel = TRUE)),
         mod_Ricker = map(data_jags, ~ jags(data = .x,
                                            parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
                                            model.file = ".\\scripts\\Ricker_changepoint.txt",
                                            n.chains = 3,
                                            n.iter = 5e4,
                                            n.burnin = 1e4,
                                            n.thin = 480, 
                                            parallel = TRUE))
  )
#saveRDS(rep_scenarios_Ricker, file = ".\\rep_scenarios_Ricker.rds")
rep_scenarios <- readRDS(file = ".\\rep_scenarios.rds")

# * Plot results -------
plot_sim <- function(scenario_file, ..., facet_var, rep_var = NULL){
  if(is.null(rep_var) == TRUE){rep_var = sample(unique(scenario_file$rep), 1, replace = FALSE)}
  
  filters <- enquos(...)
  
  scenarios <-
    scenario_file %>% 
    filter(!!!filters) %>% 
    select(scenario) %>% 
    unlist(use.names = FALSE)
  
  SR_data <- 
    scenario_file %>%
    filter(rep == rep_var, 
           scenario %in% scenarios) %>%
    select(scenario, data) %>%
    unnest(data) %>%
    mutate(regime = ifelse(sim <= 20, "1", "2"),
           source = "ass") %>%
    select(scenario, regime, source, S, R) %>%
    mutate(Y = R - S,
           lnRS = log(R /S)) %>%
    pivot_longer(cols = c(R, Y, lnRS), names_to = "param", values_to = "value")
  
  range <- 
    SR_data %>% 
    group_by(scenario) %>% 
    summarise(max_S = max(S))
  
  temp <- 
    scenario_file %>%
    filter(rep == rep_var, 
           scenario %in% scenarios) %>%
    rename(lnalpha_1_ass = lnalpha, 
           lnalpha_2_ass = lnalpha_red,
           beta_1_ass = beta,
           Smsy1 = Smsy) %>%
    mutate(beta_2_ass = beta_1_ass,
           gamma_1_ass = 1,
           gamma_2_ass = 1,
           lnalpha_1_est = map_dbl(mod, ~ .$q50$lnalpha[1]),
           lnalpha_2_est = map_dbl(mod, ~ .$q50$lnalpha[2]),
           beta_1_est = map_dbl(mod, ~ .$q50$beta[1]),
           beta_2_est = map_dbl(mod, ~ .$q50$beta[2]),
           gamma_1_est = map_dbl(mod, ~ .$q50$gamma[1]),
           gamma_2_est = map_dbl(mod, ~ .$q50$gamma[2])) %>%
    select(scenario, Smsy1, lb_pctMSY, ub_pctMSY, starts_with(c("lnalpha", "beta", "gamma"))) %>% 
    left_join(range, by = c("scenario")) %>%
    pivot_longer(cols = starts_with(c("lnalpha", "beta", "gamma")), 
                 names_to = c("param", "regime", "source"),
                 names_pattern = "(.*)_(.*)_(.*)",
                 values_to = "values") %>%
    pivot_wider(names_from = "param",
                values_from = "values")
  
  curves <- 
    temp %>%
    rowwise() %>%
    mutate(curves = list(
      data.frame(S = seq(0, max_S, by = max_S / 99)) %>%
        mutate(R = SRgamma(alpha = exp(lnalpha), beta = beta, gamma = gamma, S = S),
               Y = R - S,
               lnRS = log(R / S))))
  
  ref_points <- 
    temp %>%
    rowwise() %>%
    mutate(Smsy2 =   
             optimize(
               function(a, b, g, x){
                 SRgamma(alpha = exp(a), beta = b, gamma = g, S = x) - x
               },
               interval = c(0, 4 * 1 / beta),
               maximum = TRUE,
               a = lnalpha,
               b = beta,
               g = gamma)$maximum,
           Smsy = ifelse(regime == 1, ifelse(source == "ass", Smsy1, Smsy2), NA),
           set = ifelse(regime == 2 & gamma >= 1, (gamma - 1) / beta, NA))
  
  facet_name <- deparse(substitute(facet_var))
  labels0 <- 
    scenario_file %>% 
    filter(rep == rep_var,
           !!!filters) %>% 
    select(scenario, {{facet_var}}) %>% 
    mutate(label = paste0(facet_name, ": ", {{facet_var}}))
  label_facet <- labels0$label
  names(label_facet) <- labels0$scenario
  
  facet_scales <- 
    SR_data %>% 
    group_by(param) %>% 
    summarise(ymin = min(value),
              ymax = max(value),
              y95 = quantile(value, 0.95)) %>% 
    mutate(ymin = ifelse(param == "R", 0, ymin),
           ymax = ifelse(param == "lnRS", ymax, y95),
           n = 5) %>%
    ungroup() %>%
    rename(Panel = param)
  facet_scales <- split(facet_scales, facet_scales$Panel)
  
  scales <- lapply(facet_scales, function(x) {
    scale_y_continuous(limits = c(x$ymin, x$ymax), n.breaks = x$n)
  })
  
  curves %>%
    unnest(curves) %>%
    pivot_longer(cols = c(R, Y, lnRS), names_to = "param", values_to = "value") %>%
    ggplot(aes(x = S, y = value, color = source, shape = regime)) +
    geom_line() +
    geom_point(data = SR_data) +
    geom_rect(aes(xmin = lb_pctMSY, xmax = ub_pctMSY, ymin = -Inf, ymax = Inf), 
              inherit.aes = FALSE, alpha= 0.002) +
    scale_color_manual(values = c("orange", "black")) +
    geom_hline(aes(yintercept = value), 
               data.frame(param = c("lnRS", "R", "Y"), value = c(0, 0, 0)), color = "grey60") +
    geom_abline(aes(intercept = int, slope = slope), 
               data.frame(param = "R", int = 0, slope = 1), color = "grey60") +
    geom_vline(data = ref_points, 
               aes(xintercept = set, color = source),
               linetype = 2) +
    geom_vline(data = ref_points, 
               aes(xintercept = Smsy, color = source), 
               linetype = 3) +
    facet_grid(param ~ scenario, scales = "free", labeller = labeller(.cols = label_facet)) +
    ggh4x::facetted_pos_scales(y = scales)
}

plot_sim(rep_scenarios_Ricker, 
         sigma == 0.5, beta == 0.0001, lnalpha == 2, 
         facet_var = pct_lb, 
         rep_var = NULL)

# average % of gamma posterior > 1 across each replicate
rep_scenarios_Ricker %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario, rep) %>%
  mutate(pct_dep = map_dbl(mod_gamma, ~ mean(.$sims.list$gamma[, 2] > 1))) %>%
  ggplot(aes(x = lnalpha, group = scenario, y = pct_dep)) + 
  geom_boxplot() +
  facet_grid(sigma ~ pct_lb)

# average % of beta posterior > 0 across each replicate
rep_scenarios_Ricker %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario) %>%
  mutate(pct_dep = map_dbl(mod_gamma, ~ mean(.$sims.list$beta[, 1] > 0))) %>%
  ggplot(aes(x = scenario, group = scenario, y = pct_dep)) + 
  geom_boxplot()

# Rhat boxplot for main parameters (lnalpha, beta, gamma, sigma, y_d)
rep_scenarios_Ricker %>%
  rowwise() %>%
  mutate(Rhat_gamma = list(mod_gamma$summary[c("lnalpha[1]", "lnalpha[2]",
                                               "beta[1]", "beta[2]",
                                               "gamma[1]", "gamma[2]",
                                               "sigma[1]", "sigma[2]",
                                               "y_d"), "Rhat"]),
         Rhat_Ricker = list(mod_Ricker$summary[c("lnalpha[1]", "lnalpha[2]",
                                                 "beta[1]", "beta[2]",
                                                 "gamma[1]", "gamma[2]",
                                                 "sigma[1]", "sigma[2]",
                                                 "y_d"), "Rhat"])) %>%
  select(lnalpha, beta, sigma, phi, pct_lb, scenario, Rhat_gamma, Rhat_Ricker) %>%
  unnest(cols = c(Rhat_gamma, Rhat_Ricker)) %>%
  filter(!is.na(Rhat_Ricker), !is.na(Rhat_gamma)) %>%
  pivot_longer(cols = starts_with("Rhat"), 
               names_to = "model", 
               names_pattern = "Rhat_(.*)",
               values_to = "Rhat") %>%
  ggplot(aes(x = scenario, group = scenario, y = Rhat)) + 
  geom_boxplot() +
  facet_grid(. ~ model)

# Percent of replicates where the gamma is the preferred SR relationship
rep_scenarios_Ricker %>%
  mutate(DIC_gamma = map_dbl(mod_gamma, ~ .$DIC),
         DIC_Ricker = map_dbl(mod_Ricker, ~ .$DIC),
         DIC_logical = DIC_gamma < DIC_Ricker) %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario) %>%
  summarise(DIC_mean = mean(DIC_logical)) %>%
  ggplot(aes(x = lnalpha, group = scenario, y = DIC_mean)) + 
  geom_point() +
  facet_grid(sigma ~ pct_lb)











rep_scenarios %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario, rep) %>%
  mutate(pct_dep = map_dbl(mod, ~ mean(.$sims.list$gamma[, 2] > 1))) %>%
  ggplot(aes(x = lnalpha, group = scenario, y = pct_dep, 
             color = as.character(beta))) + 
    geom_boxplot() +
    facet_grid(sigma ~ pct_lb)

rep_scenarios %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario) %>%
  mutate(pct_dep = map_dbl(mod, ~ .$q50$gamma[2])) %>%
  ggplot(aes(x = scenario, group = scenario, y = pct_dep)) + 
  geom_boxplot()

rep_scenarios %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario) %>%
  mutate(pct_dep = map_dbl(mod, ~ mean(.$sims.list$beta[, 1] > 0))) %>%
  ggplot(aes(x = scenario, group = scenario, y = pct_dep)) + 
  geom_boxplot()

rep_scenarios %>%
  rowwise() %>%
  mutate(Rhat = list(unlist(mod$Rhat))) %>%
  select(lnalpha, beta, sigma, phi, pct_lb, scenario, Rhat) %>%
  unnest(Rhat) %>%
  filter(!is.na(Rhat)) %>%
  ggplot(aes(x = scenario, group = scenario, y = Rhat)) + 
  geom_boxplot()












# Gamma SR -------
scenarios_g <-
  expand.grid(lnalpha = c(1, 1.5, 2),
              beta = 0.0001, #c(0.001, 0.0001, 0.000001),
              sigma = c(0.25, 0.5, 0.75),
              phi = 0,
              pct_lb = c(0.5, 0.7, 0.9),
              gamma = seq(1, 1.6, length.out = 3)) %>%
  mutate(scenario = 1:n(),
         Smsy = get_Smsy(lnalpha, beta, correct = FALSE), #notice no log alpha correction for right now.
         power = lnalpha / 2 - .2,
         Smax = 1/ beta,
         Rmax = Ricker(lnalpha, beta, Smax),
         a = gamma_par(Smax, Rmax, gamma)[[1]],
         b = gamma_par(Smax, Rmax, gamma)[[2]]) %>% # very crude. 
  rowwise() %>%
  mutate(lb_pctMSY = optimise(f = get_bounds, #'true' OYP bounds
                              interval = 1:get_Smsy(lnalpha, beta),
                              lnalpha = lnalpha,
                              beta = beta,
                              pct_MSY = pct_lb,
                              correct = TRUE,
                              sigma = sigma,
                              phi = phi)$minimum,
         ub_pctMSY = optimise(f = get_bounds, #'true' OYP bounds
                              interval = get_Smsy(lnalpha, beta):(get_Smsy(lnalpha, beta)*5),
                              lnalpha = lnalpha,
                              beta = beta,
                              pct_MSY = 0.7,
                              correct = TRUE,
                              sigma = sigma,
                              phi = phi)$minimum) %>%
  mutate(lnalpha_red = lb_pctMSY * beta,
         Rmax_red = Ricker(lnalpha_red, beta, Smax),
         a_red = gamma_par(Smax, Rmax_red, gamma)[[1]]) %>% 
  ungroup ()

# Replicate each scenario and simulate SR data
rep_scenarios_gamma <-
  scenarios_g %>%
  slice(rep(1:nrow(scenarios), each = 20)) %>%
  mutate(rep = rep(1:20, times = nrow(scenarios))) %>%
  rowwise() %>%
  mutate(data = list(
    sim_SRgamma(c(rep(a, 40 / 2), rep(a_red, 40 / 2)),
               b,
               gamma = gamma,
               sigW = sigma,
               phi = 0,
               age0 = c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02),
               Sims0 = 40,
               Hfun = H_goal,
               lb_goal = lb_pctMSY,
               ub_goal = ub_pctMSY,
               power = power,
               sigF = 0.2,
               sigN = 0.2))) %>%
  ungroup() %>%
  mutate(data_jags =
           map(data, ~.x %>%
                 {list(nyrs = max(.$sim),
                       S = .$S,
                       R = .$R,
                       ar1 = 0)}),
         mod_gamma = map(data_jags, ~ jags(data = .x,
                                           parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
                                           model.file = ".\\scripts\\SRgamma_changepoint.txt",
                                           n.chains = 3,
                                           n.iter = 5e4,
                                           n.burnin = 1e4,
                                           n.thin = 480, 
                                           parallel = TRUE)),
         mod_Ricker = map(data_jags, ~ jags(data = .x,
                                     parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
                                     model.file = ".\\scripts\\Ricker_changepoint.txt",
                                     n.chains = 3,
                                     n.iter = 5e4,
                                     n.burnin = 1e4,
                                     n.thin = 480, 
                                     parallel = TRUE))
  )
#saveRDS(rep_scenarios_gamma, file = ".\\rep_scenarios_gamma.rds")
rep_scenarios_gamma <- readRDS(file = ".\\rep_scenarios_gamma.rds")


# * Plot results ----------
plot_sim(rep_scenarios_gamma, 
         sigma == 0.5, beta == 0.0001, lnalpha == 2, 
         facet_var = pct_lb, 
         rep_var = NULL)

# average % of gamma posterior > 1 across each replicate
rep_scenarios_gamma %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario, rep) %>%
  mutate(pct_dep = map_dbl(mod_gamma, ~ mean(.$sims.list$gamma[, 2] > 1))) %>%
  ggplot(aes(x = lnalpha, group = scenario, y = pct_dep, 
             color = as.character(gamma))) + 
  geom_boxplot() +
  facet_grid(sigma ~ pct_lb)

# average % of beta posterior > 0 across each replicate
rep_scenarios_gamma %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario) %>%
  mutate(pct_dep = map_dbl(mod_gamma, ~ mean(.$sims.list$beta[, 1] > 0))) %>%
  ggplot(aes(x = scenario, group = scenario, y = pct_dep)) + 
  geom_boxplot()

# Rhat boxplot for main parameters (lnalpha, beta, gamma, sigma, y_d)
rep_scenarios_gamma %>%
  rowwise() %>%
  mutate(Rhat_gamma = list(mod_gamma$summary[c("lnalpha[1]", "lnalpha[2]",
                                               "beta[1]", "beta[2]",
                                               "gamma[1]", "gamma[2]",
                                               "sigma[1]", "sigma[2]",
                                               "y_d"), "Rhat"]),
         Rhat_Ricker = list(mod_Ricker$summary[c("lnalpha[1]", "lnalpha[2]",
                                                 "beta[1]", "beta[2]",
                                                 "gamma[1]", "gamma[2]",
                                                 "sigma[1]", "sigma[2]",
                                                 "y_d"), "Rhat"])) %>%
  select(lnalpha, beta, sigma, phi, pct_lb, scenario, Rhat_gamma, Rhat_Ricker) %>%
  unnest(cols = c(Rhat_gamma, Rhat_Ricker)) %>%
  filter(!is.na(Rhat_Ricker), !is.na(Rhat_gamma)) %>%         
  pivot_longer(cols = starts_with("Rhat"), 
               names_to = "model", 
               names_pattern = "Rhat_(.*)",
               values_to = "Rhat") %>%
  ggplot(aes(x = scenario, group = scenario, y = Rhat)) + 
  geom_boxplot() +
  facet_grid(. ~ model)

# Percent of replicates where the gamma is the preferred SR relationship
rep_scenarios_gamma %>%
  mutate(DIC_gamma = map_dbl(mod_gamma, ~ .$DIC),
         DIC_Ricker = map_dbl(mod_Ricker, ~ .$DIC),
         DIC_logical = DIC_gamma < DIC_Ricker) %>%
  group_by(lnalpha, beta, sigma, phi, gamma, pct_lb, scenario) %>%
  summarise(DIC_mean = mean(DIC_logical)) %>%
  ggplot(aes(x = lnalpha, group = scenario, y = DIC_mean, 
             color = as.character(gamma))) + 
    geom_point() +
    facet_grid(sigma ~ pct_lb)




# Try changes in F---------
scenarios_RF <-
  expand.grid(lnalpha = c(1, 1.5, 2),
              beta = 0.0001, #c(0.001, 0.0001, 0.000001),
              sigma = c(0.25, 0.5, 0.75),
              phi = 0,
              pct_lb = c(0.5, 0.7, 0.9),
              sigF = c(0.2, 0.4)) %>%
  mutate(scenario = 1:n(),
         Smsy = get_Smsy(lnalpha, beta, correct = FALSE), #notice no log alpha correction for right now.
         power = lnalpha / 2 - .2) %>% # very crude. 
  rowwise() %>%
  mutate(lb_pctMSY = optimise(f = get_bounds, #'true' OYP bounds
                              interval = 1:get_Smsy(lnalpha, beta),
                              lnalpha = lnalpha,
                              beta = beta,
                              pct_MSY = pct_lb,
                              correct = TRUE,
                              sigma = sigma,
                              phi = phi)$minimum,
         ub_pctMSY = optimise(f = get_bounds, #'true' OYP bounds
                              interval = get_Smsy(lnalpha, beta):(get_Smsy(lnalpha, beta)*5),
                              lnalpha = lnalpha,
                              beta = beta,
                              pct_MSY = 0.7,
                              correct = TRUE,
                              sigma = sigma,
                              phi = phi)$minimum) %>%
  mutate(lnalpha_red = lb_pctMSY * beta) %>% 
  ungroup ()

# Replicate each scenario and simulate SR data
rep_scenarios_RF <-
  scenarios_RF %>%
  slice(rep(1:nrow(scenarios_RF), each = 20)) %>%
  mutate(rep = rep(1:20, times = nrow(scenarios_RF))) %>%
  rowwise() %>%
  mutate(data = list(
    sim_Ricker(c(rep(lnalpha, 40 / 2), rep(lnalpha_red, 40 / 2)),
               beta,
               sigW = sigma,
               phi = 0,
               age0 = c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02),
               Sims0 = 40,
               Hfun = H_goal,
               lb_goal = lb_pctMSY,
               ub_goal = ub_pctMSY,
               power = power,
               sigF = sigF,
               sigN = 0.2))) %>%
  ungroup() %>%
  mutate(data_jags =
           map(data, ~.x %>%
                 {list(nyrs = max(.$sim),
                       S = .$S,
                       R = .$R,
                       ar1 = 0)}),
         mod_gamma = map(data_jags, ~ jags(data = .x,
                                           parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
                                           model.file = ".\\scripts\\SRgamma_changepoint.txt",
                                           n.chains = 3,
                                           n.iter = 5e4,
                                           n.burnin = 1e4,
                                           n.thin = 480, 
                                           parallel = TRUE)),
         mod_Ricker = map(data_jags, ~ jags(data = .x,
                                            parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
                                            model.file = ".\\scripts\\Ricker_changepoint.txt",
                                            n.chains = 3,
                                            n.iter = 5e4,
                                            n.burnin = 1e4,
                                            n.thin = 480, 
                                            parallel = TRUE))
  )
saveRDS(rep_scenarios_RF, file = ".\\rep_scenarios_RF.rds")
#rep_scenarios <- readRDS(file = ".\\rep_scenarios.rds")

scenarios_gF <-
  expand.grid(lnalpha = c(1, 1.5, 2),
              beta = 0.0001, #c(0.001, 0.0001, 0.000001),
              sigma = 0.5, #c(0.25, 0.5, 0.75),
              phi = 0,
              pct_lb = c(0.5, 0.7, 0.9),
              gamma = seq(1, 1.6, length.out = 3),
              sigF = c(0.2, 0.4)) %>%
  mutate(scenario = 1:n(),
         Smsy = get_Smsy(lnalpha, beta, correct = FALSE), #notice no log alpha correction for right now.
         power = lnalpha / 2 - .2,
         Smax = 1/ beta,
         Rmax = Ricker(lnalpha, beta, Smax),
         a = gamma_par(Smax, Rmax, gamma)[[1]],
         b = gamma_par(Smax, Rmax, gamma)[[2]]) %>% # very crude. 
  rowwise() %>%
  mutate(lb_pctMSY = optimise(f = get_bounds, #'true' OYP bounds
                              interval = 1:get_Smsy(lnalpha, beta),
                              lnalpha = lnalpha,
                              beta = beta,
                              pct_MSY = pct_lb,
                              correct = TRUE,
                              sigma = sigma,
                              phi = phi)$minimum,
         ub_pctMSY = optimise(f = get_bounds, #'true' OYP bounds
                              interval = get_Smsy(lnalpha, beta):(get_Smsy(lnalpha, beta)*5),
                              lnalpha = lnalpha,
                              beta = beta,
                              pct_MSY = 0.7,
                              correct = TRUE,
                              sigma = sigma,
                              phi = phi)$minimum) %>%
  mutate(lnalpha_red = lb_pctMSY * beta,
         Rmax_red = Ricker(lnalpha_red, beta, Smax),
         a_red = gamma_par(Smax, Rmax_red, gamma)[[1]]) %>% 
  ungroup ()

# Replicate each scenario and simulate SR data
rep_scenarios_gF <-
  scenarios_gF %>%
  slice(rep(1:nrow(scenarios_gF), each = 20)) %>%
  mutate(rep = rep(1:20, times = nrow(scenarios_gF))) %>%
  rowwise() %>%
  mutate(data = list(
    sim_SRgamma(c(rep(a, 40 / 2), rep(a_red, 40 / 2)),
                b,
                gamma = gamma,
                sigW = sigma,
                phi = 0,
                age0 = c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02),
                Sims0 = 40,
                Hfun = H_goal,
                lb_goal = lb_pctMSY,
                ub_goal = ub_pctMSY,
                power = power,
                sigF = sigF,
                sigN = 0.2))) %>%
  ungroup() %>%
  mutate(data_jags =
           map(data, ~.x %>%
                 {list(nyrs = max(.$sim),
                       S = .$S,
                       R = .$R,
                       ar1 = 0)}),
         mod_gamma = map(data_jags, ~ jags(data = .x,
                                           parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
                                           model.file = ".\\scripts\\SRgamma_changepoint.txt",
                                           n.chains = 3,
                                           n.iter = 5e4,
                                           n.burnin = 1e4,
                                           n.thin = 480, 
                                           parallel = TRUE)),
         mod_Ricker = map(data_jags, ~ jags(data = .x,
                                            parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
                                            model.file = ".\\scripts\\Ricker_changepoint.txt",
                                            n.chains = 3,
                                            n.iter = 5e4,
                                            n.burnin = 1e4,
                                            n.thin = 480, 
                                            parallel = TRUE))
  )
saveRDS(rep_scenarios_gF, file = ".\\rep_scenarios_gF.rds")
rep_scenarios_gF <- readRDS(file = ".\\rep_scenarios_gF.rds")



# * Plot results ----------
plot_sim(rep_scenarios_gF, 
         sigma == 0.5, beta == 0.0001, lnalpha == 2, 
         facet_var = pct_lb, 
         rep_var = NULL)

# average % of gamma posterior > 1 across each replicate
rep_scenarios_gF %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario, rep) %>%
  mutate(pct_dep = map_dbl(mod_gamma, ~ mean(.$sims.list$gamma[, 2] > 1))) %>%
  ggplot(aes(x = gamma, group = scenario, y = pct_dep, 
             color = as.character(lnalpha))) + 
  geom_boxplot() +
  facet_grid(sigF ~ pct_lb)

# average % of beta posterior > 0 across each replicate
rep_scenarios_gF %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario) %>%
  mutate(pct_dep = map_dbl(mod_gamma, ~ mean(.$sims.list$beta[, 1] > 0))) %>%
  ggplot(aes(x = scenario, group = scenario, y = pct_dep)) + 
  geom_boxplot()

# Rhat boxplot for main parameters (lnalpha, beta, gamma, sigma, y_d)
rep_scenarios_gF %>%
  rowwise() %>%
  mutate(Rhat_gamma = list(mod_gamma$summary[c("lnalpha[1]", "lnalpha[2]",
                                               "beta[1]", "beta[2]",
                                               "gamma[1]", "gamma[2]",
                                               "sigma[1]", "sigma[2]",
                                               "y_d"), "Rhat"]),
         Rhat_Ricker = list(mod_Ricker$summary[c("lnalpha[1]", "lnalpha[2]",
                                                 "beta[1]", "beta[2]",
                                                 "gamma[1]", "gamma[2]",
                                                 "sigma[1]", "sigma[2]",
                                                 "y_d"), "Rhat"])) %>%
  select(lnalpha, beta, sigma, phi, pct_lb, scenario, Rhat_gamma, Rhat_Ricker) %>%
  unnest(cols = c(Rhat_gamma, Rhat_Ricker)) %>%
  filter(!is.na(Rhat_Ricker), !is.na(Rhat_gamma)) %>%         
  pivot_longer(cols = starts_with("Rhat"), 
               names_to = "model", 
               names_pattern = "Rhat_(.*)",
               values_to = "Rhat") %>%
  ggplot(aes(x = scenario, group = scenario, y = Rhat)) + 
  geom_boxplot() +
  facet_grid(. ~ model)

# Percent of replicates where the gamma is the preferred SR relationship
rep_scenarios_gF %>%
  mutate(DIC_gamma = map_dbl(mod_gamma, ~ .$DIC),
         DIC_Ricker = map_dbl(mod_Ricker, ~ .$DIC),
         DIC_logical = DIC_gamma < DIC_Ricker) %>%
  group_by(lnalpha, beta, sigma, phi, gamma, pct_lb, sigF, scenario) %>%
  summarise(DIC_mean = mean(DIC_logical)) %>%
  ggplot(aes(x = gamma, group = scenario, y = DIC_mean, 
             color = as.character(lnalpha))) +
  geom_point() +
  facet_grid(sigF ~ pct_lb)


# Try changes in productivity---------
# Replicate each scenario and simulate SR data
rep_scenarios_Rp <-
  scenarios_R %>%
  slice(rep(1:nrow(scenarios_R), each = 20)) %>%
  mutate(rep = rep(1:20, times = nrow(scenarios_R))) %>%
  rowwise() %>%
  mutate(data = list(
    sim_Ricker(c(rep(lnalpha, 40 / 2), rep(lnalpha_red, 40 / 4), rep(lnalpha_red * 0.7, 40 / 4)),
               beta,
               sigW = sigma,
               phi = 0,
               age0 = c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02),
               Sims0 = 40,
               Hfun = H_goal,
               lb_goal = lb_pctMSY,
               ub_goal = ub_pctMSY,
               power = power,
               sigF = 0.2,
               sigN = 0.2))) %>%
  ungroup() %>%
  mutate(data_jags =
           map(data, ~.x %>%
                 {list(nyrs = max(.$sim),
                       S = .$S,
                       R = .$R,
                       ar1 = 0)}),
         mod_gamma = map(data_jags, ~ jags(data = .x,
                                           parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
                                           model.file = ".\\scripts\\SRgamma_changepoint.txt",
                                           n.chains = 3,
                                           n.iter = 5e4,
                                           n.burnin = 1e4,
                                           n.thin = 480, 
                                           parallel = TRUE)),
         mod_Ricker = map(data_jags, ~ jags(data = .x,
                                            parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
                                            model.file = ".\\scripts\\Ricker_changepoint.txt",
                                            n.chains = 3,
                                            n.iter = 5e4,
                                            n.burnin = 1e4,
                                            n.thin = 480, 
                                            parallel = TRUE))
  )
#saveRDS(rep_scenarios_Rp, file = ".\\rep_scenarios_Rp.rds")
rep_scenarios_Rp <- readRDS(file = ".\\rep_scenarios_Rp.rds")

# * Plot results ----------
plot_sim(rep_scenarios_Rp, 
         sigma == 0.5, beta == 0.0001, lnalpha == 2, 
         facet_var = pct_lb, 
         rep_var = NULL)

# average % of gamma posterior > 1 across each replicate
rep_scenarios_Rp %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario, rep) %>%
  mutate(pct_dep = map_dbl(mod_gamma, ~ mean(.$sims.list$gamma[, 2] > 1))) %>%
  ggplot(aes(x = lnalpha, group = scenario, y = pct_dep)) + 
  geom_boxplot() +
  facet_grid(sigma ~ pct_lb)

# average % of beta posterior > 0 across each replicate
rep_scenarios_Rp %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario) %>%
  mutate(pct_dep = map_dbl(mod_gamma, ~ mean(.$sims.list$beta[, 1] > 0))) %>%
  ggplot(aes(x = scenario, group = scenario, y = pct_dep)) + 
  geom_boxplot()

# Rhat boxplot for main parameters (lnalpha, beta, gamma, sigma, y_d)
rep_scenarios_Rp %>%
  rowwise() %>%
  mutate(Rhat_gamma = list(mod_gamma$summary[c("lnalpha[1]", "lnalpha[2]",
                                               "beta[1]", "beta[2]",
                                               "gamma[1]", "gamma[2]",
                                               "sigma[1]", "sigma[2]",
                                               "y_d"), "Rhat"]),
         Rhat_Ricker = list(mod_Ricker$summary[c("lnalpha[1]", "lnalpha[2]",
                                                 "beta[1]", "beta[2]",
                                                 "gamma[1]", "gamma[2]",
                                                 "sigma[1]", "sigma[2]",
                                                 "y_d"), "Rhat"])) %>%
  select(lnalpha, beta, sigma, phi, pct_lb, scenario, Rhat_gamma, Rhat_Ricker) %>%
  unnest(cols = c(Rhat_gamma, Rhat_Ricker)) %>%
  filter(!is.na(Rhat_Ricker), !is.na(Rhat_gamma)) %>%         
  pivot_longer(cols = starts_with("Rhat"), 
               names_to = "model", 
               names_pattern = "Rhat_(.*)",
               values_to = "Rhat") %>%
  ggplot(aes(x = scenario, group = scenario, y = Rhat)) + 
  geom_boxplot() +
  facet_grid(. ~ model)

# Percent of replicates where the gamma is the preferred SR relationship
rep_scenarios_Rp %>%
  mutate(DIC_gamma = map_dbl(mod_gamma, ~ .$DIC),
         DIC_Ricker = map_dbl(mod_Ricker, ~ .$DIC),
         DIC_logical = DIC_gamma < DIC_Ricker) %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario) %>%
  summarise(DIC_mean = mean(DIC_logical)) %>%
  ggplot(aes(x = lnalpha, group = scenario, y = DIC_mean)) + 
  geom_col() +
  geom_hline(aes(yintercept = 0.5)) +
  facet_grid(sigma ~ pct_lb)



# Replicate each scenario and simulate SR data
rep_scenarios_gp <-
  scenarios_g %>%
  slice(rep(1:nrow(scenarios_g), each = 20)) %>%
  mutate(rep = rep(1:20, times = nrow(scenarios_g))) %>%
  rowwise() %>%
  mutate(data = list(
    sim_SRgamma(c(rep(a, 40 / 2), rep(a_red, 40 / 4), rep(a_red * 0.7, 40 / 4)),
                b,
                gamma = gamma,
                sigW = sigma,
                phi = 0,
                age0 = c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02),
                Sims0 = 40,
                Hfun = H_goal,
                lb_goal = lb_pctMSY,
                ub_goal = ub_pctMSY,
                power = power,
                sigF = 0.2,
                sigN = 0.2))) %>%
  ungroup() %>%
  mutate(data_jags =
           map(data, ~.x %>%
                 {list(nyrs = max(.$sim),
                       S = .$S,
                       R = .$R,
                       ar1 = 0)}),
         mod_gamma = map(data_jags, ~ jags(data = .x,
                                           parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
                                           model.file = ".\\scripts\\SRgamma_changepoint.txt",
                                           n.chains = 3,
                                           n.iter = 5e4,
                                           n.burnin = 1e4,
                                           n.thin = 480, 
                                           parallel = TRUE)),
         mod_Ricker = map(data_jags, ~ jags(data = .x,
                                            parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
                                            model.file = ".\\scripts\\Ricker_changepoint.txt",
                                            n.chains = 3,
                                            n.iter = 5e4,
                                            n.burnin = 1e4,
                                            n.thin = 480, 
                                            parallel = TRUE))
  )
#saveRDS(rep_scenarios_gp, file = ".\\rep_scenarios_gp.rds")
rep_scenarios_gp <- readRDS(file = ".\\rep_scenarios_gp.rds")

# Rhat boxplot for main parameters (lnalpha, beta, gamma, sigma, y_d)
rep_scenarios_gp %>%
  rowwise() %>%
  mutate(Rhat_gamma = list(mod_gamma$summary[c("lnalpha[1]", "lnalpha[2]",
                                               "beta[1]", "beta[2]",
                                               "gamma[1]", "gamma[2]",
                                               "sigma[1]", "sigma[2]",
                                               "y_d"), "Rhat"]),
         Rhat_Ricker = list(mod_Ricker$summary[c("lnalpha[1]", "lnalpha[2]",
                                                 "beta[1]", "beta[2]",
                                                 "gamma[1]", "gamma[2]",
                                                 "sigma[1]", "sigma[2]",
                                                 "y_d"), "Rhat"])) %>%
  select(lnalpha, beta, sigma, phi, pct_lb, scenario, Rhat_gamma, Rhat_Ricker) %>%
  unnest(cols = c(Rhat_gamma, Rhat_Ricker)) %>%
  filter(!is.na(Rhat_Ricker), !is.na(Rhat_gamma)) %>%         
  pivot_longer(cols = starts_with("Rhat"), 
               names_to = "model", 
               names_pattern = "Rhat_(.*)",
               values_to = "Rhat") %>%
  ggplot(aes(x = scenario, group = scenario, y = Rhat)) + 
  geom_boxplot() +
  geom_hline(aes(yintercept = 1.01)) +
  coord_cartesian(ylim = c(1, 1.25)) +
  facet_grid(. ~ model)

# 
rep_scenarios_gp %>%
  mutate(lnalpha_est_R1 = map_dbl(mod_Ricker, ~ .$q50$lnalpha[1]),
         lnalpha_est_R2 = map_dbl(mod_Ricker, ~ .$q50$lnalpha[2]),
         lnalpha_est_g1 = map_dbl(mod_gamma, ~ .$q50$lnalpha[1]),
         lnalpha_est_g2 = map_dbl(mod_gamma, ~ .$q50$lnalpha[2])) %>%
  pivot_longer(cols = starts_with("lnalpha_est"), 
               names_to = c("model", "regime"), 
               names_pattern = "lnalpha_est_([Rg])([12])") %>%
  select(scenario, rep, lnalpha, lnalpha_red, a, a_red, gamma, model, regime, value) %>%
  filter(model == "g", regime == 1) %>%
  ggplot(aes(x = log(a),  y = value, group = a)) + 
  geom_jitter() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  facet_grid(gamma ~ .)

rep_scenarios_gp %>%
  mutate(gamma_est_R1 = map_dbl(mod_Ricker, ~ .$q50$gamma[1]),
         gamma_est_R2 = map_dbl(mod_Ricker, ~ .$q50$gamma[2]),
         gamma_est_g1 = map_dbl(mod_gamma, ~ .$q50$gamma[1]),
         gamma_est_g2 = map_dbl(mod_gamma, ~ .$q50$gamma[2])) %>%
  pivot_longer(cols = starts_with("gamma_est"), 
               names_to = c("model", "regime"), 
               names_pattern = "gamma_est_([Rg])([12])") %>%
  select(scenario, rep, lnalpha, lnalpha_red, a, a_red, gamma, model, regime, value) %>%
  filter(model == "g", regime == 2) %>%
  ggplot(aes(x = gamma,  y = value, group = gamma)) + 
  geom_jitter() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  facet_grid(gamma ~ .)

rep_scenarios_gp %>%
  mutate(lnalpha_est_g1 = map_dbl(mod_gamma, ~ .$q50$lnalpha[1]),
         lnalpha_est_g2 = map_dbl(mod_gamma, ~ .$q50$lnalpha[2]),
         gamma_est_g1 = map_dbl(mod_gamma, ~ .$q50$gamma[1]),
         gamma_est_g2 = map_dbl(mod_gamma, ~ .$q50$gamma[2])) %>%
  select(scenario, rep, lnalpha, lnalpha_red, a, a_red, gamma, lnalpha_est_g1, lnalpha_est_g2, gamma_est_g1, gamma_est_g2)
  mutate(gamma_est_R1 = map_dbl(mod_Ricker, ~ .$q50$gamma[1]),
         gamma_est_R2 = map_dbl(mod_Ricker, ~ .$q50$gamma[2]),
         gamma_est_g1 = map_dbl(mod_gamma, ~ .$q50$gamma[1]),
         gamma_est_g2 = map_dbl(mod_gamma, ~ .$q50$gamma[2])) %>%
  pivot_longer(cols = starts_with("gamma_est"), 
               names_to = c("model", "regime"), 
               names_pattern = "gamma_est_([Rg])([12])") %>%
  select(scenario, rep, lnalpha, lnalpha_red, a, a_red, gamma, model, regime, value) %>%
  filter(model == "g", regime == 2) %>%
  ggplot(aes(x = gamma,  y = value, group = gamma)) + 
  geom_jitter() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  facet_grid(gamma ~ .)

# * Plot results ----------
plot_sim(rep_scenarios_gp, 
         sigma == 0.5, beta == 0.0001, lnalpha == 2, 
         facet_var = pct_lb, 
         rep_var = NULL)

# average % of gamma posterior > 1 across each replicate
rep_scenarios_gp %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario, rep) %>%
  mutate(pct_dep = map_dbl(mod_gamma, ~ mean(.$sims.list$gamma[, 2] > 1))) %>%
  ggplot(aes(x = gamma, group = scenario, y = pct_dep, 
             color = as.character(lnalpha))) + 
  geom_boxplot() +
  facet_grid(sigma ~ pct_lb)

# average % of beta posterior > 0 across each replicate
rep_scenarios_gp %>%
  group_by(lnalpha, beta, sigma, phi, pct_lb, scenario) %>%
  mutate(pct_dep = map_dbl(mod_gamma, ~ mean(.$sims.list$beta[, 1] > 0))) %>%
  ggplot(aes(x = scenario, group = scenario, y = pct_dep)) + 
  geom_boxplot()



# Percent of replicates where the gamma is the preferred SR relationship
rep_scenarios_gp %>%
  mutate(DIC_gamma = map_dbl(mod_gamma, ~ .$DIC),
         DIC_Ricker = map_dbl(mod_Ricker, ~ .$DIC),
         DIC_logical = DIC_gamma < DIC_Ricker) %>%
  group_by(lnalpha, beta, sigma, phi, gamma, pct_lb, scenario) %>%
  summarise(DIC_mean = mean(DIC_logical)) %>%
  ggplot(aes(x = gamma, group = scenario, y = DIC_mean, 
             color = as.character(lnalpha))) +
  geom_point() +
  facet_grid(sigma ~ pct_lb)
