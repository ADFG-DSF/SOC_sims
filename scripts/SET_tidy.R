# SET threshold

# Author: Adam Reimer
# Version: 2025-02-13

# Packages
packs <- c("tidyverse", "jagsUI", "mgcv", "ggpubr", "flextable")
lapply(packs, require, character.only = TRUE)

# source functions
function_files <- list.files(path=".\\functions")
lapply(function_files, function(x) source(paste0(".\\functions\\", x)))

#The concept of using the peak of the ln(R/S) curve from a estimated gamma SR explored in
# the SET.R script seems like it might work. This script is using a tidy workflow so it is 
# different parameter combinations.
# --------------------------------------------------------------------------------------
# motivation
# Simulate a dataset with time varying productivity to demonstrate 
# connection between S and our inability to estimate certain parameters.
# Also demonstrate how a depensatory pattern can arise from non-depensatory dynamics
# when time varying productivity is present.
  
sim_tv <- 
  sim_Ricker(c(rep(1.75, 300), rep(1.5, 300), rep(1.25, 300), rep(1, 300), rep(0.75, 300), rep(0.375, 300)),
             0.0001,
             sigW = 0.25,
             phi = 0,
             age0 = c('3' = 0.1, '4' = 0.38, '5' = 0.3, '6' = 0.2, '7' = 0.02),
             Sims0 = 1800,
             Hfun = H_goal,
             lb_goal = 3893,
             ub_goal = 10529,
             power = 0.6,
             sigF = 0,
             sigN = 0) %>%
  mutate(lnRS = log(R/S))
#plot the data and the data generating model
colors <- scales::hue_pal()(6)
plot_simtv <- 
  sim_tv %>% 
  ggplot(aes(x = S, y = lnRS, color = as.character(lnalpha))) + 
  geom_point() +
  stat_function(fun = function(x) log(Ricker(1.75, 0.0001, x) / x),
                linewidth = 1.25,
                color = colors[6]) +
  stat_function(fun = function(x) log(Ricker(1.5, 0.0001, x) / x),
                linewidth = 1.25,
                color = colors[5]) +
  stat_function(fun = function(x) log(Ricker(1.25, 0.0001, x) / x),
                linewidth = 1.25,
                color = colors[4]) +
  stat_function(fun = function(x) log(Ricker(1, 0.0001, x) / x),
                linewidth = 1.25,
                color = colors[3]) +
  stat_function(fun = function(x) log(Ricker(0.75, 0.0001, x) / x),
                linewidth = 1.25,
                color = colors[2]) +
  stat_function(fun = function(x) log(Ricker(0.375, 0.0001, x) / x),
                linewidth = 1.25,
                color = colors[1]) +
  theme_bw() +
  annotate("rect", xmin = 3893, xmax = 10529, ymin = -Inf, ymax = Inf, fill = 'gray', alpha = .4)
plot_simtv
# Smooth through all of the data
# Temporal varibility destroy our ability to estimate beta.
# You could argue we are seeing depensation but it's masked by lack of compensation.
plot_simtv +
  stat_smooth(aes(x = S, y = lnRS),
              se = FALSE,
              inherit.aes = FALSE, 
              color = "black")
# if we smooth through the most and least productive halves of the data we can 
# estimate compensation and depensation. Notice depensation is not intrinsic to the stock 
# but is an artifact of our observation process. I'm arguing this does not matter
# wrt management.
plot_simtv +
  stat_smooth(data = sim_tv[sim_tv$lnalpha >= 1.25, ], 
              aes(x = S, y = lnRS),
              se = FALSE,
              inherit.aes = FALSE, 
              color = "black") +
  stat_smooth(data = sim_tv[sim_tv$lnalpha < 1.25, ],
              aes(x = S, y = lnRS),
              se = FALSE,
              inherit.aes = FALSE, 
              color = "black")

# Base behavior ------------------------------------
SET_point <-
  expand.grid(lnalpha = c(1, 1.5, 2),
              beta = c(0.001, 0.0001, 0.00001, 0.000001),
              sigma = c(0.25, 0.5, 0.75),
              gamma = seq(1, 1.6, 0.1),
              pct_MSY = c(0.5, 0.6, 0.7, 0.8, 0.9)) %>%
  mutate(Smsy_Ricker = get_Smsy(lnalpha, beta),
         Smax = 1 / beta,
         Rmax = Ricker(lnalpha, beta, Smax),
         a = gamma_par(Smax, Rmax, gamma)[[1]],
         b = gamma_par(Smax, Rmax, gamma)[[2]],
         SET = (gamma - 1) / b) %>%
  rowwise() %>%
  mutate(lb_Ricker = optimise(f = get_bounds, #'true' OYP bounds
                              interval = 1:Smsy_Ricker,
                              lnalpha = lnalpha,
                              beta = beta,
                              pct_MSY = pct_MSY,
                              correct = FALSE)$minimum,
         ub_Ricker = optimise(f = get_bounds, #'true' OYP bounds
                              interval = Smsy_Ricker:(Smsy_Ricker*3),
                              lnalpha = lnalpha,
                              beta = beta,
                              pct_MSY = 0.7,
                              correct = FALSE)$minimum,
         Smsy_gamma =   
           optimize(
             function(a, b, g, x){
               SRgamma(alpha = a, beta = b, gamma = g, S = x) - x
             },
             interval = c(0, 4 * 1 / beta),
             maximum = TRUE,
             a = a,
             b = b,
             g = gamma)$maximum,
         MSY_gamma = SRgamma(alpha = a, beta = b, gamma = gamma, S = Smsy_gamma) - Smsy_gamma) %>%
  mutate(lb_gamma = 
           optimize(
             function(a, b, g, pct_MSY, MSY_gamma, x){
               ((SRgamma(alpha = a, beta = b, gamma = g, S = x) - x) - pct_MSY * MSY_gamma)^2
             },
             interval = c(Smsy_gamma / 10, Smsy_gamma),
             a = a,
             b = b,
             g = gamma,
             pct_MSY = pct_MSY,
             MSY = MSY_gamma)$minimum,
         ub_gamma = 
           optimize(
             function(a, b, g, pct_MSY, MSY_gamma, x){
               ((SRgamma(alpha = a, beta = b, gamma = g, S = x) - x) - pct_MSY * MSY_gamma)^2
             },
             interval = c(Smsy_gamma, Smsy_gamma * 3),
             a = a,
             b = b,
             g = gamma,
             pct_MSY = 0.7,
             MSY = MSY_gamma)$minimum)
SET_point %>%
  mutate(pct_Smsy = Smsy_gamma / Smsy_Ricker) %>%
  ggplot(aes(x = gamma, y = pct_Smsy, color = as.character(lnalpha))) +
  geom_line()

SET_point %>%
  mutate(pct_lb = lb_gamma / lb_Ricker) %>%
  ggplot(aes(x = gamma, y = pct_lb, color = as.character(lnalpha))) +
  geom_line() +
  facet_grid(. ~ pct_MSY)

SET_point %>%
  mutate(pct_ub = ub_gamma / ub_Ricker) %>%
  ggplot(aes(x = gamma, y = pct_ub, color = as.character(lnalpha))) +
  geom_line()

SET_point %>%
  mutate(pct_lbRicker = SET / lb_Ricker,
         pct_lbgamma = SET / lb_gamma) %>%
  pivot_longer(cols = c("pct_lbRicker", "pct_lbgamma"), 
               names_to = "Reference", 
               values_to = "Percentage") %>%
  ggplot(aes(x = gamma, y = Percentage, color = as.character(lnalpha), linetype = Reference)) +
    geom_line() +
  geom_hline(aes(yintercept = 1)) +
    facet_grid(. ~ pct_MSY) #beta makes no difference

SET_point %>%
  filter(pct_MSY %in% c(0.5, 0.9),
         gamma %in% c(1, 1.3, 1.6),
         lnalpha == 1.5,
         beta == 0.0001,
         sigma == 0.25) %>%
  mutate(curves = list(data.frame(S = rep(seq(0, 15000, 100), times = 3),
                                  R = SRgamma(a, b, gamma, rep(seq(0, 15000, 100), times = 3))))) %>% 
  unnest(curves) %>%
  ggplot(aes(x = S, y = R)) +
  geom_line() +
  geom_rect(aes(xmin = lb_gamma, xmax = ub_gamma, ymin = -Inf, ymax = Inf), alpha = 0.005, fill = "green") +
  geom_vline(aes(xintercept = SET), color = "red", linetype = 2) +
  geom_abline(aes(intercept = 0, slope = 1))+
  facet_grid(paste0("\u03B3: ", ifelse(gamma == 1, "1 (Ricker)", gamma)) ~ paste0("%MSR @ lb: ", pct_MSY))

# Define scenarios to simulate -------
scenarios_g <-
  expand.grid(lnalpha = c(1, 1.5, 2),
              beta = 0.0001, #c(0.001, 0.0001, 0.000001),
              sigma = c(0.25, 0.5),
              phi = 0,
              pct_lb = c(0.5, 0.9),
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
                              correct = TRUE, #but used the log alpha correction here!
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
  mutate(lnalpha_red = lb_pctMSY * 0.5  * beta,
         Rmax_red = Ricker(lnalpha_red, beta, Smax),
         a_red = gamma_par(Smax, Rmax_red, gamma)[[1]]) %>% 
  ungroup()

#  * Scenarios table ------------------------------------------------------
# for each scenario we have 50 replicate datasets and SR models
scenarios_g %>% 
  select(scenario, gamma, pct_lb, sigma, lnalpha) %>%
  flextable() %>%
  set_header_labels(
    scenario = "Scenario",
    lnalpha = "ln(\u03b1)",
    sigma = "\u03c3",
    gamma = "\u03b3",
    pct_lb = "MSY % @ EG \n lower bound"
  ) %>%
  merge_v(j = c("lnalpha", "sigma", "gamma", "pct_lb")) %>%
  valign(j = c("lnalpha", "sigma", "gamma", "pct_lb"), valign = "top") %>%
  autofit()

#Simulate data and estimate parameters --------
rep_scenarios_gp <-
  scenarios_g %>%
  slice(rep(1:nrow(scenarios_g), each = 50)) %>%
  mutate(rep = rep(1:50, times = nrow(scenarios_g))) %>%
  rowwise() %>%
  mutate(data = list(
    sim_SRgamma(c(rep(a, 40 / 2), rep(a_red, 40 / 2)),
                b,
                gamma = gamma,
                sigW = sigma,
                phi = 0,
                age0 = c('3' = 0.1, '4' = 0.38, '5' = 0.3, '6' = 0.2, '7' = 0.02),
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
                       lnRS = log(.$R / .$S),
                       S = .$S,
                       ar1 = 0)}),
         mod_gamma = map(data_jags, ~ jags(data = .x,
                                           parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
                                           model.file = ".\\scripts\\gammRS_changepoint.txt",
                                           n.chains = 3,
                                           n.iter = 5e4,
                                           n.burnin = 1e4,
                                           n.thin = 480, 
                                           parallel = TRUE))
         # mod_Ricker = map(data_jags, ~ jags(data = .x,
         #                                    parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "lambda"),
         #                                    model.file = ".\\scripts\\Ricker_changepoint.txt",
         #                                    n.chains = 3,
         #                                    n.iter = 5e4,
         #                                    n.burnin = 1e4,
         #                                    n.thin = 480, 
         #                                    parallel = TRUE))
  )
#Simulate data and estimate parameters --------
rep_scenarios_gp2 <-
  rep_scenarios_gp %>%
  mutate(mod_Ricker = map(data_jags, ~ jags(data = .x,
                                            parameters.to.save = c("lnalpha", "beta", "sigma", "y_d", "lambda"),
                                            model.file = ".\\scripts\\Ricker_changepoint.txt",
                                            n.chains = 3,
                                            n.iter = 5e4,
                                            n.burnin = 1e4,
                                            n.thin = 480,
                                            parallel = TRUE)))
saveRDS(rep_scenarios_gp2, file = ".\\rep_scenarios_gp2.rds")
rep_scenarios_gp <- readRDS(file = ".\\rep_scenarios_gp2.rds")

# * Rhat boxplot ---------
# gamma model
# for main parameters (lnalpha, beta, gamma, sigma, y_d)
rep_scenarios_gp %>%
  rowwise() %>%
  mutate(Rhat_gamma = list(mod_gamma$summary[c("lnalpha[1]", "lnalpha[2]",
                                               "beta[1]", "beta[2]",
                                               "gamma[1]", "gamma[2]",
                                               "sigma[1]", "sigma[2]",
                                               "y_d"), "Rhat"])) %>%
  select(lnalpha, beta, sigma, phi, pct_lb, scenario, Rhat_gamma) %>%
  unnest(cols = c(Rhat_gamma)) %>%
  filter(!is.na(Rhat_gamma)) %>%         
  pivot_longer(cols = starts_with("Rhat"), 
               names_to = "model", 
               names_pattern = "Rhat_(.*)",
               values_to = "Rhat") %>%
  ggplot(aes(x = scenario, group = scenario, y = Rhat)) + 
  geom_boxplot() +
  geom_hline(aes(yintercept = 1.01)) +
  coord_cartesian(ylim = c(1, 1.25))

#Ricker model
rep_scenarios_gp %>%
  rowwise() %>%
  mutate(Rhat_Ricker = list(mod_Ricker$summary[c("lnalpha[1]", "lnalpha[2]",
                                               "beta[1]", "beta[2]",
                                               "sigma[1]", "sigma[2]",
                                               "y_d"), "Rhat"])) %>%
  select(lnalpha, beta, sigma, phi, pct_lb, scenario, Rhat_Ricker) %>%
  unnest(cols = c(Rhat_Ricker)) %>%
  filter(!is.na(Rhat_Ricker)) %>%         
  pivot_longer(cols = starts_with("Rhat"), 
               names_to = "model", 
               names_pattern = "Rhat_(.*)",
               values_to = "Rhat") %>%
  ggplot(aes(x = scenario, group = scenario, y = Rhat)) + 
  geom_boxplot() +
  geom_hline(aes(yintercept = 1.01)) +
  coord_cartesian(ylim = c(1, 1.25))

#  * SET and EG lb --------------------------------------
# median SET for 50 replicated under each scenario
# black line is the SET from the data generating model
# red star is the lower bound of the escapement goal (50% or 90% of MSY)
# Note that SET > lb for many 50% of MSY scenarios but i think that is extreme...
# The egegik lb is ~ 50% of MSY based on the BBSRI top 20 prior 
# (which gave a lower Smsy than thier mean prior).
rep_scenarios_gp %>%
  select(scenario, rep, lnalpha, lnalpha_red, beta, sigma, gamma, a, a_red, b, pct_lb, lb_pctMSY, mod_gamma) %>%
  mutate(set = (gamma - 1) / b) %>%
  mutate(set_est = map_dbl(mod_gamma, function(x) (x$q50$gamma[2] - 1) / x$q50$beta[2])) %>%
  ggplot(aes(x = scenario, y = set_est, group = scenario)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = set)) +
  geom_point(aes(y = lb_pctMSY), shape = 8, color = "red") +
  facet_grid(. ~ gamma, scales = "free_x")

# Simulated values for the SET in each scenario (all replicates)
# black line is the SET from the data generating model
# dashed line is the lower bound of the escapement goal (50% or 90% of MSY)
rep_scenarios_gp %>%
  mutate(set = (gamma - 1) / b, 
         set_sims = map(mod_gamma, function(x) (x$sims.list$gamma[,2] - 1) / x$sims.list$beta[,2])) %>%
  select(lnalpha, lnalpha_red, beta, sigma, pct_lb, a, a_red, b, gamma, lb_pctMSY, scenario, rep, set, set_sims) %>%
  unnest(set_sims) %>%
  ggplot(aes(x = set_sims, fill = as.character(rep), group = rep)) + 
  geom_histogram() +
  geom_vline(aes(xintercept = set)) +
  geom_vline(aes(xintercept = lb_pctMSY), linetype = 2) +
  scale_x_continuous(limits = c(0, 12000)) +
  facet_wrap(. ~ scenario)

#  * Plot data and fit ----------------------------------------------------
# randomly selects on replicate to display.
# user is responsible for identifying a singe scenario in the dots.
plot_fit <- function(scenario_file, ...){
  filters <- enquos(...)
  
  data0 <- 
    scenario_file %>%
    filter(!!!filters)
  
  rep_var = sample(unique(data0$rep), 1, replace = FALSE)
  
  data <- 
    data0 %>%
    filter(rep == rep_var)
  
  maxS <- round(max(data$data_jags[[1]]$S), -2)
  
  true <- 
    data %>%
    select(scenario, rep, a, a_red, b, gamma, lb_pctMSY, ub_pctMSY) %>%
    mutate(set = (gamma - 1)/ b)
  
  estimates <- 
    data %>%
    select(scenario, rep, mod_gamma) %>%
    mutate(a_1 = map_dbl(mod_gamma, ~exp(.$q50$lnalpha[1])),
           a_2 = map_dbl(mod_gamma, ~exp(.$q50$lnalpha[2])),
           b_1 = map_dbl(mod_gamma, ~.$q50$beta[1]),
           b_2 = map_dbl(mod_gamma, ~.$q50$beta[2]),
           gamma_1 = map_dbl(mod_gamma, ~.$q50$gamma[1]),
           gamma_2 = map_dbl(mod_gamma, ~.$q50$gamma[2]),
           set = (gamma_2 - 1) / b_2)
  
  plot1 <- 
    data %>%
    mutate(lnRS = map(data_jags, function(x) data.frame(year = 1:40, lnRS = x$lnRS, S = x$S)),
           cut = map_dbl(mod_gamma, ~.$q50$y_d)) %>%
    select(scenario, rep, a, a_red, b, gamma, lb_pctMSY, ub_pctMSY, cut, lnRS) %>%
    unnest(lnRS) %>%
    mutate(regime = ifelse(year <= cut, "normal", "reduced")) %>%
    ggplot(aes(x = S, y = lnRS)) +
    geom_point(aes(color = regime)) +
    geom_vline(data = true, aes(xintercept = set), color = "red") +
    geom_vline(data = estimates, aes(xintercept = set), linetype = 2, color = "red") +
    geom_rect(data = true, aes(xmin = lb_pctMSY, xmax = ub_pctMSY, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, 
              alpha= .25) +
    geom_function(fun = function(x) log(SRgamma(true$a, true$b, true$gamma, x) / x) ,
                  color = "black") +
    geom_function(fun = function(x) log(SRgamma(true$a_red, true$b, true$gamma, x) / x),
                  color = "red") +
    geom_function(fun = function(x) log(SRgamma(estimates$a_1, estimates$b_1, estimates$gamma_1, x) / x),
                  color = "black",
                  linetype = 2) +
    geom_function(fun = function(x) log(SRgamma(estimates$a_2, estimates$b_2, estimates$gamma_2, x) / x),
                  col = "red",
                  linetype = 2) +
    scale_color_manual(values = c("black", "red")) +
    scale_x_continuous(limits = c(0, 1.1 * maxS)) +
    theme_bw()
  
  plot2 <- 
    data %>%
    mutate(S = map(data_jags, function(x) data.frame(year = 1:40, S = x$S)),
           cut = map_dbl(mod_gamma, ~.$q50$y_d)) %>%
    unnest(S) %>%
    mutate(regime = ifelse(year <= cut, "normal", "reduced")) %>%
    select(scenario, rep, year, b, gamma, lb_pctMSY, regime, S) %>% 
    mutate(set = (gamma - 1)/ b) %>% 
    ggplot(aes(x = year, y = S)) + 
    geom_line() +
    geom_point(aes(color = regime)) +
    geom_hline(data = true, aes(yintercept = set), color = "red") +
    geom_hline(data = estimates, aes(yintercept = set), linetype = 2, color = "red") +
    geom_rect(data = true, aes(ymin = lb_pctMSY, ymax = ub_pctMSY, xmin = -Inf, xmax = Inf),
              inherit.aes = FALSE, 
              alpha= .25) +
    scale_color_manual(values = c("black", "red")) +
    scale_y_continuous(limits = c(0, 1.1 * maxS)) +
    theme_bw()
  
  ggarrange(plot1, plot2, nrow = 2, heights = c(6,4))
}
plot_fit(rep_scenarios_gp, 
         lnalpha == 2,
         sigma == .5,
         gamma == 1.3,
         pct_lb == 0.5)


# * Horsetail ------------------------------------
plot_horse <- function(scenario_file, scenario){
  coeflines1 <-
    data.frame(lnalpha = exp(scenario_file$mod_gamma[[26]]$sims.list$lnalpha[, 1]),
               beta = scenario_file$mod_gamma[[26]]$sims.list$beta[, 1],
               gamma = scenario_file$mod_gamma[[26]]$sims.list$gamma[, 1]) %>%
    dplyr::sample_n(50) %>%
    as.matrix() %>%
    plyr::alply(1, function(coef) {
      ggplot2::stat_function(fun = function(x) log(coef[1] * x^(coef[3] - 1) * exp(-coef[2]*x)), 
                             colour="grey",
                             alpha = 0.5)})
  coeflines2 <-
    data.frame(lnalpha = exp(scenario_file$mod_gamma[[26]]$sims.list$lnalpha[, 2]),
               beta = scenario_file$mod_gamma[[26]]$sims.list$beta[, 2],
               gamma = scenario_file$mod_gamma[[26]]$sims.list$gamma[, 2]) %>%
    dplyr::sample_n(50) %>%
    as.matrix() %>%
    plyr::alply(1, function(coef) {
      ggplot2::stat_function(fun = function(x) log(coef[1] * x^(coef[3] - 1) * exp(-coef[2]*x)), 
                             colour="grey",
                             alpha = 0.5)})
  
  data.frame(S = scenario_file$data_jags[[17]]$S,
             lnRS = scenario_file$data_jags[[17]]$lnRS) %>%
    ggplot(aes(x = S, y = lnRS)) +
    geom_point() +
    geom_rect(data = scenario_file[scenario_file$scenario == scenario & scenario_file$rep == 1, ],
              aes(xmin = lb_pctMSY, 
                  xmax = ub_pctMSY, 
                  ymin = -Inf, 
                  ymax = Inf),
              inherit.aes = FALSE,
              fill = "grey",
              alpha= .25) +
    theme_bw() +
    coeflines1 + 
    coeflines2
}
plot_horse(rep_scenarios_gp, 20)

# * Plot curves -----------------------
# This was an early version... not sure it's helpful
# Does show assumed (grey) and estimated (orange) goal ranges, SET, and Smsy.
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
    rename(gamma_1_ass = gamma,
           Smsy1 = Smsy) %>%
    mutate(lnalpha_1_ass = log(a), 
           lnalpha_2_ass = log(a_red),
           beta_1_ass = gamma_1_ass * beta,
           beta_2_ass = beta_1_ass,
           gamma_2_ass = gamma_1_ass,
           lnalpha_1_est = map_dbl(mod_gamma, ~ .$q50$lnalpha[1]),
           lnalpha_2_est = map_dbl(mod_gamma, ~ .$q50$lnalpha[2]),
           beta_1_est = map_dbl(mod_gamma, ~ .$q50$beta[1]),
           beta_2_est = map_dbl(mod_gamma, ~ .$q50$beta[2]),
           gamma_1_est = map_dbl(mod_gamma, ~ .$q50$gamma[1]),
           gamma_2_est = map_dbl(mod_gamma, ~ .$q50$gamma[2])) %>%
    select(-lnalpha, -lnalpha_red, -beta) %>%
    select(scenario, Smsy1, pct_lb, lb_pctMSY, ub_pctMSY, starts_with(c("lnalpha", "beta", "gamma"))) %>% 
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
           MSY = SRgamma(exp(lnalpha), beta, gamma, Smsy) - Smsy,
           set = ifelse(regime == 2 & gamma >= 1, (gamma - 1) / beta, NA))
  
  eg_bounds <- function(S, lnalpha, beta, gamma, pct_MSY, MSY){
    ((SRgamma(exp(lnalpha), beta, gamma, S) - S) - pct_MSY * MSY)^2
  }
  
  bounds <-
    ref_points %>%
    filter(regime == 1, source =="est") %>%
    mutate(lb = optimise(f = eg_bounds, #'true' OYP bounds
                         interval = 1:Smsy,
                         lnalpha = lnalpha,
                         beta = beta,
                         gamma = gamma,
                         MSY = MSY,
                         pct_MSY = pct_lb)$minimum,
           ub = optimise(f = eg_bounds, #'true' OYP bounds
                         interval = Smsy:(Smsy*5),
                         lnalpha = lnalpha,
                         beta = beta,
                         gamma = gamma,
                         MSY = MSY,
                         pct_MSY = 0.7)$minimum) %>%
    select(lb, ub, scenario)
  
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
              inherit.aes = FALSE, alpha= 0.002, fill = "orange") +
    geom_rect(data = bounds,
              aes(xmin = lb, xmax = ub, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, alpha= .2) +
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
    ggh4x::facetted_pos_scales(y = scales) +
    theme_bw()
}
plot_sim(rep_scenarios_gp, 
         sigma == 0.5, beta == 0.0001, lnalpha == 2, gamma == 1.3, 
         facet_var = pct_lb, 
         rep_var = NULL)

# How often did the correct model have the lower dic?
#Ricker model
rep_scenarios_gp %>%
  rowwise() %>%
  mutate(dic_Ricker = mod_Ricker$DIC,
         dic_gamma = mod_gamma$DIC) %>%
  select(lnalpha, beta, sigma, phi, pct_lb, scenario, starts_with("dic")) %>%
  mutate(model = ifelse(dic_Ricker < dic_gamma, "Ricker", "gamma")) %>%
  ggplot(aes(x = scenario, fill = model)) + 
  geom_bar() +
  geom_hline(aes(yintercept = 25))

# SET stable across variation in lnalpha and gamma?
data.frame(lnalpha = rep_scenarios_gp$mod_gamma[[36]]$sims.list$lnalpha[,2], 
           gamma = rep_scenarios_gp$mod_gamma[[36]]$sims.list$gamma[,2],
           beta = rep_scenarios_gp$mod_gamma[[36]]$sims.list$beta[,2]) %>%
  mutate(set = (gamma - 1) / beta,
         set_cut = (cut(set, quantile(set, c(0, 0.05, 0.2, 0.8, 0.95, 1)))),
         beta_cut = (cut(beta, quantile(beta, c(0, 0.05, 0.2, 0.8, 0.95, 1))))) %>%
  ggplot(aes(x = lnalpha, y = gamma, color = set_cut)) +
  geom_point()
