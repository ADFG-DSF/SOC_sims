# SET threshold

# Author: Adam Reimer
# Version: 2025-02-13

# Packages
packs <- c("tidyverse", "mgcv", "ggpubr", "flextable", "ggforce")
lapply(packs, require, character.only = TRUE)

# source functions
function_files <- list.files(path=".\\functions")
lapply(function_files, function(x) source(paste0(".\\functions\\", x)))

#The concept of using the peak of the ln(R/S) curve from a estimated gamma SR explored in
# the SET.R script seems like it might work. This script is using a tidy workflow so it is 
# different parameter combinations.
# --------------------------------------------------------------------------------------

# Gamma SR behavior ------------------------------------
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
             MSY = MSY_gamma)$minimum,
         dep_gamma = 
           optimize(
             function(a, b, g, x){
               ((SRgamma(alpha = a, beta = b, gamma = g, S = x) - x)^2)
             },
             interval = c(1, SET),
             a = a,
             b = b,
             g = gamma)$minimum)


#  * ref points vrs Smax ---------------------------------------------------------
pal <- scales::hue_pal()(5)
col_values <- c("0.9" = pal[1], "0.8" = pal[2], "0.7" = pal[3], "0.6" = pal[4], "0.5" = pal[5])
plot_Smax <- 
  SET_point %>%
  mutate(pct_Smsy = Smsy_gamma / Smax,
         pct_lb = lb_gamma / Smax,
         pct_set = SET / Smax,
         pct_dep = dep_gamma / Smax) %>%
  pivot_longer(c(pct_set, pct_lb, pct_Smsy, pct_dep), #pct_ub, 
               names_to = c("stat", "metric"), names_sep = "_",
               values_to = "pct") 
plot_Smax %>%
  filter(metric == "lb") %>%
  ggplot(aes(x = gamma, y = pct, color = as.character(pct_MSY), linetype = "lb")) +
  geom_line(linewidth = 1) +
  geom_line(aes(linetype = "set"), data = plot_Smax[plot_Smax$metric == "set", ], color = "black", linewidth = 1) +
  geom_line(aes(linetype = "Smsy"), data = plot_Smax[plot_Smax$metric == "Smsy", ], color = "black", linewidth = 1) +
  geom_line(aes(linetype = "dep"), data = plot_Smax[plot_Smax$metric == "dep", ], color = "black", linewidth = 2) +
  scale_color_manual(name = "Prop. MSY @ EG lb", values = col_values) +
  scale_linetype_manual(name = "Metric", 
                        breaks = c("dep", "lb", "set", "Smsy"),
                        labels = c("dep"," EG lower bound", "SET", bquote(S[msy])),
                        values = c("dep" = 5, "Smsy" = 3, "set" = 1, "lb" = 2)) +
  scale_y_continuous(name = bquote("Proportion of S"[max])) +
  facet_grid(. ~ paste0("ln(", "\u03B1", "): ",lnalpha))


#  * ref points vrs lb at 90%MSY ---------------------------------------------
plot_lb90 <- 
  SET_point %>%
  filter(pct_MSY == 0.9) %>%
  mutate(pct_Smsy = Smsy_gamma / lb_gamma,
         pct_lb = lb_gamma / lb_gamma,
         pct_set = SET / lb_gamma,
         pct_dep = dep_gamma / lb_gamma) %>%
  pivot_longer(c(pct_set, pct_lb, pct_Smsy, pct_dep), #pct_ub, 
               names_to = c("stat", "metric"), names_sep = "_",
               values_to = "pct") 
plot_lb90 %>%
  filter(metric == "lb") %>%
  ggplot(aes(x = gamma, y = pct, color = as.character(pct_MSY), linetype = "lb")) +
  geom_line(linewidth = 1) +
  geom_line(aes(linetype = "set"), data = plot_lb90[plot_lb90$metric == "set", ], color = "black", linewidth = 1) +
  geom_line(aes(linetype = "Smsy"), data = plot_lb90[plot_lb90$metric == "Smsy", ], color = "black", linewidth = 1) +
  geom_line(aes(linetype = "dep"), data = plot_lb90[plot_lb90$metric == "dep", ], color = "black", linewidth = 2) +
  scale_color_manual(name = "Prop. MSY @ EG lb", values = col_values) +
  scale_linetype_manual(name = "Metric", 
                        breaks = c("dep", "lb", "set", "Smsy"),
                        labels = c("dep"," EG lower bound", "SET", bquote(S[msy])),
                        values = c("dep" = 5, "Smsy" = 3, "set" = 1, "lb" = 2)) +
  scale_y_continuous(name = bquote("Proportion of EG lb")) +
  facet_grid(. ~ paste0("ln(", "\u03B1", "): ",lnalpha))


#  * ref points vrs Smsy ---------------------------------------------------------
plot_Smsy <- 
  SET_point %>%
  mutate(pct_Smsy = Smsy_gamma / Smsy_gamma,
         pct_lb = lb_gamma / Smsy_gamma,
         pct_set = SET / Smsy_gamma,
         pct_dep = dep_gamma / Smsy_gamma) %>%
  pivot_longer(c(pct_set, pct_lb, pct_Smsy, pct_dep), #pct_ub, 
               names_to = c("stat", "metric"), names_sep = "_",
               values_to = "pct") 
plot_Smsy %>%
  filter(metric == "lb") %>%
  ggplot(aes(x = gamma, y = pct, color = as.character(pct_MSY), linetype = "lb")) +
  geom_line(linewidth = 1) +
  geom_line(aes(linetype = "set"), data = plot_Smsy[plot_Smsy$metric == "set", ], color = "black", linewidth = 1) +
  geom_line(aes(linetype = "Smsy"), data = plot_Smsy[plot_Smsy$metric == "Smsy", ], color = "black", linewidth = 1) +
  geom_line(aes(linetype = "dep"), data = plot_Smsy[plot_Smsy$metric == "dep", ], color = "black", linewidth = 2) +
  scale_color_manual(name = "Prop. MSY @ EG lb", values = col_values) +
  scale_linetype_manual(name = "Metric", 
                        breaks = c("dep", "lb", "set", "Smsy"),
                        labels = c("dep"," EG lower bound", "SET", bquote(S[msy])),
                        values = c("dep" = 5, "Smsy" = 3, "set" = 1, "lb" = 2)) +
  scale_y_continuous(name = bquote("Proportion of S"[msy])) +
  facet_grid(. ~ paste0("ln(", "\u03B1", "): ", lnalpha))

plot_Smsy %>%
  filter(metric == "lb") %>%
  filter(gamma %in% c(1, 1.3, 1.6)) %>%
  ggplot(aes(x = lnalpha, y = pct, color = as.character(pct_MSY), linetype = "lb")) +
  geom_line(linewidth = 1) +
  geom_line(aes(linetype = "set"), 
            data = plot_Smsy[plot_Smsy$metric == "set" & plot_Smsy$gamma %in% c(1, 1.3, 1.6), ], 
            color = "black", 
            linewidth = 1) +
  geom_line(aes(linetype = "Smsy"), 
            data = plot_Smsy[plot_Smsy$metric == "Smsy" & plot_Smsy$gamma %in% c(1, 1.3, 1.6), ], 
            color = "black", 
            linewidth = 1) +
  geom_line(aes(linetype = "dep"), 
            data = plot_Smsy[plot_Smsy$metric == "dep" & plot_Smsy$gamma %in% c(1, 1.3, 1.6), ], 
            color = "black", 
            linewidth = 2) +
  scale_color_manual(name = "Prop. MSY @ EG lb", values = col_values) +
  scale_linetype_manual(name = "Metric", 
                        breaks = c("dep", "lb", "set", "Smsy"),
                        labels = c("dep"," EG lower bound", "SET", bquote(S[msy])),
                        values = c("dep" = 5, "Smsy" = 3, "set" = 1, "lb" = 2)) +
  scale_y_continuous(name = bquote("Proportion of S"[msy])) +
  facet_grid(. ~ paste0("\u03B3: ", gamma))


#  * SET vrs. lb ----------------------------------------------------------
# SET vrs. 90%MSY lb for Ricker and Gamma goals
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

# Example goal ranges and SETs for 3 levels of gamma and 2 levels of pct_MSY
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






# Simulations -------
# * Define scenarios to simulate --------
scenarios <-
  expand.grid(lnalpha_1 = c(1, 1.5, 2),
              beta = 0.0001, #c(0.001, 0.0001, 0.000001),
              sigma = c(0.6),
              phi = c(0, 0.5),
              gamma = seq(1, 1.6, length.out = 3),
              harvest_rate = c(0.18, 0.39),
              sigN = c(0.05, 0.4)) %>%
  arrange(lnalpha_1, sigma, gamma, harvest_rate, sigN) %>%
  mutate(scenario = 1:n(),
         Smax = 1/ beta,
         Rmax = Ricker(lnalpha_1, beta, Smax),
         a_1 = gamma_par(Smax, Rmax, gamma)[[1]],
         b = gamma_par(Smax, Rmax, gamma)[[2]]) %>%
  rowwise() %>%
  mutate(lb = optimise(f = get_bounds, #'true' OYP bounds
                       interval = 1:get_Smsy(lnalpha_1, beta),
                       lnalpha = lnalpha_1,
                       beta = beta,
                       pct_MSY = 0.9,
                       correct = FALSE,
                       sigma = sigma,
                       phi = phi)$minimum) %>%
  mutate(lnalpha_2 = lb * beta,
         Rmax_2 = Ricker(lnalpha_2, beta, Smax),
         a_2 = gamma_par(Smax, Rmax_2, gamma)[[1]]) %>%
  select(-Smax, -Rmax, -Rmax_2, lb) %>%
  select(scenario, lnalpha_1, lnalpha_2, a_1, a_2, beta, b, gamma, sigma, phi, harvest_rate, sigN) %>%
  ungroup()

#  * Scenarios table ------------------------------------------------------
# for each scenario we have 50 replicate datasets and SR models
# scenarios %>% 
#   select(scenario, lnalpha_1, phi, gamma, harvest_rate) %>%
#   arrange(lnalpha_1, phi, gamma, harvest_rate) %>%
#   flextable() %>%
#   set_header_labels(
#     scenario = "Scenario",
#     lnalpha_1 = "ln(\u03b1)",
#     phi = "\u03d5",
#     gamma = "\u03b3",
#     harvest_rate = "regime 2 \n Seq / lb"
#   ) %>%
#   merge_v(j = c("lnalpha_1", "phi", "gamma", "harvest_rate")) %>%
#   valign(j = c("lnalpha_1", "phi", "gamma", "harvest_rate"), valign = "top") %>%
#   autofit()


#Simulate data and estimate parameters --------
get_SOC2 <- function(x, misses_in, makes_out, window = 5){
  SOC <- character()
  SOC[1:min(length(x), (window - 1))] <- NA
  from_noconcern <- function(x){
    if(sum(x, na.rm = TRUE) >= misses_in){"SOC"} else("No concern")
  }
  from_SOC <- function(x){
    if(sum(x, na.rm = TRUE) > (window - makes_out)){"SOC"} else("No concern")
  }
  if(length(x) >= window){
    for(i in window:length(x)){
      SOC[i] <- switch(as.character(SOC[i-1]),
                       'NA' = from_noconcern(x[((i-1)-(window - 1)):(i-1)]),
                       "No concern" = from_noconcern(x[((i-1)-(window - 1)):(i-1)]),
                       "SOC" = from_SOC(x[((i-1)-(window - 1)):(i-1)]))
    }
  }
  if(length(x) < window){
    for(i in length(x)){
      SOC[i] <- NA
    }
  }
  out <- ifelse(SOC == "SOC", TRUE, ifelse(SOC == "No concern", FALSE, NA))
  
  return(out)
}

rep_scenarios <-
  scenarios %>%
  slice(rep(1:nrow(scenarios), each = 50)) %>%
  mutate(rep = rep(1:50, times = nrow(scenarios))) %>%
  select(scenario, rep, lnalpha_1, lnalpha_2, a_1, a_2, beta, b, gamma, sigma, phi, harvest_rate, sigN) %>%
  rowwise() %>%
  mutate(data = list(
    sim_SRgamma(c(rep(a_1, 20), rep(a_2, 30)),
                b,
                gamma = gamma,
                sigW = sigma,
                phi = 0,
                age0 = c('3' = 0.1, '4' = 0.38, '5' = 0.3, '6' = 0.2, '7' = 0.02),
                Sims0 = 50,
                Hfun = "rate",
                harvest_rate = harvest_rate,
                sigF = 0.2,
                sigN = sigN))) %>%
  ungroup()

rep_scenarios %>% 
  mutate(contrast = map_dbl(data, function(x) max(x$S[1:20]) / min(x$S[1:20]))) %>%
  group_by(lnalpha_1, phi, gamma, harvest_rate, sigN) %>%
  summarise(n = n(),
            median = median(contrast),
            min = min(contrast),
            max = max(contrast)) %>%
  print(n = 100)

pct_goals <- function(harvest_rate, contrast, sigN){
  tier <- 
    ifelse(harvest_rate < 0.4, 
           ifelse(contrast > 8, 
                  ifelse(sigN <= 0.15, "Tier 1", "Tier 2"),
           "Tier 3"),
    NA)
  
  pct <- data.frame(
    tier = c("Tier 1", "Tier 2", "Tier 3", NA),
    lb = c(0.2, 0.15, 0.05, NA),
    ub = c(0.6, 0.65, 0.65, NA))
  
  pct[pct$tier %in% tier, ]
}

#test <- 
rep_scenarios %>% 
  mutate(contrast = map_dbl(data, function(x) max(x$S[1:20]) / min(x$S[1:20]))) %>%
  filter(!is.na(contrast)) %>%
  mutate(sim = map(data, function(x) x$sim),
         S = map(data, function(x) x$S)) %>%
  unnest(c(sim, S)) %>%
  group_by(scenario, rep) %>%
  mutate(
    pct_tier_lb = pct_goals(harvest_rate, contrast, sigN)[, 2],
    pct_tier_ub = pct_goals(harvest_rate, contrast, sigN)[, 3],
    tier_lb = mapply(function(i, j) quantile(S[1:i], probs = pct_tier_lb[j]), i = seq_along(S), j = row_number()),
    tier_ub = mapply(function(i, j) quantile(S[1:i], probs = pct_tier_ub[j]), i = seq_along(S), j = row_number())) %>% 
  select(scenario, rep, sim, tier_lb, tier_ub) %>% 
  print(n = 100)



  unnest(cols = c(sim, S)) %>%

  mutate()
  mutate(lb_goal = zoo::rollapply(S, width = 5, 
                                       FUN = pct_goals, 
                                       harvest_rate = unique(harvest_rate), 
                                       contrast = contrast, 
                                       sigN = unique(sigN),
                                       bound = "lb"))



test %>% 
  select(scenario, rep, S, lb_goal) %>% 
  filter(scenario == 1, rep == 1)


  print(n = 100)
  

%>%
  mutate(sim = map(data, function(x) x$sim),
         S = map(data, function(x) x$S),
         N = map(data, function(x) x$N)) %>%
  unnest(c(sim, S, N))
  group_by(scenario, rep) %>%
  mutate()
  mutate(yr_SOC = map_dbl(data, function(x) which.min(x$SOC != "Management")))

# ** SOC listings ---------------------------------------------------------
# First SOC finding
table(rep_scenarios$yr_SOC)
plot_SOCyr <- 
  rep_scenarios %>%
    ggplot(aes(x = yr_SOC, fill = as.character(lnalpha_1))) +
    geom_bar()
plot_SOCyr + 
  facet_wrap_paginate(~ paste0("%MSY @ lb: ", pct_MSY) +
                      paste0("%lb @ Seq: ", pct_lb) +
                      paste0("\u03B3: ", gamma), 
                      ncol = 3, nrow = 2, page = 1)
plot_SOCyr + 
  facet_wrap_paginate(~ paste0("%MSY @ lb: ", pct_MSY) +
                        paste0("%lb @ Seq: ", pct_lb) +
                        paste0("\u03B3: ", gamma), 
                      ncol = 3, nrow = 2, page = 2)  

# ** q50 S vrs. year. ---------  
plot_S <-
  rep_scenarios %>%
  mutate(sim = map(data, function(x) x$sim),
         S = map(data, function(x) x$S)) %>%
  unnest(c(sim, S)) %>%
  group_by(scenario, sim, lnalpha_1, lnalpha_2, a_1, a_2, b, gamma, sigma, pct_MSY, pct_lb, lb) %>%
  summarize(S = median(S)) %>%
  mutate(SET = (gamma - 1) / b) %>%
  ggplot(aes(x = sim, y = S, color = as.character(lnalpha_1))) +
  geom_line() +
  geom_hline(aes(yintercept = lb, color = as.character(lnalpha_1))) +
  geom_hline(aes(yintercept = SET), linetype = 2)  
plot_S + 
  facet_wrap_paginate(~ paste0("%MSY @ lb: ", pct_MSY) +
                        paste0("%lb @ Seq: ", pct_lb) +
                        paste0("\u03B3: ", gamma), 
                      ncol = 3, nrow = 2, page = 1)
plot_S + 
  facet_wrap_paginate(~ paste0("%MSY @ lb: ", pct_MSY) +
                        paste0("%lb @ Seq: ", pct_lb) +
                        paste0("\u03B3: ", gamma), 
                      ncol = 3, nrow = 2, page = 2)

# Analysis of sim data ----------------------------------------------------

#Model notes
#gamma_RS_change or Ricker _RS_change: SR regression on the log(R/S) scale with a change point 
#     to detect regime shifts. The change point is not very sophisticated but sufficient 
#     for simulate data. For real data the regime detection could be different.
#gamma_RS_change_set: Set prior on the set instead of log(alpha). Avoiding a prior on log(alpha[2]) 
#     bc there appears to be a near linear relationship between log(alpha) and gamma 
#     in the gamma SR relationship.
#gamma_RS_change_scale: Add a parameter scale (which is ln(Rmax)) so I can set a prior on it 
#     and not on alpha[2]. Avoiding a prior on log(alpha[2]) bc there is a 
#     linear relationship between log(alpha) and gamma in the gamma SR relationship.
#     log(alpha) = log(rmax) - gamma(log(Smax) - 1)
#gamma_RS_change_scale_hier: Tried a hierarcical prior on gamma[2]. No sure it helped much.

#Posterior notes
#In general the posterior is named after the model.
#rep_scenarios_seqlb: An early run based on a version of gamma_RS_change (changes are to the priors)   
#rep_scenarios_seqlb_scale_SOC: gamma_RS_change_scale but 50 year simulated  datasets that are 
#     truncated before analysis to the first year when an SOC would have been listed 
#     (4 out of 5 years with N < lb)
rep_scenarios_seqlb_scale_SOC <-
  rep_scenarios %>%
  filter(yr_SOC != 1, !is.na(yr_SOC)) %>%
  mutate(data_jags =
            map2(data, yr_SOC, ~.x[1:.y,] %>%
                 {list(nyrs = max(.$sim),
                       lnRS = log(.$R / .$S),
                       S = .$S,
                       ar1 = 0)}),
         mod_gamma = map(data_jags, ~ jags(data = .x,
                                           parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "scale"),
                                           model.file = ".\\scripts\\gamma_RS_change_scale.txt",
                                           n.chains = 3,
                                           n.iter = 5e4,
                                           n.burnin = 1e4,
                                           n.thin = 480, 
                                           parallel = TRUE))
  )
saveRDS(rep_scenarios_seqlb_scale_SOC, file = ".\\rep_scenarios_seqlb_scale_SOC.rds")
rep_scenarios_seqlb_scale_SOC <-
  rep_scenarios_seqlb_scale_SOC %>%
  mutate(mod_Ricker = map(data_jags, ~ jags(data = .x,
                                           parameters.to.save = c("lnalpha", "beta", "sigma", "y_d"),
                                           model.file = ".\\scripts\\Ricker_RS_change.txt",
                                           n.chains = 3,
                                           n.iter = 5e4,
                                           n.burnin = 1e4,
                                           n.thin = 480, 
                                           parallel = TRUE))
  )
saveRDS(rep_scenarios_seqlb_scale_SOC, file = ".\\rep_scenarios_seqlb_scale_SOC.rds")
rep_scenarios_seqlb_scale_SOC <- readRDS(file = ".\\rep_scenarios_seqlb_scale_SOC.rds")


# Results -----------------------------------------------------------------
# * Rhat boxplot ---------
# gamma model
# for main parameters (lnalpha, beta, gamma, sigma, y_d)
rep_scenarios_seqlb_scale_SOC %>%
  rowwise() %>%
  mutate(Rhat_gamma = list(mod_gamma$summary[c("lnalpha[1]", "lnalpha[2]",
                                               "beta[1]", "beta[2]",
                                               "gamma[1]", "gamma[2]",
                                               "sigma[1]", "sigma[2]",
                                               "y_d"), "Rhat"])) %>%
  select(scenario, lnalpha_1, beta, sigma, phi, pct_MSY, pct_lb, Rhat_gamma) %>%
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
rep_scenarios_seqlb_scale_SOC %>%
  rowwise() %>%
  mutate(Rhat_Ricker = list(mod_Ricker$summary[c("lnalpha[1]", "lnalpha[2]",
                                               "beta[1]", "beta[2]",
                                               "sigma[1]", "sigma[2]",
                                               "y_d"), "Rhat"])) %>%
  select(scenario, lnalpha_1, beta, sigma, phi, pct_MSY, pct_lb, Rhat_Ricker) %>%
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

#  * SET q50s --------------------------------------
# median SET for 50 replicates under each scenario
# black line is the SET from the data generating model
# Hard to see the effect of pct_lb
rep_scenarios_seqlb_scale_SOC %>%
  select(scenario, rep, lnalpha_1, lnalpha_2, beta, sigma, gamma, a_1, a_2, b, pct_lb, pct_MSY, lb, mod_gamma) %>%
  mutate(set = (gamma - 1) / b) %>%
  mutate(set_est = map_dbl(mod_gamma, function(x) (x$q50$gamma[2] - 1) / x$q50$beta[2])) %>%
  ggplot(aes(x = as.character(sigma), y = set_est, color = as.character(lnalpha_1), group = scenario)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = set)) +
  #geom_point(aes(y = lb), shape = 8, color = "red") +
  facet_grid(paste0("\u03B3: ", gamma) ~ 
               paste0("% of MSY: ", pct_MSY) +
               paste0("% of lb: ", pct_lb), scales = "free_x")


#  * SET simulated values -------------------------------------------------
# Simulated values for the SET in each scenario (all replicates)
# black line is the SET from the data generating model
# dashed line is the lower bound of the escapement goal (50% or 90% of MSY)
plot_SETsims <- 
  rep_scenarios_seqlb_scale_SOC %>%
  mutate(set = (gamma - 1) / b, 
         set_sims = map(mod_gamma, function(x) (x$sims.list$gamma[,2] - 1) / x$sims.list$beta[,2])) %>%
  select(scenario, rep, lnalpha_1, lnalpha_2, beta, sigma, gamma, a_1, a_2, b, pct_lb, pct_MSY, lb, set, set_sims) %>%
  unnest(set_sims) %>%
  ggplot(aes(x = set_sims, fill = as.character(rep), group = rep)) + 
  geom_histogram() +
  geom_vline(aes(xintercept = set)) +
  geom_vline(aes(xintercept = lb), linetype = 2) +
  scale_x_continuous(limits = c(0, 12000)) +
  theme(legend.position="none")
plot_SETsims +
  facet_grid_paginate(paste0("\u03B3: ", gamma) ~ paste0("ln(\u03B1): ", lnalpha_1) +
                        paste0("% of MSY: ", pct_MSY) + 
                        paste0("\u03B4: ", sigma) +
                        paste0("% of lb: ", pct_lb), nrow = 3, ncol = 4, page = 1)
plot_SETsims +
  facet_grid_paginate(paste0("\u03B3: ", gamma) ~ paste0("ln(\u03B1): ", lnalpha_1) +
                        paste0("% of MSY: ", pct_MSY) + 
                        paste0("\u03B4: ", sigma) +
                        paste0("% of lb: ", pct_lb), nrow = 3, ncol = 4, page = 2)
plot_SETsims +
  facet_grid_paginate(paste0("\u03B3: ", gamma) ~ paste0("ln(\u03B1): ", lnalpha_1) +
                        paste0("% of MSY: ", pct_MSY) + 
                        paste0("\u03B4: ", sigma) +
                        paste0("% of lb: ", pct_lb), nrow = 3, ncol = 4, page = 3)
plot_SETsims +
  facet_grid_paginate(paste0("\u03B3: ", gamma) ~ paste0("ln(\u03B1): ", lnalpha_1) +
                        paste0("% of MSY: ", pct_MSY) + 
                        paste0("\u03B4: ", sigma) +
                        paste0("% of lb: ", pct_lb), nrow = 3, ncol = 4, page = 4)
plot_SETsims +
  facet_grid_paginate(paste0("\u03B3: ", gamma) ~ paste0("ln(\u03B1): ", lnalpha_1) +
                        paste0("% of MSY: ", pct_MSY) + 
                        paste0("\u03B4: ", sigma) +
                        paste0("% of lb: ", pct_lb), nrow = 3, ncol = 4, page = 5)
plot_SETsims +
  facet_grid_paginate(paste0("\u03B3: ", gamma) ~ paste0("ln(\u03B1): ", lnalpha_1) +
                        paste0("% of MSY: ", pct_MSY) + 
                        paste0("\u03B4: ", sigma) +
                        paste0("% of lb: ", pct_lb), nrow = 3, ncol = 4, page = 6)

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
    select(scenario, rep, a_1, a_2, b, gamma, lb, ub) %>%
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
    mutate(lnRS = map(data_jags, function(x) data.frame(year = 1:x$nyrs, lnRS = x$lnRS, S = x$S)),
           cut = map_dbl(mod_gamma, ~.$q50$y_d)) %>%
    select(scenario, rep, a_1, a_2, b, gamma, lb, ub, cut, lnRS) %>%
    unnest(lnRS) %>%
    mutate(regime = ifelse(year <= cut, "normal", "reduced")) %>%
    ggplot(aes(x = S, y = lnRS)) +
    geom_point(aes(color = regime)) +
    geom_vline(data = true, aes(xintercept = set), color = "red") +
    geom_vline(data = estimates, aes(xintercept = set), linetype = 2, color = "red") +
    geom_rect(data = true, aes(xmin = lb, xmax = ub, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, 
              alpha= .25) +
    geom_function(fun = function(x) log(SRgamma(true$a_1, true$b, true$gamma, x) / x) ,
                  color = "black") +
    geom_function(fun = function(x) log(SRgamma(true$a_2, true$b, true$gamma, x) / x),
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
    mutate(S = map(data_jags, function(x) data.frame(year = 1:x$nyrs, S = x$S)),
           cut = map_dbl(mod_gamma, ~.$q50$y_d)) %>%
    unnest(S) %>%
    mutate(regime = ifelse(year <= cut, "normal", "reduced")) %>%
    select(scenario, rep, year, b, gamma, lb, regime, S) %>% 
    mutate(set = (gamma - 1)/ b) %>% 
    ggplot(aes(x = year, y = S)) + 
    geom_line() +
    geom_point(aes(color = regime)) +
    geom_hline(data = true, aes(yintercept = set), color = "red") +
    geom_hline(data = estimates, aes(yintercept = set), linetype = 2, color = "red") +
    geom_rect(data = true, aes(ymin = lb, ymax = ub, xmin = -Inf, xmax = Inf),
              inherit.aes = FALSE, 
              alpha= .25) +
    scale_color_manual(values = c("black", "red")) +
    scale_y_continuous(limits = c(0, 1.1 * maxS)) +
    theme_bw() +
    labs(caption = paste0(data$scenario, ", ", data$rep))
  
  ggarrange(plot1, plot2, nrow = 2, heights = c(6,4))
}
plot_fit(rep_scenarios_seqlb_scale_SOC, 
         lnalpha_1 == 1.5,
         sigma == .5,
         gamma == 1,
         pct_MSY == 0.5,
         pct_lb == 1)
#Use this code to check SOC designations (I have checked several and they look good).
rep_scenarios_seqlb_scale_SOC %>% 
  filter(scenario == 38, rep == 12) %>% 
  select(data) %>% 
  unnest(data) %>% 
  select(sim, S, N, lb_goal, SOC) %>% 
  tail(n = 20)

# * Horsetail ------------------------------------
plot_horse <- function(scenario_file, ...){
  filters <- enquos(...)
  
  data0 <- 
    scenario_file %>%
    filter(!!!filters)
  
  rep_var = sample(unique(data0$rep), 1, replace = FALSE)
  
  data <- 
    data0 %>%
    filter(rep == rep_var)
  
  coeflines1 <-
    data.frame(lnalpha = exp(data$mod_gamma[[1]]$sims.list$lnalpha[, 1]),
               beta = data$mod_gamma[[1]]$sims.list$beta[, 1],
               gamma = data$mod_gamma[[1]]$sims.list$gamma[, 1]) %>%
    dplyr::sample_n(50) %>%
    as.matrix() %>%
    plyr::alply(1, function(coef) {
      ggplot2::stat_function(fun = function(x) log(coef[1] * x^(coef[3] - 1) * exp(-coef[2]*x)), 
                             colour="grey",
                             alpha = 0.5)})
  coeflines2 <-
    data.frame(lnalpha = exp(data$mod_gamma[[1]]$sims.list$lnalpha[, 2]),
               beta = data$mod_gamma[[1]]$sims.list$beta[, 2],
               gamma = data$mod_gamma[[1]]$sims.list$gamma[, 2]) %>%
    dplyr::sample_n(50) %>%
    as.matrix() %>%
    plyr::alply(1, function(coef) {
      ggplot2::stat_function(fun = function(x) log(coef[1] * x^(coef[3] - 1) * exp(-coef[2]*x)), 
                             colour="grey",
                             alpha = 0.5)})
  
  data.frame(S = data$data_jags[[1]]$S,
             lnRS = data$data_jags[[1]]$lnRS) %>%
    ggplot(aes(x = S, y = lnRS)) +
    geom_point() +
    geom_rect(data = data,
              aes(xmin = lb, 
                  xmax = ub, 
                  ymin = -Inf, 
                  ymax = Inf),
              inherit.aes = FALSE,
              fill = "grey",
              alpha= .25) +
    scale_x_continuous(limits = c(0, 1.05 * max(data$data_jags[[1]]$S))) +
    theme_bw() +
    coeflines1 + 
    coeflines2
}
plot_horse(rep_scenarios_seqlb_scale_SOC,
           lnalpha_1 == 1.5,
           sigma == .5,
           gamma == 1,
           pct_MSY == 0.5,
           pct_lb == 1)


#  * DIC compare ----------------------------------------------------------
# How often did the gamma model have the lower dic?
rep_scenarios_seqlb_scale_SOC %>%
  rowwise() %>%
  mutate(dic_score = mod_gamma$DIC < mod_Ricker$DIC) %>%
  select(scenario, lnalpha_1, beta, sigma, phi, gamma, pct_MSY, pct_lb, starts_with("dic")) %>%
  group_by(scenario, lnalpha_1, sigma, gamma, pct_MSY, pct_lb) %>%
  summarise(dic_mean = mean(dic_score)) %>%
  ggplot(aes(x = as.character(pct_lb), 
             y = dic_mean, 
             color = as.character(lnalpha_1),
             shape = as.character(sigma))) + 
  geom_point() +
  geom_hline(yintercept = 0.5) +
  facet_grid(paste0("\u03B3: ", gamma) ~ paste0("ln(\u03B1): ", lnalpha_1) +
               paste0("% of MSY: ", pct_MSY)) # +
               #paste0("% of lb: ", pct_lb))



#  * gamma simulated values ------------------------------------------------------------
#gamma mostly spans prior. 
# don't really believe some of the prior (1.6-2)
# If gamma could take those values we should know it.
plot_gammasims <- 
  rep_scenarios_seqlb_scale_SOC %>%
  mutate(gamma_sims = map(mod_gamma, function(x) x$sims.list$gamma[,2])) %>%
  select(scenario, rep, lnalpha_1, lnalpha_2, beta, sigma, gamma, a_1, a_2, b, pct_lb, pct_MSY, lb, gamma_sims) %>%
  unnest(gamma_sims) %>%
  ggplot(aes(x = gamma_sims, fill = as.character(rep), group = rep)) + 
  geom_histogram() +
  geom_vline(aes(xintercept = gamma)) +
  theme(legend.position="none")
plot_gammasims +
  facet_grid_paginate(paste0("\u03B3: ", gamma) ~ paste0("ln(\u03B1): ", lnalpha_1) +
                        paste0("% of MSY: ", pct_MSY) + 
                        paste0("\u03B4: ", sigma) +
                        paste0("% of lb: ", pct_lb), nrow = 3, ncol = 4, page = 1)
plot_gammasims +
  facet_grid_paginate(paste0("\u03B3: ", gamma) ~ paste0("ln(\u03B1): ", lnalpha_1) +
                        paste0("% of MSY: ", pct_MSY) + 
                        paste0("\u03B4: ", sigma) +
                        paste0("% of lb: ", pct_lb), nrow = 3, ncol = 4, page = 2)
plot_gammasims +
  facet_grid_paginate(paste0("\u03B3: ", gamma) ~ paste0("ln(\u03B1): ", lnalpha_1) +
                        paste0("% of MSY: ", pct_MSY) + 
                        paste0("\u03B4: ", sigma) +
                        paste0("% of lb: ", pct_lb), nrow = 3, ncol = 4, page = 3)
plot_gammasims +
  facet_grid_paginate(paste0("\u03B3: ", gamma) ~ paste0("ln(\u03B1): ", lnalpha_1) +
                        paste0("% of MSY: ", pct_MSY) + 
                        paste0("\u03B4: ", sigma) +
                        paste0("% of lb: ", pct_lb), nrow = 3, ncol = 4, page = 4)
plot_gammasims +
  facet_grid_paginate(paste0("\u03B3: ", gamma) ~ paste0("ln(\u03B1): ", lnalpha_1) +
                        paste0("% of MSY: ", pct_MSY) + 
                        paste0("\u03B4: ", sigma) +
                        paste0("% of lb: ", pct_lb), nrow = 3, ncol = 4, page = 5)
plot_gammasims +
  facet_grid_paginate(paste0("\u03B3: ", gamma) ~ paste0("ln(\u03B1): ", lnalpha_1) +
                        paste0("% of MSY: ", pct_MSY) + 
                        paste0("\u03B4: ", sigma) +
                        paste0("% of lb: ", pct_lb), nrow = 3, ncol = 4, page = 6)


#  * gamma q50s -----------------------------------------------------------
# point estimates are not terrible
# I think this was helped greatly by "scale" prior parameterzation
rep_scenarios_seqlb_scale_SOC %>%
  select(scenario, rep, lnalpha_1, lnalpha_2, beta, sigma, gamma, a_1, a_2, b, pct_lb, pct_MSY, lb, mod_gamma) %>%
  mutate(set = (gamma - 1) / b) %>%
  mutate(gamma_est = map_dbl(mod_gamma, function(x) x$q50$gamma[2])) %>%
  ggplot(aes(x = as.character(sigma), y = gamma_est, color = as.character(lnalpha_1), group = scenario)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = gamma)) +
  facet_grid(paste0("\u03B3: ", gamma) ~ 
               paste0("% of MSY: ", pct_MSY) +
               paste0("% of lb: ", pct_lb), scales = "free_x")

#notice that variation in beta counteracted variation in gamma to stabilize the SET.
# y axis is (gamma - 1)/gamma of the SET as a percentage of Smax 
# x axis in Smax
# sometimes!!!
plot_SETcomponents <- function(scenario_file, ...){
  filters <- enquos(...)
  
  data0 <- 
    scenario_file %>%
    filter(!!!filters)
  
  rep_var = sample(unique(data0$rep), 9, replace = FALSE)
  
  data0 %>%
    filter(rep %in% rep_var) %>%
    mutate(gamma_est = map_dbl(mod_gamma, function(x) x$q50$gamma[2]),
           beta_est = map_dbl(mod_gamma, function(x) x$q50$beta[2]),
           gamma_sims = map(mod_gamma, function(x) x$sims.list$gamma[, 2]),
           beta_sims = map(mod_gamma, function(x) x$sims.list$beta[, 2])) %>%
    unnest(c(gamma_sims, beta_sims)) %>%
    group_by(rep) %>%
    mutate(gamma_sims_1 = gamma_sims - 1,
           set = (gamma - 1) / b,
           set_est = (gamma_est - 1) / beta_est,
           set_sims = gamma_sims_1 / beta_sims,
           set_est_cut = cut(set_sims, 
                             quantile(set_sims, c(0, 0.5, 1)),
                             labels = c("lower 50%", "upper 50%"))) %>%
    ungroup() %>%
    ggplot(aes(x = beta_sims, y = gamma_sims_1)) +
    geom_point() +
    geom_smooth(method = "lm") +
    geom_point(aes(color = set_est_cut)) +
    geom_abline(aes(intercept = 0, slope = set), linewidth = 1) +
    geom_abline(aes(intercept = 0, slope = set_est), color = "red", linewidth = 1) +
    facet_wrap(.~rep)
}
plot_SETcomponents(
  rep_scenarios_seqlb_scale_SOC, 
  lnalpha_1 == 1,
  sigma == .5,
  gamma == 1.3,
  pct_MSY == 0.9,
  pct_lb == 1)
  
