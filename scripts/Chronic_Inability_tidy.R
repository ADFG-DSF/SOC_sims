# Chronic Inability

# Author: Adam Reimer
# Version: 2025-05-12

# Packages
packs <- c("tidyverse", "ggforce", "RcppRoll", "flextable", "zoo")
lapply(packs, require, character.only = TRUE)

# source functions
function_files <- list.files(path=".\\functions")
lapply(function_files, function(x) source(paste0(".\\functions\\", x)))

#Created under the assumption SOC recommendations in SOC.qmd are adopted.
# Work began after the MAY in person EGPIT meeting as we had a staff consensus and some concern 
# over limited progress on the chronic inability definition. 
# --------------------------------------------------------------------------------------
# Simulations -------
# * Define scenarios to simulate --------
scenarios <-
  expand.grid(lnalpha_1 = c(1, 1.5, 2),
              beta = 0.0001, #c(0.001, 0.0001, 0.000001),
              sigma = c(0.25, 0.5),
              phi = 0,
              pct_MSY = c(0.5, 0.9),
              pct_lb = c(0.5, 1), 
              gamma = seq(1, 1.6, length.out = 3)) %>%
  arrange(lnalpha_1, sigma, gamma, pct_MSY) %>%
  mutate(scenario = 1:n(),
         power = lnalpha_1 / 2 - .2,
         Smax = 1/ beta,
         Rmax = Ricker(lnalpha_1, beta, Smax),
         a_1 = gamma_par(Smax, Rmax, gamma)[[1]],
         b = gamma_par(Smax, Rmax, gamma)[[2]]) %>% # very crude. 
  rowwise() %>%
  mutate(lb = optimise(f = get_bounds, #'true' OYP bounds
                       interval = 1:get_Smsy(lnalpha_1, beta),
                       lnalpha = lnalpha_1,
                       beta = beta,
                       pct_MSY = pct_MSY,
                       correct = FALSE,
                       sigma = sigma,
                       phi = phi)$minimum,
         ub = optimise(f = get_bounds, #'true' OYP bounds
                       interval = get_Smsy(lnalpha_1, beta):(get_Smsy(lnalpha_1, beta)*5),
                       lnalpha = lnalpha_1,
                       beta = beta,
                       pct_MSY = 0.7,
                       correct = FALSE,
                       sigma = sigma,
                       phi = phi)$minimum) %>%
  mutate(lnalpha_2 = lb * pct_lb * beta,
         Rmax_2 = Ricker(lnalpha_2, beta, Smax),
         a_2 = gamma_par(Smax, Rmax_2, gamma)[[1]]) %>%
  select(-Smax, -Rmax, -Rmax_2) %>%
  select(scenario, lnalpha_1, lnalpha_2, a_1, a_2, beta, b, gamma, sigma, phi, pct_MSY, pct_lb, lb, ub, power) %>%
  ungroup()

#  * Scenarios table ------------------------------------------------------
# for each scenario we have 50 replicate datasets and SR models
scenarios %>% 
  select(scenario, lnalpha_1, sigma, gamma, pct_MSY, pct_lb) %>%
  arrange(lnalpha_1, sigma, gamma, pct_MSY, pct_lb) %>%
  flextable() %>%
  set_header_labels(
    scenario = "Scenario",
    lnalpha_1 = "ln(\u03b1)",
    sigma = "\u03c3",
    gamma = "\u03b3",
    pct_MSY = "MSY % @ EG \n lower bound",
    pct_lb = "regime 2 \n Seq / lb"
  ) %>%
  merge_v(j = c("lnalpha_1", "sigma", "gamma", "pct_MSY")) %>%
  valign(j = c("lnalpha_1", "sigma", "gamma", "pct_MSY", "pct_lb"), valign = "top") %>%
  autofit()


#Simulate data and estimate parameters --------
rep_scenarios <-
  scenarios %>%
  slice(rep(1:nrow(scenarios), each = 100)) %>%
  mutate(rep = rep(1:100, times = nrow(scenarios))) %>%
  select(scenario, rep, lnalpha_1, lnalpha_2, a_1, a_2, beta, b, gamma, sigma, phi, pct_MSY, pct_lb, lb, ub, power) %>%
  rowwise() %>%
  mutate(data = list(
    sim_SRgamma(c(rep(a_1, 1), rep(a_2, 234)),
                b,
                gamma = gamma,
                sigW = sigma,
                phi = 0,
                age0 = c('3' = 0.1, '4' = 0.38, '5' = 0.3, '6' = 0.2, '7' = 0.02),
                Sims0 = 235,
                Hfun = H_goal,
                lb_goal = lb,
                ub_goal = ub,
                power = power,
                sigF = 0.2,
                sigN = 0.2))) %>%
  ungroup()

# ** SOC listings ---------------------------------------------------------
get_SOC2 <- function(x, years){
  SOC <- character()
  SOC[1:4] <- NA
  from_noconcern <- function(x){
    if(sum(x, na.rm = TRUE) >= years){"SOC"} else("No concern")
  }
  from_SOC <- function(x){
    if(sum(x, na.rm = TRUE) > 5 - years){"SOC"} else("No concern")
  }
  if(length(x) >= 5){
    for(i in 5:length(x)){
      SOC[i] <- switch(as.character(SOC[i-1]),
                       'NA' = from_noconcern(x[(i-4):i]),
                       "No concern" = from_noconcern(x[(i-4):i]),
                       "SOC" = from_SOC(x[(i-4):i]))
    }
  }
  out <- ifelse(SOC == "SOC", TRUE, ifelse(SOC == "No concern", FALSE, NA))
  
  return(out)
}

dat_SOC <- 
  rep_scenarios %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, pct_lb == 1) %>%
  mutate(
    data_trunc = map(data, function(x) x[x$sim > 30, c("sim", "N", "lb_goal")]),
    roll5 = map(data_trunc, function(x) rollmean(x$N, k = 5, align = "right", fill = NA)),
    SOC_roll5 = map2(roll5, data_trunc, function(x, y) x < y$lb_goal),
    mean_roll5 = map_dbl(SOC_roll5, function(x) mean(x, na.rm = TRUE)),
    lb_N = map(data_trunc, function(x) x$N < x$lb_goal),
    SOC_N45 = map(lb_N, function(x) get_SOC2(x, 4)),
    mean_45 = map_dbl(SOC_N45, function(x) mean(x, na.rm = TRUE)),
    SOC_N55 = map(lb_N, function(x) get_SOC2(x, 5)),
    mean_55 = map_dbl(SOC_N55, function(x) mean(x, na.rm = TRUE))) %>%
  group_by(scenario, rep) 

# Look for differences
# Have to manually enter the sim numbers in the last line
# 1) Run lines 1-5 to find sim numbers were criteria differ 
# 2) comment out line 5 to see why the criteria differ

# when N45 more sensitive....
# rolling mean provides the protection DCF was looking for. SOC designations do not occur when 
# several small misses are surrounded by 1 large make

# when roll5 more sensitive....
# it is pulled down by large misses surrounded by several close makes
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, pct_MSY == 0.9, gamma == 1, rep == 1) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, data_trunc, SOC_N45, roll5, SOC_roll5) %>%
  unnest(c(data_trunc, SOC_N45, roll5, SOC_roll5)) %>%
  #filter(SOC_N45 != SOC_roll5) #%>%
  #filter(sim >= (71-6) & sim <= 76) #%>% #N45 more sensitive
  filter(sim >= (43-6) & sim <= 49) #roll5 more sensitive 
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, pct_MSY == 0.9, gamma == 1, rep == 5) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, data_trunc, SOC_N45, roll5, SOC_roll5) %>%
  unnest(c(data_trunc, SOC_N45, roll5, SOC_roll5)) %>%
  #filter(SOC_N45 != SOC_roll5) #%>%
  #filter(sim >= (115-6) & sim <= 125) #%>% #N45 more sensitive
  filter(sim >= (222-6) & sim <= 227) #roll5 more sensitive

# that said differences are minor
plot_SOC <-
  dat_SOC %>%
  pivot_longer(dplyr::starts_with("mean_"), names_to = "Chronic_I", values_to = "pct_SOC") %>%
  ggplot(aes(x = Chronic_I, y  = pct_SOC)) + 
  geom_boxplot() +
  guides(x =  guide_axis(angle = -20))

plot_SOC + 
  facet_wrap_paginate(~ paste0("%MSY @ lb: ", pct_MSY) +
                        #paste0("%lb @ Seq: ", pct_lb) +
                        paste0("\u03B3: ", gamma), 
                      ncol = 3, nrow = 2, page = 1)


rep_scenarios %>%
  select(-beta, -gamma, -phi) %>%
  unnest(data) %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, pct_lb == 1) %>%
  group_by(scenario, pct_MSY, gamma) %>%
  mutate(SOC_l = ifelse(SOC %in% c("Management", "Conservation"), TRUE, FALSE)) %>%
  summarize(SOC_pct = mean(SOC_l))