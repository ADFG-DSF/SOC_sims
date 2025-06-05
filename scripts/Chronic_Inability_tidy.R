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
              sigma = c(0.5, 0.8),
              phi = c(0, 0.4, 0.8),
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
# NEEDS TO BE UPDATED AFTER SCENARIOS FINALIZED #####
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
    sim_SRgamma(c(rep(a_1, 100), rep(a_2, 130)),
                b,
                gamma = gamma,
                sigW = sigma,
                phi = 0,
                age0 = c('3' = 0.1, '4' = 0.38, '5' = 0.3, '6' = 0.2, '7' = 0.02),
                Sims0 = 230,
                Hfun = H_goal,
                lb_goal = lb,
                ub_goal = ub,
                power = power,
                sigF = 0.2,
                sigN = 0.2))) %>%
  ungroup()
#saveRDS(rep_scenarios, ".//rep_scenarios_ChronicInability.rds")
rep_scenarios <- readRDS(".//rep_scenarios_ChronicInability.rds")

# SOC listings ---------------------------------------------------------
# function to calculate SOC status
# window: # of years under consideration when evaluating "chronic inability"
# misses_in: # of missed in window that result in a SOC listing
# makes_out: # of makes in window that result in a SOC delisting.
# originally written to work inside a list-column although it seems to work in groups too.
get_SOC2 <- function(x, misses_in, makes_out, window = 5){
  SOC <- character()
  SOC[1:(window - 1)] <- NA
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
  out <- ifelse(SOC == "SOC", TRUE, ifelse(SOC == "No concern", FALSE, NA))
  
  return(out)
}

# SOC status for different definitions of chronic inability
# N### naming convention is {misses_in}{makes_out}{window}
dat_SOC <- 
  rep_scenarios %>%
  filter(scenario != 64, rep != 17) %>%
  filter(pct_lb == 1) %>%
  mutate(data_trunc = map(data, function(x) x[x$sim > 30, c("sim", "N", "lb_goal")])) %>%
  unnest(data_trunc) %>%
  select(scenario:pct_lb, power, sim:lb_goal) %>%
  filter(sim < 100 | sim >= 130) %>%
  mutate(regime = ifelse(sim < 100, "Historic", ifelse(sim >= 130, "Low", NA))) %>%
  group_by(scenario, regime, rep) %>%
  mutate(
    roll5 = rollmean(N, k = 5, align = "right", fill = NA),
    SOC_roll5 = roll5 < lb_goal, #rolling mean to account for the size of the miss/make
    mean_roll5 = mean(SOC_roll5, na.rm = TRUE),
    lb_N = N < lb_goal,
    SOC_N155 = get_SOC2(lb_N, 1, 5, 5), #most sensitive
    mean_155 = mean(SOC_N155, na.rm = TRUE),
    SOC_N445 = get_SOC2(lb_N, 4, 4, 5), #ADF&G most common?
    mean_445 = mean(SOC_N445, na.rm = TRUE),
    SOC_N515 = get_SOC2(lb_N, 5, 1, 5), #lease sensitive
    mean_515 = mean(SOC_N515, na.rm = TRUE)) %>%
  group_by(scenario, rep) 

# demonstrate regime and SOC entry/exit criteria differences for one group of scenarios
# Hypothesis test for 445
# NULL: High productivity
# false positives (Look at Historic regime): non-existent, i.e. alpha = 0
# false negatives  (Look at Low regime): 
#     gamma=1: Fail to reject Null often, beta ~ .35,  
#     gamma>1: Reject null occasionally in "safest" situation, beta > .75
# Entry/Exit criteria: would need a very sensitive test to get high power w Ricker 
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, phi == 0) %>%
  pivot_longer(dplyr::starts_with("mean_"), names_to = "Chronic_I", values_to = "pct_SOC") %>%
  ggplot(aes(x = Chronic_I, y  = pct_SOC, color = regime)) + 
  geom_boxplot() +
  guides(x =  guide_axis(angle = -20)) +
  scale_color_discrete(name = "Regime (statistic)", 
                     labels = c("Historic" = "Historic (\u03B1)", "Low" = "Low (\u03B2)")) +
  labs(x = "Chronic Inability Criteria",
       y = "Probability") +
  facet_grid(paste0("%MSY: ", pct_MSY) ~ paste0("\u03B3: ", gamma))

####
# Look for sensitivity wrt 445 and roll5 across SR parameter combinations
# pct_MSY, phi and gamma
# phi has little effect
dat_SOC %>%
  filter(regime == "Low", sigma == 0.5, lnalpha_1 == 1.5) %>%
  pivot_longer(dplyr::starts_with("mean_"), names_to = "Chronic_I", values_to = "pct_SOC") %>%
  filter(Chronic_I %in% c("mean_445", "mean_roll5")) %>%
  ggplot(aes(x = Chronic_I, y  = pct_SOC, color = as.character(phi))) + 
  geom_boxplot() +
  labs(title = paste0("ln(\u03B1): ", 1.5, ", \u03C3: ", 0.5)) +
  facet_grid(paste0("%MSY: ", pct_MSY) ~ paste0("\u03B3: ", gamma))

# Look for sensitivity wrt 445 and roll5
# pct_MSY, lnalpha and gamma
# increasing lnalpha decreases sensitivity w depensatory dynamics
dat_SOC %>%
  filter(regime == "Low", phi == 0, sigma == 0.5) %>%
  pivot_longer(dplyr::starts_with("mean_"), names_to = "Chronic_I", values_to = "pct_SOC") %>%
  filter(Chronic_I %in% c("mean_445", "mean_roll5")) %>%
  ggplot(aes(x = Chronic_I, y  = pct_SOC, color = as.character(lnalpha_1))) + 
  geom_boxplot() +
  labs(title = paste0("%MSY: ", 0.9, ", \u03C3: ", 0.5)) +
  facet_grid(paste0("%MSY: ", pct_MSY) ~ paste0("\u03B3: ", gamma))

# Look for sensitivity wrt 445 and roll5
# pct_MSY, sigma and gamma
# increased sigma decreases sensitivity
dat_SOC %>%
  filter(regime == "Low", phi == 0, lnalpha_1 == 1.5) %>%
  pivot_longer(dplyr::starts_with("mean_"), names_to = "Chronic_I", values_to = "pct_SOC") %>%
  filter(Chronic_I %in% c("mean_445", "mean_roll5")) %>%
  ggplot(aes(x = Chronic_I, y  = pct_SOC, color = as.character(sigma))) + 
  geom_boxplot() +
  labs(title = paste0("phi: ", 0, ", ln(\u03B1): ", 1.5)) +
  facet_grid(paste0("%MSY: ", pct_MSY) ~ paste0("\u03B3: ", gamma))








# Look for differences between SOC 445 and SOC_roll5
# Have to manually enter the sim numbers in the last line
# 1) Run lines 1-4 to find sim numbers were criteria differ
# 2) comment out line 4 and either 5 or 6 to see why the criteria differ

# when N45 more sensitive....
# rolling mean provides the protection DCF was looking for. SOC designations do not occur when
# several small misses are surrounded by 1 large make

# when roll5 more sensitive....
# it is pulled down by large misses surrounded by several close makes
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, phi == 0, pct_MSY == 0.9, gamma == 1, rep == 1) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, sim, N, lb_goal, SOC_N445, roll5, SOC_roll5) %>%
  filter(SOC_N445 != SOC_roll5) %>% print(n = 100)
 
col <- c("#00CC00", "red")
#N45 more sensitive
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, phi == 0, pct_MSY == 0.9, gamma == 1, rep == 1) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, sim, N, lb_goal, SOC_N445, roll5, SOC_roll5) %>%
  filter(sim >= (178-5) & sim <= 190) %>%
  ggplot(aes(x = as.character(sim), y = N, fill = SOC_N445)) +
  geom_bar(stat = "identity", alpha = 0.5) +
  geom_hline(aes(yintercept = lb_goal)) +
  geom_line(aes(y = roll5, group = 1, color = as.numeric(SOC_roll5)), 
            linewidth = 3, 
            show.legend = FALSE) +
  scale_color_gradient(low = col[1], high = col[2]) +
  scale_fill_manual(values = col) +
  labs(x = "Simulation Year", 
       y = "Total Run", 
       fill = "Stock of Concern?", 
       title = paste0("ln(\u03B1)= ", 1.5, ", \u03C3= ", 0.5, ", \u03B3= ", 1, ", \u03A6= ", 0, ", %MSY= ", 90))

#roll5 more sensitive 
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, phi == 0, pct_MSY == 0.9, gamma == 1, rep == 1) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, sim, N, lb_goal, SOC_N445, roll5, SOC_roll5) %>%
  filter(sim >= (154-5) & sim <= 159) %>%
  ggplot(aes(x = as.character(sim), y = N, fill = SOC_N445)) +
  geom_bar(stat = "identity", alpha = 0.5) +
  geom_hline(aes(yintercept = lb_goal)) +
  geom_line(aes(y = roll5, group = 1, color = as.numeric(SOC_roll5)), 
            linewidth = 3, 
            show.legend = FALSE) +
  scale_color_gradient(low = col[1], high = col[2]) +
  scale_fill_manual(values = col) +
  labs(x = "Simulation Year", 
       y = "Total Run", 
       fill = "Stock of Concern?", 
       title = paste0("ln(\u03B1)= ", 1.5, ", \u03C3= ", 0.5, ", \u03B3= ", 1, ", \u03A6= ", 0, ", %MSY= ", 90))

#roll5 more sensitive
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.8, phi == 0, pct_MSY == 0.9, gamma == 1.3, rep == 2) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, sim, N, lb_goal, SOC_N445, roll5, SOC_roll5) %>%
  filter(sim >= (134-5) & sim <= 150) %>%
  ggplot(aes(x = as.character(sim), y = N, fill = SOC_N445)) +
  geom_bar(stat = "identity", alpha = 0.5) +
  geom_hline(aes(yintercept = lb_goal)) +
  geom_line(aes(y = roll5, group = 1, color = as.numeric(SOC_roll5)), 
            linewidth = 3, 
            show.legend = FALSE) +
  scale_color_gradient(low = col[1], high = col[2]) +
  scale_fill_manual(values = col) +
  labs(x = "Simulation Year", 
       y = "Total Run", 
       fill = "Stock of Concern?", 
       title = paste0("ln(\u03B1)= ", 1.5, ", \u03C3= ", 0.8, ", \u03B3= ", 1.3, ", \u03A6= ", 0, ", %MSY= ", 90))
