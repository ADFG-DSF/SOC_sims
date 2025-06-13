# Chronic Inability

# Author: Adam Reimer
# Version: 2025-05-12

# Packages
packs <- c("tidyverse", "ggforce", "RcppRoll", "flextable", "zoo")
lapply(packs, require, character.only = TRUE)

# source functions
function_files <- list.files(path=".\\functions")
lapply(function_files, function(x) source(paste0(".\\functions\\", x)))

col <- c("#00CC00", "red")

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
  arrange(lnalpha_1, sigma, phi, gamma, pct_MSY) %>%
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
# 216 scenarios
# for each scenario we have 100 replicate datasets and SR models
scenarios %>% 
  select(scenario, lnalpha_1, sigma, phi, gamma, pct_MSY, pct_lb) %>%
  arrange(lnalpha_1, sigma, phi, gamma, pct_MSY, pct_lb) %>%
  flextable() %>%
  set_header_labels(
    scenario = "Scenario",
    lnalpha_1 = "ln(\u03b1)",
    sigma = "\u03c3",
    phi = "\u03A6",
    gamma = "\u03b3",
    pct_MSY = "MSY % @ EG \n lower bound",
    pct_lb = "regime 2 \n Seq / lb"
  ) %>%
  merge_v(j = c("lnalpha_1", "sigma", "phi", "gamma", "pct_MSY")) %>%
  valign(j = c("lnalpha_1", "sigma", "phi", "gamma", "pct_MSY", "pct_lb"), valign = "top") %>%
  autofit()


#Simulate data and estimate parameters --------
# 216 scenarios 
# 100 reps
# 230 years per scenario: 100 historic regime, 130 low regime (first 30 reps of Low regime discarded)
rep_scenarios0 <-
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

# Because of the assumptions he have made some of our decompensatory scenarios go to extinction.
# This shows how many fail to make it at least 40 reps into the low productivity regime
# recall we discard the first 30 reps of the low productivity regime for burn in
# All gamma = 1.6, more common when the lower bound is targeting 50% of Smsy.
# Censor these. 
# Note if we change the criteria to 150 we would censor a lot more.
rep_scenarios0 %>% 
  mutate(sim_y = map_dbl(data, function(x) max(x$sim))) %>% 
  select(scenario, rep, lnalpha_1, sigma, phi, gamma, pct_MSY, pct_lb, sim_y) %>% 
  filter(sim_y < 140) %>%
  group_by(scenario, lnalpha_1, sigma, phi, gamma, pct_MSY, pct_lb) %>%
  summarise(n = length(scenario)) %>% 
  print(n = 150)

rep_scenarios <-
  rep_scenarios0 %>% 
  mutate(sim_y = map_dbl(data, function(x) max(x$sim))) %>% 
  filter(sim_y >= 140) %>%
  select(-sim_y)
#saveRDS(rep_scenarios, ".//rep_scenarios_ChronicInability.rds")
rep_scenarios <- readRDS(".//rep_scenarios_ChronicInability.rds")

# example time series of run sizes
plot_ts <- 
  rep_scenarios %>%
  filter(pct_lb == 1) %>%
  mutate(data_trunc = map(data, function(x) x[x$sim > 0, c("sim", "N", "lb_goal", "ub_goal")])) %>%
  unnest(data_trunc) %>%
  select(scenario:pct_lb, power, sim:ub_goal) %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, phi == 0)
plot_ts_median <- 
  plot_ts %>%
  group_by(scenario, gamma, pct_MSY, sim) %>%
  summarise(median_N = median(N))
test <- 
  plot_ts %>% 
  group_by(gamma, pct_MSY) %>% 
  summarise(lb_goal = unique(lb_goal), 
            ub_goal = unique(ub_goal))
plot_ts %>%
  ggplot(aes(x = sim, y = N)) +
  geom_line(aes(group = rep), color = "grey") +
  annotate("rect", xmin = 101, xmax = 107, ymin = 0, ymax = Inf, alpha = .5, fill = "red") +
  geom_rect(data = test, aes(xmin = 0, xmax = Inf, ymin = lb_goal, ymax = ub_goal), 
            inherit.aes = FALSE, 
            alpha = .2,) +
  geom_line(data = plot_ts_median, aes(y = median_N)) +
  coord_cartesian(ylim = c(0, 35000)) +
  facet_grid(paste0("%MSY: ", pct_MSY) ~ paste0("\u03B3: ", gamma)) +
  theme_bw() +
  labs(x = "Simulation Year", 
       y = "Total Run", 
       title = paste0("ln(\u03B1)= ", 1.5, ", \u03C3= ", 0.5, ", \u03B3= ", 1, ", \u03A6= ", 0))

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
  filter(pct_lb == 1) %>%
  mutate(data_trunc = map(data, function(x) x[x$sim > 0, c("sim", "N", "lb_goal", "ub_goal")])) %>%
  unnest(data_trunc) %>%
  select(scenario:pct_lb, power, sim:ub_goal) %>%
  filter(sim < 100 | sim >= 107) %>%
  mutate(regime = ifelse(sim < 100, "Historic", ifelse(sim >= 107, "Low", NA))) %>%
  group_by(scenario, regime, rep) %>%
  mutate(
    roll5 = rollmean(N, k = 5, align = "right", fill = NA),
    SOC_roll5 = roll5 < lb_goal, #rolling mean to account for the size of the miss/make
    mean_roll5 = mean(SOC_roll5, na.rm = TRUE),
    lb_N = N < lb_goal,
    SOC_N155 = get_SOC2(lb_N, 1, 5, 5), #most sensitive
    mean_N155 = mean(SOC_N155, na.rm = TRUE),
    SOC_N445 = get_SOC2(lb_N, 4, 4, 5), #ADF&G most common?
    mean_N445 = mean(SOC_N445, na.rm = TRUE),
    SOC_N515 = get_SOC2(lb_N, 5, 1, 5), #lease sensitive
    mean_N515 = mean(SOC_N515, na.rm = TRUE)) %>%
  group_by(scenario, rep)

#demonstrate how each criteria works
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, phi == 0, pct_MSY == 0.9, gamma == 1, rep == 1) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, sim, N, lb_goal, starts_with("SOC"), roll5) %>%
  pivot_longer(starts_with("SOC"), 
               names_to = "criteria", 
               names_pattern =  "SOC_(.*)", 
               values_to = "SOC") %>%
  mutate(roll5 = ifelse(criteria == "roll5", roll5, NA)) %>%
  filter(sim >= (202-9) & sim <= 220) %>%
  arrange(criteria) %>% # WARNING: needed for the coloration of the roll 5 line to lean the bars 
  ggplot(aes(x = as.character(sim), y = N, fill = SOC)) +
  geom_bar(stat = "identity", alpha = 0.5) +
  geom_hline(aes(yintercept = lb_goal)) +
  geom_line(aes(y = roll5, group = 1, color = as.numeric(lead(SOC))),
            linewidth = 2,
            show.legend = FALSE) +
  scale_color_gradient(low = col[1], high = col[2]) +
  scale_fill_manual(values = col) +
  facet_grid(criteria ~ .) +
  labs(x = "Simulation Year", 
       y = "Total Run", 
       fill = "Stock of Concern?", 
       title = paste0("ln(\u03B1)= ", 1.5, ", \u03C3= ", 0.5, ", \u03B3= ", 1, ", \u03A6= ", 0, ", %MSY= ", 90))


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
  pivot_longer(dplyr::starts_with("mean_"), 
               names_to = "Chronic_I",
               names_pattern = "mean_(.*)",
               values_to = "pct_SOC") %>%
  ggplot(aes(x = Chronic_I, y  = pct_SOC, color = regime)) + 
  geom_boxplot() +
  scale_color_discrete(name = "Regime", 
                     labels = c("Historic" = "Historic", "Low" = "Low")) +
  labs(x = "Chronic Inability Criteria",
       y = "Probability") +
  facet_grid(paste0("%MSY: ", pct_MSY) ~ paste0("\u03B3: ", gamma)) +
  labs(title = paste0("ln(\u03B1): ", 1.5, ", \u03C3: ", 0.5, ", \u03A6: ", 0))

####
# Look for sensitivity wrt 445 and roll5 across SR parameter combinations
# pct_MSY, phi and gamma
# phi has little effect
dat_SOC %>%
  filter(regime == "Low", sigma == 0.5, lnalpha_1 == 1.5) %>%
  pivot_longer(dplyr::starts_with("mean_"), 
               names_to = "Chronic_I", 
               names_pattern = "mean_(.*)",
               values_to = "pct_SOC") %>%
  filter(Chronic_I %in% c("N445", "roll5")) %>%
  ggplot(aes(x = Chronic_I, y  = pct_SOC, color = as.character(phi))) + 
  geom_boxplot() +
  labs(title = paste0("ln(\u03B1): ", 1.5, ", \u03C3: ", 0.5),
       color = "Autocorrelation: \u03A6",
       x = "Chronic Inability Criteria",
       y = "Probability") +
  facet_grid(paste0("%MSY: ", pct_MSY) ~ paste0("\u03B3: ", gamma))

# Look for sensitivity wrt 445 and roll5
# pct_MSY, lnalpha and gamma
# increasing lnalpha decreases sensitivity w depensatory dynamics
dat_SOC %>%
  filter(regime == "Low", phi == 0, sigma == 0.5) %>%
  pivot_longer(dplyr::starts_with("mean_"), 
               names_to = "Chronic_I", 
               names_pattern = "mean_(.*)",
               values_to = "pct_SOC") %>%
  filter(Chronic_I %in% c("N445", "roll5")) %>%
  ggplot(aes(x = Chronic_I, y  = pct_SOC, color = as.character(lnalpha_1))) + 
  geom_boxplot() +
  labs(title = paste0("%MSY: ", 0.9, ", \u03C3: ", 0.5),
       color = "Productivity: ln(\u03B1)",
       x = "Chronic Inability Criteria",
       y = "Probability") +
  facet_grid(paste0("%MSY: ", pct_MSY) ~ paste0("\u03B3: ", gamma))

# Look for sensitivity wrt 445 and roll5
# pct_MSY, sigma and gamma
# increased sigma decreases sensitivity
dat_SOC %>%
  filter(regime == "Low", phi == 0, lnalpha_1 == 1.5) %>%
  pivot_longer(dplyr::starts_with("mean_"), 
               names_to = "Chronic_I", 
               names_pattern = "mean_(.*)",
               values_to = "pct_SOC") %>%
  filter(Chronic_I %in% c("N445", "roll5")) %>%
  ggplot(aes(x = Chronic_I, y  = pct_SOC, color = as.character(sigma))) + 
  geom_boxplot() +
  labs(title = paste0("phi: ", 0, ", ln(\u03B1): ", 1.5), 
       color = "process error: \u03C3",
       x = "Chronic Inability Criteria",
       y = "Probability") +
  facet_grid(paste0("%MSY: ", pct_MSY) ~ paste0("\u03B3: ", gamma))








# Demonstrate differences between SOC 445 and SOC_roll5
# My conclusion is that thy are comparable
# roll 5 may be better at targeting actual periods of low abundance...
# but it is hard to manage to with SEG

# when N445 more sensitive....
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, phi == 0.8, pct_MSY == 0.9, gamma == 1, rep == 3) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, sim, N, lb_goal, SOC_N445, roll5, SOC_roll5) %>%
  filter(SOC_N445 > SOC_roll5) %>% print(n = 100)
 
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, phi == .8, pct_MSY == 0.9, gamma == 1, rep == 3) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, sim, N, 
         lb_goal, ub_goal, starts_with("SOC"), roll5) %>%
  pivot_longer(starts_with("SOC"), 
               names_to = "criteria", 
               names_pattern =  "SOC_(.*)", 
               values_to = "SOC") %>%
  mutate(roll5 = ifelse(criteria == "roll5", roll5, NA)) %>%
  filter(sim >= (218-5) & sim <= 226, criteria %in% c("N445", "roll5")) %>%
  arrange(criteria) %>% 
  ggplot(aes(x = as.character(sim), y = N, fill = SOC)) +
  geom_bar(stat = "identity", alpha = 0.5) +
 annotate("rect", xmin = 0, xmax = Inf, ymin = 4000, ymax = 10000, alpha = .2) +
  #geom_rect(aes(ymin = lb_goal, ymax = ub_goal, xmin = 0, xmax = Inf), fill = "lightgrey", alpha = 0.1) +
  geom_line(aes(y = roll5, group = 1, color = as.numeric(lead(SOC))), 
            linewidth = 3, 
            show.legend = FALSE) +
  scale_color_gradient(low = col[1], high = col[2]) +
  scale_fill_manual(values = col) +
  facet_grid(criteria ~ .) +
  theme_bw() +
  labs(x = "Simulation Year", 
       y = "Total Run", 
       fill = "Stock of Concern?", 
       title = paste0("ln(\u03B1)= ", 1.5, ", \u03C3= ", 0.5, ", \u03B3= ", 1, ", \u03A6= ", 0.8, ", %MSY= ", 90))

# when roll5 more sensitive....
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, phi == 0.8, pct_MSY == 0.9, gamma == 1, rep == 1) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, sim, N, lb_goal, SOC_N445, roll5, SOC_roll5) %>%
  filter(SOC_N445 < SOC_roll5) %>% print(n = 100)

#N45 more sensitive
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.5, phi == .8, pct_MSY == 0.9, gamma == 1, rep == 1) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, sim, N, 
         lb_goal, ub_goal, starts_with("SOC"), roll5) %>%
  pivot_longer(starts_with("SOC"), 
               names_to = "criteria", 
               names_pattern =  "SOC_(.*)", 
               values_to = "SOC") %>%
  mutate(roll5 = ifelse(criteria == "roll5", roll5, NA)) %>%
  filter(sim >= (224-5) & sim <= 230, criteria %in% c("N445", "roll5")) %>%
  arrange(criteria) %>% 
  ggplot(aes(x = as.character(sim), y = N, fill = SOC)) +
  geom_bar(stat = "identity", alpha = 0.5) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 3840, ymax = 10347, alpha = .2) +
  #geom_rect(aes(ymin = lb_goal, ymax = ub_goal, xmin = 0, xmax = Inf), fill = "lightgrey", alpha = 0.1) +
  geom_line(aes(y = roll5, group = 1, color = as.numeric(lead(SOC))), 
            linewidth = 3, 
            show.legend = FALSE) +
  scale_color_gradient(low = col[1], high = col[2]) +
  scale_fill_manual(values = col) +
  facet_grid(criteria ~ .) +
  theme_bw() +
  labs(x = "Simulation Year", 
       y = "Total Run", 
       fill = "Stock of Concern?", 
       title = paste0("ln(\u03B1)= ", 1.5, ", \u03C3= ", 0.5, ", \u03B3= ", 1, ", \u03A6= ", 0.8, ", %MSY= ", 90))

# confused result....
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.8, phi == 0, pct_MSY == 0.9, gamma == 1.3, rep == 20) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, sim, N, lb_goal, SOC_N445, roll5, SOC_roll5) %>%
  filter(SOC_N445 < SOC_roll5) %>% print(n = 100)

dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.8, phi == 0, pct_MSY == 0.9, gamma == 1.3, rep == 20) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, sim, N, 
         lb_goal, ub_goal, starts_with("SOC"), roll5) %>%
  pivot_longer(starts_with("SOC"), 
               names_to = "criteria", 
               names_pattern =  "SOC_(.*)", 
               values_to = "SOC") %>%
  mutate(roll5 = ifelse(criteria == "roll5", roll5, NA)) %>%
  filter(sim >= (165-5) & sim <= 182, criteria %in% c("N445", "roll5")) %>%
  arrange(criteria) %>% 
  ggplot(aes(x = as.character(sim), y = N, fill = SOC)) +
  geom_bar(stat = "identity", alpha = 0.5) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 4000, ymax = 10000, alpha = .2) +
  #geom_rect(aes(ymin = lb_goal, ymax = ub_goal, xmin = 0, xmax = Inf), fill = "lightgrey", alpha = 0.1) +
  geom_line(aes(y = roll5, group = 1, color = as.numeric(lead(SOC))), 
            linewidth = 3, 
            show.legend = FALSE) +
  scale_color_gradient(low = col[1], high = col[2]) +
  scale_fill_manual(values = col) +
  facet_grid(criteria ~ .) +
  theme_bw() +
  labs(x = "Simulation Year", 
       y = "Total Run", 
       fill = "Stock of Concern?", 
       title = paste0("ln(\u03B1)= ", 1.5, ", \u03C3= ", 0.8, ", \u03B3= ", 1.3, ", \u03A6= ", 0, ", %MSY= ", 90))


# confused result....
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.8, phi == 0, pct_MSY == 0.9, gamma == 1.6, rep == 3) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, sim, N, lb_goal, SOC_N445, roll5, SOC_roll5) %>%
  filter(SOC_N445 != SOC_roll5) %>% print(n = 100)
dat_SOC %>%
  filter(lnalpha_1 == 1.5, sigma == 0.8, phi == 0, pct_MSY == 0.9, gamma == 1.6, rep == 3) %>%
  select(scenario, rep, lnalpha_1, sigma, gamma, pct_MSY, sim, N, 
         lb_goal, ub_goal, starts_with("SOC"), roll5) %>%
  pivot_longer(starts_with("SOC"), 
               names_to = "criteria", 
               names_pattern =  "SOC_(.*)", 
               values_to = "SOC") %>%
  mutate(roll5 = ifelse(criteria == "roll5", roll5, NA)) %>%
  filter(sim >= (135-5) & sim <= 156, criteria %in% c("N445", "roll5")) %>%
  arrange(criteria) %>% 
  ggplot(aes(x = as.character(sim), y = N, fill = SOC)) +
  geom_bar(stat = "identity", alpha = 0.5) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 4000, ymax = 10000, alpha = .2) +
  #geom_rect(aes(ymin = lb_goal, ymax = ub_goal, xmin = 0, xmax = Inf), fill = "lightgrey", alpha = 0.1) +
  geom_line(aes(y = roll5, group = 1, color = as.numeric(lead(SOC))), 
            linewidth = 3, 
            show.legend = FALSE) +
  scale_color_gradient(low = col[1], high = col[2]) +
  scale_fill_manual(values = col) +
  facet_grid(criteria ~ .) +
  theme_bw() +
  labs(x = "Simulation Year", 
       y = "Total Run", 
       fill = "Stock of Concern?", 
       title = paste0("ln(\u03B1)= ", 1.5, ", \u03C3= ", 0.8, ", \u03B3= ", 1.6, ", \u03A6= ", 0, ", %MSY= ", 90))
