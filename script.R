# Simulate SOC listing criteria under a variety of productivity regimes

# Author: Adam Reimer
# Version: 2023-09-15

# Packages
packs <- c("tidyverse", "ggforce") #, "RcppRoll", "knitr", "ggforce")
lapply(packs, require, character.only = TRUE)

# source functions
function_files <- list.files(path=".\\functions")
lapply(function_files, function(x) source(paste0(".\\functions\\", x)))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Inputs -----------------------------------------------------------
# Is there concern when a stock falls below the escapement goal wo a change in SR parameters?
# How often does that occur?

# age-at-maturity 3-7 years 
#Kenai Chinook posterior for finding reasonable population parameter estimates
# post <- readRDS(file = "S:\\RTS\\Reimer\\KenaiSRA\\posts\\post_run_test3.rds")
# lnalpha.y <- post$q50$lnalpha.vec[,2]
# median(lnalpha.y[(length(lnalpha.y)-4):length(lnalpha.y)])
# median(lnalpha.y[(length(lnalpha.y)-9):length(lnalpha.y)])

# Grid of SR parameter values
vec_lnalpha <- seq(1, 2, length.out = 3)
vec_sigW <- seq(0.25, 0.75, length.out = 3)
vec_phi <- seq(0, .8, length.out = 3)
beta <- 0.0001
df_power <- data.frame(lnalpha = vec_lnalpha,
                       power = c(0.5, 0.6, 0.8))
input <- 
  expand.grid(lnalpha = vec_lnalpha, sigW = vec_sigW, phi = vec_phi) %>%
  mutate(beta = beta, 
         Smsy = get_Smsy(lnalpha, beta, TRUE, sigma = sigW, phi = phi),
         lb_p = Smsy * 0.8,
         ub_p = Smsy * 1.6) %>%
  rowwise() %>%
  mutate(lb = optimise(get_bounds, 
                       1:Smsy, 
                       lnalpha = lnalpha, 
                       beta = beta,
                       pct_MSY = 0.8,
                       correct = TRUE,
                       sigma = sigW,
                       phi = phi)$minimum,
          ub = optimise(get_bounds, 
                        Smsy:(Smsy*5), 
                        lnalpha = lnalpha, 
                        beta = beta,
                        pct_MSY = 0.8,
                        correct = TRUE,
                        sigma = sigW,
                        phi = phi)$minimum) %>%
  left_join(df_power, by = "lnalpha")

Chinook_age <- c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02)

input_sims = 1500


# Simulation params -------------------------------------------------------
Ricker_plots <- list()
for(i in 1:length(unique(input$lnalpha))){
Ricker_plots[[i]] <- 
  input %>% 
    slice(rep(1:n(), each = 151)) %>% 
    mutate(S = rep(seq(0, 3 * 1 / max(beta), by = 200), times = nrow(input)),
           R = Ricker(lnalpha = lna_p, beta = beta, S = S),
           S_msy = lb_p / 0.8,
           R_msy = Ricker(lna_p, beta, S_msy),
           S_max = 1 / beta,
           R_max = Ricker(lnalpha = lna_p, beta = beta, S = S_max)) %>%
    ggplot(aes(x = S, y = R)) +
      geom_rect(aes(xmin = lb_p, xmax = ub_p, ymin = -Inf, ymax = Inf), fill = "grey95") +
      geom_line() +
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      geom_point(aes(x = S_max, y = R_max)) +
      geom_segment(aes(x = S_msy, xend = S_msy, y = S_msy, yend = R_msy), linetype = 3) +
      theme_bw() +
      facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                          #scales = "free_y",
                          ncol = 3, nrow = 3, page = i)
}
Ricker_plots

# Base case -------------------------------------------------------------------
#First simulation sigN = 0
#Not very realistic
# Simulation results
sim_base_list <- mapply(FUN = simSR_goal,
                        lnalpha = input$lnalpha,
                        sigW = input$sigW,
                        phi = input$phi,
                        lb_goal = input$lb_p,
                        ub_goal = input$ub_p,
                        lb_manage = input$lb_p,
                        ub_manage = input$ub_p,
                        MoreArgs = list(beta = beta,
                                        age0 = Chinook_age,
                                        Sims = input_sims,
                                        sigN = 0.2,
                                        sigF = 0,
                                        Hfun = H_goal),
                        SIMPLIFY = FALSE)
saveRDS(sim_base_list, file = ".\\sim_base_list.R")
# sim_base_list <- readRDS(file = ".\\sim_base_list.R")

# * Create dataframe --------------------------------------------------------
sim_base_df <-
  lapply(sim_base_list, as.data.frame) %>%
  lapply(function(x){mutate(x, sim = row_number())}) %>%
  do.call("rbind", .) %>%
  filter(R != 0) %>% #remove rows where population perished
  mutate(resid = S - lb_goal) %>%  
  select(-starts_with("N_age")) %>%
  group_by(lnalpha, sigW, phi) %>%
  mutate(resid_roll = roll_mean(x = resid, 5, align = "right", fill = NA) / lb_goal)

# * Plot time series --------------------------------------------------------
#writing function here since it has dubious use outside of this single application
plot_ts <- function(dat, lnalpha0){ 
  df <- 
    dat %>%
    select(lnalpha, sigW, phi, sim, lb_goal, ub_goal, S, N) %>%
    filter(lnalpha == lnalpha0) %>%
    pivot_longer(c(S, N), names_to = "stat", values_to = "value")
  
  limit <- quantile(df$value[df$stat == "N"], .95, na.rm = TRUE)
  
      ggplot(df, aes(x = sim, y = value, color = stat)) +
        geom_line(alpha = 0.4) +
        geom_hline(aes(yintercept = lb_goal), linetype = 2) +
        geom_hline(aes(yintercept = ub_goal), linetype = 2) +
        scale_y_continuous(limits = c(0, limit)) +
        scale_x_continuous(name = "Simulation year") +
        facet_grid(paste0("\u03C6: ", phi) ~ paste0("\u03C3: ", sigW)) +
        ggtitle(label = bquote(log(alpha): .(lnalpha0)))
}
# because of the colors used teal below the lower bound represents times we fished below the lower bound w abundant fish... call the a yield concern
# while grey below the lower bound represent when the fish were not there to make the goal ... call that a management concern
# When N is far below the lb it is a conservation concern
# Note that some phi and sigma W combinations are not sustainable

plot_ts(sim_base_df, 1)
plot_ts(sim_base_df, 1.5)
plot_ts(sim_base_df, 2)



# * Table of concern criteria -----------------------------------------------
sim_base_df %>%
  summarise(N = length(sim),
            below_lb_fish = sum(cc, na.rm = TRUE) / N,
            below_lb_harvest = sum(mc, na.rm = TRUE)/ N,
            above_ub = sum(yc, na.rm = TRUE)/ N,
            cc4.5 = sum(SOC == "Conservation", na.rm = TRUE)/ N,
            mc4.5 = sum(SOC == "Management", na.rm = TRUE)/ N,
            yc4.5 = sum(SOC == "Yield", na.rm = TRUE)/ N,
            resid = sum(resid_roll < 0, na.rm = TRUE)/ N) %>%
  arrange(phi, sigW, lnalpha) %>%
  kable(digits = 2,
        col.names = c("ln(\u03B1)",
                      "\u03C3",
                      "\u03C6",
                      "N",
                      "P(N < lower bound / 2)",
                      "P(N < lb)",
                      "P(S < lb, N > lb)",
                      "Conservation concern",
                      "Management concern",
                      "Yield concern",
                      "Residual"))

# Criteria occurrence
sim_base_df %>%
  summarise(cc = sum(SOC == "Conservation"),
            mc = sum(SOC == "Management"),
            yc = sum(SOC == "Yield")) %>%
  pivot_longer(cols = c("cc", "mc", "yc"), 
              names_to = "SOC", 
              values_to = "N") %>%
  ggplot(aes(x = sigW, y = N, color = SOC)) +
  geom_line() +
  facet_grid(phi ~ lnalpha, labeller = label_bquote(rows = phi: .(phi), cols = log(alpha): .(lnalpha)))


# * Average 5 year rolling residual as a diagnostic -------------------------
plot_resid <- function(dat, lnalpha0){
  dat %>%
    filter(lnalpha == lnalpha0, !is.na(SOC)) %>%
    mutate(resid_pct = resid / lb_goal) %>%
    ggplot(aes(x = resid_pct, fill = SOC)) +
    geom_histogram() +
    scale_x_continuous(limits = c(-2, 2)) +
    geom_vline(aes(xintercept = 0)) +
    facet_grid(paste0("\u03C6: ", phi) ~ paste0("\u03C3: ", sigW),
               scales = "free_y") +
    ggtitle(label = bquote(log(alpha): .(lnalpha0))) #.(lnalpha) pulls from the parent frame instead of the function's environment
}

#Very few concerns and when they do occur residual and "miss" diagnostics mostly agree.
plot_resid(sim_base_df, 1)
plot_resid(sim_base_df, 1.5)
plot_resid(sim_base_df, 2)


# Reduced alpha: lb = Seq  -------------------------------------
# #What if productivity changes but our goal remains based on historic productivity?
# #Note S_eq = lnalpha/beta so lnalpha = goal_lb*beta is a productivity change that would reduce yield to 0 at the lower bound.
# # lower bound goal which is 0.8% of Smsy
lnalpha <- 1.5
sigW <- 0.5
lb <- lnalpha/beta*(0.5 - 0.07 * lnalpha) * 0.8
# #Assume goal set using corrected lnalpha
lnalpha_p <- lnalpha + sigW * sigW / 2
lb_p <- lnalpha_p/beta*(0.5 - 0.07 * lnalpha_p) * 0.8
ub_p <- lnalpha_p/beta*(0.5 - 0.07 * lnalpha_p) * 1.6

# example
S<- 1:30000
plot(S, S*exp(lnalpha - beta*S))
abline(a = 0, b = 1)
polygon(x = c(lb_p, lb_p, ub_p, ub_p), 
        y = c(0, 17000, 17000, 0), 
        col = scales::alpha("grey", .3))
lines(S, S*exp(lb_p*beta - beta*S), col = "red")
# 33% reduction in productivity in this case
lb_p*beta / lnalpha

plot(S, log(S*exp(lnalpha - beta*S) / S))
lines(S, log(S*exp(lb_p*beta - beta*S) / S), col = "red")
polygon(x = c(lb_p, lb_p, ub_p, ub_p), 
        y = c(-1.5, 1.5, 1.5, -1.5), 
        col = scales::alpha("grey", .3))
abline(a = 0, b = 0, lt = 2)

# * Simulation params -------------------------------------------------------
Ricker_plots_Seq <- list()
for(i in 1:length(unique(input$lnalpha))){
  Ricker_plots_Seq[[i]] <- 
    input %>% 
    slice(rep(1:n(), each = 151)) %>% 
    mutate(S = rep(seq(0, 3 * 1 / max(beta), by = 200), times = nrow(input)),
           R = Ricker(lnalpha = lna_p, beta = beta, S = S),
           S_msy = lb_p / 0.8,
           R_msy = Ricker(lna_p, beta, S_msy),
           S_max = 1 / beta,
           R_max = Ricker(lnalpha = lna_p, beta = beta, S = S_max),
           R_Seq = Ricker(lnalpha = lb_p * beta, beta = beta, S = S),
           R_max_Seq = Ricker(lnalpha = lb_p * beta, beta = beta, S = S_max)) %>%
    ggplot(aes(x = S, y = R)) +
    geom_rect(aes(xmin = lb_p, xmax = ub_p, ymin = -Inf, ymax = Inf), fill = "grey95") +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_point(aes(x = S_max, y = R_max)) +
    geom_segment(aes(x = S_msy, xend = S_msy, y = S_msy, yend = R_msy), linetype = 3) +
    geom_line(aes(y = R_Seq), color = "red") +
    geom_point(aes(x = S_max, y = R_max_Seq), color = "red") +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        #scales = "free_y",
                        ncol = 3, nrow = 3, page = i)
}
Ricker_plots_Seq[[1]] + coord_cartesian(ylim = c(0, 10000))
Ricker_plots_Seq[[2]] + coord_cartesian(ylim = c(0, 10000))
Ricker_plots_Seq[[3]] + coord_cartesian(ylim = c(0, 10000))


# * simulation reduced alpha -------------------------------------------------------------------
#grid simulation using same lower bound as above (average productivity lb) but reduced productivity
#average productivity
input$lna_p
#reduced productivity
input$lb_p  * beta
# % of average productivity ~ 14-37%
input$lb_p * beta / input$lna_p

#average productivity lower bound (used in simulation)
input$lb_p
#reduced productivity lower bound (not used in simulation)
input$lb_p*beta/beta * (0.5 - 0.07 * input$lb_p*beta) *.8
# % of average productivity lower bound ~ 38%
input$lb_p*beta/beta * (0.5 - 0.07 * input$lb_p*beta) *.8 / input$lb_p

#simulation results
# Simulation results
# sim_Seq_N1_list <- mapply(FUN = simSR_goal,
#                           lnalpha = input$lb_p*beta,  #reduced ln_alpha
#                           sigW = input$sigW,
#                           phi = input$phi,
#                           lb_goal = input$lb_p,
#                           ub_goal = input$ub_p,
#                           lb_manage = input$lb_p,
#                           ub_manage = input$ub_p,
#                           MoreArgs = list(beta = beta,
#                                           age0 = Chinook_age,
#                                           Sims = input_sims,
#                                           sigN = 0.1,
#                                           sigF = 0,
#                                            Hfun = H_goal),
#                           SIMPLIFY = FALSE)
# saveRDS(sim_Seq_N1_list, file = ".\\sim_Seq_N1_list.R")
sim_Seq_N1_list <- readRDS(file = ".\\sim_Seq_N1_list.R")


# * Create dataframe --------------------------------------------------------
sim_Seq_N1_df <-
  lapply(sim_Seq_N1_list, as.data.frame) %>%
  lapply(function(x) rename(x, lnalpha_red = lnalpha)) %>%
  lapply(function(x){mutate(x, sim = row_number())}) %>%
  do.call("rbind", .) %>%
  mutate(lnalpha = rep(input$lnalpha, each = input_sims)) %>%
  filter(R != 0) %>% #remove rows where population perished
  mutate(resid = S - lb_goal) %>%  
  select(-starts_with("N_age")) %>%
  group_by(lnalpha, sigW, phi) %>%
  mutate(resid_roll = roll_mean(x = resid, 5, align = "right", fill = NA) / lb_goal)



# * Plot time series --------------------------------------------------------
# because of the colors used teal below the lower bound represents times we fished below the lower bound w abundant fish... call the a management concern
# teal above the upper bound represents underutilized yield... call that a yield concern
# while grey below the lower bound represent when the fish were not there to make the goal ... call that a conservation concern
# Note that some phi and sigma W combinations are not sustainable
plot_ts(sim_Seq_N1_df, 1)
plot_ts(sim_Seq_N1_df, 1.5)
plot_ts(sim_Seq_N1_df, 2)



# * Table of concern criteria -----------------------------------------------
# estimation error for N introduces lb misses based on fishing and increased ub misses.
sim_Seq_N1_df %>%
  summarise(N = length(sim),
            below_lb_fish = sum(cc, na.rm = TRUE) / N,
            below_lb_harvest = sum(mc, na.rm = TRUE) / N,
            above_ub = sum(yc, na.rm = TRUE) / N,
            cc4.5 = sum(SOC == "Conservation", na.rm = TRUE) / N,
            mc4.5 = sum(SOC == "Management", na.rm = TRUE) / N,
            yc4.5 = sum(SOC == "Yield", na.rm = TRUE) / N,
            resid = sum(resid_roll < 0, na.rm = TRUE)) %>%
  arrange(phi, sigW, lnalpha)  %>%
  kable(digits = 2,
        col.names = c("ln(\u03B1)",
                      "\u03C3",
                      "\u03C6",
                      "N",
                      "P(N < lower bound / 2)",
                      "P(N < lb)",
                      "P(S < lb, N > lb)",
                      "Conservation concern",
                      "Management concern",
                      "Yield concern",
                      "Residual"))

# * Criteria occurrence ####
sim_Seq_N1_df %>%
  summarise(cc = sum(SOC == "Conservation"),
            mc = sum(SOC == "Management"),
            yc = sum(SOC == "Yield")) %>%
  pivot_longer(cols = c("cc", "mc", "yc"), 
               names_to = "SOC", 
               values_to = "N") %>%
  ggplot(aes(x = sigW, y = N, color = SOC)) +
  geom_line() +
  facet_grid(phi ~ lnalpha, labeller = label_bquote(rows = phi: .(phi), cols = log(alpha): .(lnalpha)))

# * Average 5 year rolling residual as a diagnostic -------------------------
#Many concerns.
plot_resid(sim_Seq_N1_df, 1)
plot_resid(sim_Seq_N1_df, 1.5)
plot_resid(sim_Seq_N1_df, 2)





# Rebuild -----------------------------------------------------------------
# could be used to simulate a rebuild
sim_Seq_rebuild_list <- mapply(FUN = simSR_goal, 
                          lnalpha = input$lb_p*beta,  #reduced ln_alpha 
                          sigW = input$sigW, 
                          phi = input$phi,
                          lb_goal = input$lb_p,
                          ub_goal = input$ub_p,
                          lb_manage = input$ub_p,
                          ub_manage = input$ub_p,
                          MoreArgs = list(beta = beta, 
                                          age0 = Chinook_age,
                                          Sims = input_sims, 
                                          sigN = 0.1,
                                          sigF = 0,
                                          Hfun = H_soc),
                          SIMPLIFY = FALSE)

# * Create dataframe --------------------------------------------------------
sim_Seq_rebuild_df <-
  lapply(sim_Seq_rebuild_list, as.data.frame) %>%
  lapply(function(x) rename(x, lnalpha_red = lnalpha)) %>%
  lapply(function(x){mutate(x, sim = row_number())}) %>%
  do.call("rbind", .) %>%
  mutate(lnalpha = rep(input$lnalpha, each = input_sims)) %>%
  filter(R != 0) %>% #remove rows where population perished
  mutate(resid = S - lb_goal) %>%  
  select(-starts_with("N_age")) %>%
  group_by(lnalpha, sigW, phi) %>%
  mutate(resid_roll = roll_mean(x = resid, 5, align = "right", fill = NA) / lb_goal)



# * Plot time series --------------------------------------------------------
# because of the colors used teal below the lower bound represents times we fished below the lower bound w abundant fish... call the a management concern
# teal above the upper bound represents underutilized yield... call that a yield concern
# while grey below the lower bound represent when the fish were not there to make the goal ... call that a conservation concern
# Note that some phi and sigma W combinations are not sustainable
plot_ts(sim_Seq_rebuild_df, 1)
plot_ts(sim_Seq_rebuild_df, 1.5)
plot_ts(sim_Seq_rebuild_df, 1.5) + coord_cartesian(xlim = c(600, 625))
plot_ts(sim_Seq_rebuild_df, 2)


# * Table of concern criteria -----------------------------------------------
sim_Seq_rebuild_df %>%
  summarise(N = length(sim),
            below_lb_fish = sum(cc, na.rm = TRUE) / N,
            below_lb_harvest = sum(mc, na.rm = TRUE) / N,
            above_ub = sum(yc, na.rm = TRUE) / N,
            cc4.5 = sum(SOC == "Conservation", na.rm = TRUE) / N,
            mc4.5 = sum(SOC == "Management", na.rm = TRUE) / N,
            yc4.5 = sum(SOC == "Yield", na.rm = TRUE) / N,
            resid = sum(resid_roll < 0, na.rm = TRUE)) %>%
  arrange(phi, sigW, lnalpha)  %>%
  kable(digits = 2,
        col.names = c("ln(\u03B1)",
                      "\u03C3",
                      "\u03C6",
                      "N",
                      "P(N < lb)",
                      "P(S < lb, N > lb)",
                      "P(S > ub)",
                      "Conservation concern",
                      "Management concern",
                      "Yield concern",
                      "Residual"))

#an attempt at summary statistics for statuis Quo and Rebuilding fisheries under reduced productivity.
temp1 <- sim_Seq_N1_df %>%
  mutate(scenario = "StatusQuo")
temp2 <- sim_Seq_rebuild_df %>%
  mutate(scenario = "Rebuild")
# temp3 <- 
#   rbind(temp2, temp1) %>%
#   group_by(lnalpha, sigW, phi, scenario) %>%
#   select(-starts_with(c("lb", "ub"))) %>% 
#   mutate(change = case_when(SOC != lag(SOC) ~ TRUE, TRUE ~ FALSE),
#          n_change = cumsum(change),
#          H = N - S) %>%
#   group_by(lnalpha, sigW, phi, scenario, n_change, SOC) %>%
#   summarise(length = length(n_change),
#             deviation_lb = sum(resid),
#             miss = sum(resid < 0),
#             yield = sum(H),
#             min_S = min(S)) %>%
#   group_by(lnalpha, sigW, phi, scenario, SOC) %>%
#   summarise(N = length(n_change),
#             mean_years = median(length),
#             sd_years = sd(length),
#             mean_dev = median(deviation_lb),
#             sd_dev = sd(deviation_lb),
#             mean_yield = median(yield),
#             sd_yield = sd(yield),
#             mean_minS = median(min_S),
#             sd_minS = sd(min_S),
#             length = sum(length), 
#             mean_miss = sum(miss) / length,
#             sd_miss = mean_miss * (1 - mean_miss) / N, 
#             n_yield = sum(yield > 1))
#   
# 
# scenario_plots <- function(dat, var){
#   plots <- list()
#   for(i in 1:length(unique(input$lnalpha))){
#     plots[[i]] <-
#       dat %>%
#       mutate(upper = get(paste0("mean_", var)) + get(paste0("sd_", var)),
#              lower = get(paste0("mean_", var)) - get(paste0("sd_", var))) %>%
#       ggplot(aes(x = SOC, y = get(paste0("mean_", var)), weigth = N, fill = scenario)) +
#       geom_bar(stat = "identity", position = "dodge") +
#       geom_errorbar(aes(ymin = lower, ymax = upper), position = "dodge") +
#       theme_bw() +
#       facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
#                           ncol = 3, nrow = 3, page = i)
#   }
#   plots
# }
# length_plots <- scenario_plots(temp3, "years")
# length_plots[[1]]
# length_plots[[2]]
# length_plots[[3]]
# 
# dev_plots <- scenario_plots(temp3, "dev")
# dev_plots[[1]]
# dev_plots[[2]]
# dev_plots[[3]]
# 
# temp3 %>%
#   filter(SOC == "Conservation") %>%
#   select(lnalpha, sigW, phi, scenario, SOC, ends_with("dev")) %>%
#   pivot_wider(id_cols = c(lnalpha, sigW, phi, SOC), names_from = scenario,  values_from = ends_with("dev")) %>%
#   mutate(smaller =  mean_dev_StatusQuo - mean_dev_Rebuild) %>%
#   print(n = 100)
# 
# yield_plots <- scenario_plots(temp3, "yield")
# yield_plots[[1]]
# yield_plots[[2]]
# yield_plots[[3]]
# 
# minS_plots <- scenario_plots(temp3, "minS")
# minS_plots[[1]]
# minS_plots[[2]]
# minS_plots[[3]]
# 
# miss_plots <- scenario_plots(temp3, "miss")
# miss_plots[[1]]
# miss_plots[[2]]
# miss_plots[[3]]

#Similar but displays variability using box plots.
temp3 <- 
  rbind(temp2, temp1) %>%
  group_by(lnalpha, sigW, phi, scenario) %>%
  #select(-starts_with(c("lb", "ub"))) %>% 
  mutate(change = case_when(SOC != lag(SOC) ~ TRUE, TRUE ~ FALSE),
         n_change = cumsum(change),
         H = N - S,
         deviation_lb = resid / lb_goal) %>%
  group_by(lnalpha, sigW, phi, scenario, n_change, SOC) %>%
  summarise(length = length(n_change),
            deviation_lb = mean(deviation_lb),
            miss = sum(resid < 0) / length,
            yield = sum(H),
            min_S = min(S))


scenario_plots2 <- function(dat, var){
  plots <- list()
  for(i in 1:length(unique(input$lnalpha))){
    plots[[i]] <-
      dat %>%
      ggplot(aes(x = SOC, y = get(var), fill = scenario)) +
      geom_boxplot(position = "dodge") +
      theme_bw() +
      facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                          ncol = 3, nrow = 3, page = i)
  }
  plots
}
length_plots <- scenario_plots2(temp3, "length")
length_plots[[1]]
length_plots[[2]]
length_plots[[3]]

dev_plots <- scenario_plots2(temp3[temp3$SOC == "Conservation", ], "deviation_lb")
dev_plots[[1]]
dev_plots[[2]]
dev_plots[[3]]

yield_plots <- scenario_plots2(temp3, "yield")
yield_plots[[1]]
yield_plots[[2]]
yield_plots[[3]]

minS_plots <- scenario_plots2(temp3[temp3$SOC == "Conservation", ], "min_S")
minS_plots[[1]]
minS_plots[[2]]
minS_plots[[3]]

miss_plots <- scenario_plots2(temp3[temp3$SOC == "Conservation", ], "miss")
miss_plots[[1]]
miss_plots[[2]]
miss_plots[[3]]



#Shnute params
#Less parameter correlation during estimation
#Can. J. Fish. Aquat. Sci. 53: 12811293 (1996).
#h: harvest rate
#C: MSY

#to get h and C from alpha and beta
h <- rep(NA, 101)
alpha <- rep(NA, 100)
c <- rep(NA, 100)

for(i in 1:100){
  h[1] <- 0.5
  alpha[i] <- exp(h[i])/(1 - h[i])
  h[i + 1] <- h[i] + (1 - h[i])/(2 - h[i])*log(exp(1.5)/alpha[i]) #alpha = 1.5
}
h[100] #h
h[100]^2 / ( 1 - h[100]) / 0.0001 #C, beta = 0.0001

#Plot equivalence
ricker <- function(s) s*exp(1.5 - 0.0001*s)
ricker_alt <- function(s, h = 0.64, c = 8484) s/(1 - h) * exp(h - h^2/(1 - h)*s/c)
s <- 1:30000
plot(s, ricker(s), type = 'l') #alpha, beta Ricker
abline(0, 1)
Smsy <- 1.5 / 0.0001 * (0.5 - 0.07 * 1.5) #Smsy for the ricker above
Rmsy <- ricker(Smsy)
C <- Rmsy - Smsy
h <- C/Rmsy
lines(ricker_alt(s, h, C), col = "red") #h, C Ricker
#Find a value of h that provides negligible yielde for the same beta
f <- function(h, beta = 0.0001, C) {h^2 / ( 1 - h) / beta - C}
uniroot(f, c(0, 1), C = 100)
lines(ricker_alt(s, 0.095, 100), col = "green")

#Curve changes with constant C and varied h
test <-
  data.frame(C = 8757,
           h = seq(.4, .9, length.out = 4)) %>%
  mutate(a = exp(h) / (1 - h),
         b = h^2 / ((1 - h)*C))
test
plot(s, ricker_alt(s, test$h[1], test$C[1]))
abline(0, 1)
lines(s, ricker_alt(s, test$h[2], test$C[2]))
lines(s, ricker_alt(s, test$h[3], test$C[3]))
lines(s, ricker_alt(s, test$h[4], test$C[4]))

#Curve changes with constant h and varied C
test2 <-
  data.frame(C = c(100, 3000, 8757, 12000),
             h = .57) %>%
  mutate(a = exp(h) / (1 - h),
         b = h^2 / ((1 - h)*C))
test2

plot(s, ricker_alt(s, test2$h[4], test2$C[4]), col = "red")
abline(0, 1)
lines(s, ricker_alt(s, test2$h[1], test2$C[1]), col = "red")
lines(s, ricker_alt(s, test2$h[2], test2$C[2]), col = "red")
lines(s, ricker_alt(s, test2$h[3], test2$C[3]), col = "red")

#Neither of the above patterns seem realistic
#Curve changes with constant beta and varied C
#So can describe what I think is the most likely situation
#Would be interesting to see how it works in estimation.
test3 <-
  data.frame(C = c(100, 3000, 8757, 12000),
             h = c(uniroot(f, c(0, 1), C = 100)[[1]],
                   uniroot(f, c(0, 1), C = 3000)[[1]],
                   uniroot(f, c(0, 1), C = 8757)[[1]],
                   uniroot(f, c(0, 1), C = 12000)[[1]])) %>%
  mutate(a = exp(h) / (1 - h),
         b = h^2 / ((1 - h)*C))
test3

plot(s, ricker_alt(s, test3$h[4], test3$C[4]), col = "green")
abline(0, 1)
lines(s, ricker_alt(s, test3$h[1], test3$C[1]), col = "green")
lines(s, ricker_alt(s, test3$h[2], test3$C[2]), col = "green")
lines(s, ricker_alt(s, test3$h[3], test3$C[3]), col = "green")
