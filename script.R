# Simulate SOC listing criteria under a variety of productivity regimes

# Author: Adam Reimer
# Version: 2023-09-15

# Packages
packs <- c("tidyverse", "ggforce", "RcppRoll", "knitr") #, "ggforce")
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
#Note high combinations for sigW and phi lead to a population crash (herein defined as < 100 individuals)
vec_lnalpha <- seq(1, 2, length.out = 3)
vec_sigW <- seq(0.25, 0.75, length.out = 3)
vec_phi <- seq(0, .7, length.out = 3)
beta <- 0.0001
df_power <- data.frame(lnalpha = vec_lnalpha,
                       power = c(0.5, 0.6, 0.8)) #how to pick these. High for SOC sims. Should be lower for OPY sims.
input <- 
  expand.grid(lnalpha = vec_lnalpha, sigW = vec_sigW, phi = vec_phi) %>%
  mutate(beta = beta, 
         Smsy = get_Smsy(lnalpha, beta, TRUE, sigma = sigW, phi = phi),
         lb_eggers = Smsy * 0.8, #Egger bounds
         ub_eggers = Smsy * 1.6) %>%
  rowwise() %>%
  mutate(lb_pctMSY = optimise(get_bounds, #'true' OYP bounds
                              1:Smsy, 
                              lnalpha = lnalpha, 
                              beta = beta,
                              pct_MSY = 0.9,
                              correct = TRUE,
                              sigma = sigW,
                              phi = phi)$minimum,
          ub_pctMSY = optimise(get_bounds, 
                              Smsy:(Smsy*5), 
                              lnalpha = lnalpha, 
                              beta = beta,
                              pct_MSY = 0.7,
                              correct = TRUE,
                              sigma = sigW,
                              phi = phi)$minimum) %>%
  ungroup() %>%
  left_join(df_power, by = "lnalpha")

Chinook_age <- c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02)

input_sims = 5000


# Simulation params -------------------------------------------------------
Ricker_plots <- list()
for(i in 1:length(unique(input$lnalpha))){
Ricker_plots[[i]] <- 
  input %>% 
    slice(rep(1:n(), each = 251)) %>% 
    mutate(S = rep(seq(0, 5 * 1 / max(beta), by = 200), times = nrow(input)),
           R = Ricker(lnalpha = get_lnalpha_p(lnalpha, beta, sigW, phi), beta = beta, S = S),
           Rmsy = Ricker(get_lnalpha_p(lnalpha, beta, sigW, phi), beta, Smsy),
           Smax = 1 / beta,
           Rmax = Ricker(lnalpha = get_lnalpha_p(lnalpha, beta, sigW, phi), beta = beta, S = Smax)) %>%
    ggplot(aes(x = S, y = R)) +
      geom_rect(aes(xmin = lb_pctMSY, xmax = ub_pctMSY, ymin = -Inf, ymax = Inf), fill = "grey95") +
      geom_line() +
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      geom_point(aes(x = Smax, y = Rmax)) +
      geom_segment(aes(x = Smsy, xend = Smsy, y = Smsy, yend = Rmsy), linetype = 3) +
      theme_bw() +
      facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                          #scales = "free_y",
                          ncol = 3, nrow = 3, page = i)
}
Ricker_plots[[1]] + coord_cartesian(xlim = c(0, 20000), ylim = c(0, 20000))
Ricker_plots[[2]] + coord_cartesian(xlim = c(0, 30000), ylim = c(0, 30000))
Ricker_plots[[3]] + coord_cartesian(xlim = c(0, 50000), ylim = c(0, 50000))


# Base case -------------------------------------------------------------------
# Simulation results
sim_base_list <- mapply(FUN = simSR_goal,
                        lnalpha = input$lnalpha,
                        sigW = input$sigW,
                        phi = input$phi,
                        lb_goal = input$lb_pctMSY,
                        ub_goal = input$ub_pctMSY,
                        power = input$power,
                        MoreArgs = list(beta = beta,
                                        age0 = Chinook_age,
                                        Sims0 = input_sims,
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
            yc4.5 = sum(SOC == "Yield", na.rm = TRUE)/ N) %>%
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
                      "Yield concern"))

# * Average 5 year rolling residual as a diagnostic -------------------------
# Promising idea to explore later
# plot_resid <- function(dat, lnalpha0){
#   dat %>%
#     filter(lnalpha == lnalpha0, !is.na(SOC)) %>%
#     mutate(resid_pct = resid / lb_goal) %>%
#     ggplot(aes(x = resid_pct, fill = SOC)) +
#     geom_histogram() +
#     scale_x_continuous(limits = c(-2, 2)) +
#     geom_vline(aes(xintercept = 0)) +
#     facet_grid(paste0("\u03C6: ", phi) ~ paste0("\u03C3: ", sigW),
#                scales = "free_y") +
#     ggtitle(label = bquote(log(alpha): .(lnalpha0))) #.(lnalpha) pulls from the parent frame instead of the function's environment
# }
# 
# #Very few concerns and when they do occur residual and "miss" diagnostics mostly agree.
# plot_resid(sim_base_df, 1)
# plot_resid(sim_base_df, 1.5)
# plot_resid(sim_base_df, 2)


# Reduced alpha: lb = Seq  -------------------------------------
# #What if productivity changes but our goal remains based on historic productivity?
# #Note S_eq = lnalpha/beta so lnalpha = goal_lb*beta is a productivity change that would reduce yield to 0 at the lower bound.
# #consider this a model for current Chinook productivity regimes
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
    slice(rep(1:n(), each = 251)) %>% 
    mutate(S = rep(seq(0, 5 * 1 / max(beta), by = 200), times = nrow(input)),
           R = Ricker(lnalpha = get_lnalpha_p(lnalpha, beta, sigW, phi), beta = beta, S = S),
           Rmsy = Ricker(get_lnalpha_p(lnalpha, beta, sigW, phi), beta, Smsy),
           Smax = 1 / beta,
           Rmax = Ricker(lnalpha = get_lnalpha_p(lnalpha, beta, sigW, phi), beta = beta, S = Smax),
           RSeq = Ricker(lnalpha = lb_pctMSY * beta, beta = beta, S = S),
           Rmax_Seq = Ricker(lnalpha = lb_pctMSY * beta, beta = beta, S = Smax)) %>%
    ggplot(aes(x = S, y = R)) +
    geom_rect(aes(xmin = lb_pctMSY, xmax = ub_pctMSY, ymin = -Inf, ymax = Inf), fill = "grey95") +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_point(aes(x = Smax, y = Rmax)) +
    geom_segment(aes(x = Smsy, xend = Smsy, y = Smsy, yend = Rmsy), linetype = 3) +
    geom_line(aes(y = RSeq), color = "red") +
    geom_point(aes(x = Smax, y = Rmax_Seq), color = "red") +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        #scales = "free_y",
                        ncol = 3, nrow = 3, page = i)
}
Ricker_plots_Seq[[1]] + coord_cartesian(xlim = c(0, 20000), ylim = c(0, 20000))
Ricker_plots_Seq[[2]] + coord_cartesian(xlim = c(0, 30000), ylim = c(0, 30000))
Ricker_plots_Seq[[3]] + coord_cartesian(xlim = c(0, 50000), ylim = c(0, 50000))


# * simulation reduced alpha -------------------------------------------------------------------
#grid simulation using same lower bound as above (average productivity lb) but reduced productivity
# % of average productivity ~ 15-20%
input %>% 
  mutate(lnalpha_p = get_lnalpha_p(lnalpha, beta, sigW, phi),
         lnalpha_p_red = lb_pctMSY  * beta,
         pct_red = lnalpha_p_red / lnalpha_p) %>%
  print(n = 100)

# % of average productivity lower bound ~ .47%
input %>% 
  mutate(lnalpha_p_red = lb_pctMSY  * beta,
         lb_pctMSY_red = get_Smsy(lnalpha_p_red, beta, FALSE),
         pct_red = lb_pctMSY_red / lb_pctMSY) %>%
  print(n = 100)

#simulation results
# Simulation results
sim_Seq_N1_list <- mapply(FUN = simSR_goal,
                          lnalpha = input$lb_pctMSY*beta,  #reduced ln_alpha
                          sigW = input$sigW,
                          phi = input$phi,
                          lb_goal = input$lb_pctMSY,
                          ub_goal = input$ub_pctMSY,
                          MoreArgs = list(beta = beta,
                                          age0 = Chinook_age,
                                          Sims0 = input_sims,
                                          sigN = 0.1,
                                          sigF = 0,
                                           Hfun = H_goal),
                          SIMPLIFY = FALSE)
saveRDS(sim_Seq_N1_list, file = ".\\sim_Seq_N1_list.R")
#sim_Seq_N1_list <- readRDS(file = ".\\sim_Seq_N1_list.R")


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
            yc4.5 = sum(SOC == "Yield", na.rm = TRUE) / N) %>%
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
                      "Yield concern"))

# * Average 5 year rolling residual as a diagnostic -------------------------
#Many concerns.
# plot_resid(sim_Seq_N1_df, 1)
# plot_resid(sim_Seq_N1_df, 1.5)
# plot_resid(sim_Seq_N1_df, 2)


# Rebuild -----------------------------------------------------------------
# Can managing to a higher goal during management or conservation concerns reduce the duration of concern designations or 
# increase escapement during concern designations
sim_Seq_rebuild_list <- mapply(FUN = simSR_goal, 
                          lnalpha = input$lb_pctMSY*beta,  #reduced ln_alpha 
                          sigW = input$sigW, 
                          phi = input$phi,
                          lb_goal = input$lb_pctMSY,
                          ub_goal = input$ub_pctMSY,
                          lb_manage = input$ub_pctMSY,
                          ub_manage = input$ub_pctMSY,
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
            yc4.5 = sum(SOC == "Yield", na.rm = TRUE) / N) %>%
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
                      "Yield concern"))


#Summary statistics for status quo and Rebuilding fisheries under reduced productivity.
#No evidence a rebuilding plan does much, likely because the Ricker is so flat with reduced productivity.
temp_sq <- sim_Seq_N1_df %>%
  mutate(scenario = "StatusQuo")
temp_rebuild <- sim_Seq_rebuild_df %>%
  mutate(scenario = "Rebuild")
sim_combined <- 
  rbind(temp_rebuild, temp_sq) %>%
  group_by(lnalpha, sigW, phi, scenario) %>%
  select(-F, -U, -N) %>% 
  mutate(change = case_when(SOC != lag(SOC) ~ TRUE, TRUE ~ FALSE),
         n_change = cumsum(change),
         deviation_lb = resid / lb_goal) %>%
  group_by(lnalpha, sigW, phi, scenario, n_change, SOC) %>%
  summarise(length = length(n_change),
            deviation_lb = mean(deviation_lb),
            miss = sum(resid < 0) / length,
            mean_S = mean(S))

# How many SOC designations
# maybe a weak trend towards fewer designation while in rebuilding status
delta_SOC <- list()
for(i in 1:length(unique(input$lnalpha))){
  delta_SOC[[i]] <-
    sim_combined %>% 
    group_by(lnalpha, sigW, phi, scenario, SOC) %>% 
    filter(SOC != "No concern") %>%
    summarize(delta_SOC = length(n_change)) %>% 
    ggplot(aes(x = scenario, y  = delta_SOC, fill = SOC)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        ncol = 3, nrow = 3, page = i)
}
#delta_SOC[[1]]
delta_SOC[[2]]
#delta_SOC[[3]]

# Percent of time spend in each designation
# I don't see a patterns relative to status qou or rebuild but there are patterns relative the SOC parameters.
pct_SOC <- list()
for(i in 1:length(unique(input$lnalpha))){
  pct_SOC[[i]] <-
    sim_combined %>% 
    group_by(lnalpha, sigW, phi, scenario, SOC) %>% 
    summarize(pct_SOC = sum(length)/input_sims) %>% 
    ggplot(aes(x = scenario, y  = pct_SOC, fill = SOC)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        ncol = 3, nrow = 3, page = i)
}
#pct_SOC[[1]]
pct_SOC[[2]]
#pct_SOC[[3]]

#Percent of lb_goal (the goal lb based on the historic dataset) met while in conservation or management concern (i.e. while managing to a rebuilding goal)
#Unsurprisingly we do archive larger escapements but it's not huge
S_soc <- list()
for(i in 1:length(unique(input$lnalpha))){
  S_soc[[i]] <-
    sim_combined %>% 
    group_by(lnalpha, sigW, phi, scenario, SOC) %>% 
    filter(SOC %in% c("Conservation", "Management")) %>%
    ggplot(aes(x = scenario, y  = deviation_lb)) + 
    geom_boxplot() +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        ncol = 3, nrow = 3, page = i)
}
S_soc[[1]] + coord_cartesian(ylim = c(-0.5, 0.5))
S_soc[[2]] + coord_cartesian(ylim = c(-0.5, 0.5))
S_soc[[3]] + coord_cartesian(ylim = c(-0.5, 0.5))

#Plot of point estimates for each Ricker parameter combination
tab_rebuild <- sim_Seq_rebuild_df %>%
  group_by(lnalpha, sigW, phi, lnalpha_red) %>% 
  summarise(S = mean(S), N = mean(N), R = mean(R)) %>%
  mutate(scenario = "Rebuild") %>%
  pivot_longer(cols = c(S, N, R), names_to = "stat", values_to = "value")
tab_statusquo <- sim_Seq_N1_df %>%
  group_by(lnalpha, sigW, phi, lnalpha_red) %>% 
  summarise(S = mean(S), N = mean(N), R = mean(R)) %>%
  mutate(scenario = "Status_Quo") %>%
  pivot_longer(cols = c(S, N, R), names_to = "stat", values_to = "value")
bind_rows(tab_rebuild, tab_statusquo) %>%
  pivot_wider(names_from = scenario, values_from = value) %>%
  ggplot(aes(x = Status_Quo, y = Rebuild, color = stat)) +
    geom_jitter(width = 50, height = 50) +
    geom_abline(slope = 1, intercept = 0)

# Length of SOC designations (Conservation and Management)
# I don't see much here
length_cc <- list()
for(i in 1:length(unique(input$lnalpha))){
  length_cc[[i]] <-
    sim_combined %>% 
    group_by(lnalpha, sigW, phi, scenario, SOC) %>% 
    filter(SOC %in% c("Conservation", "Management")) %>%
    ggplot(aes(x = scenario, y  = length)) + 
    geom_boxplot() +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        ncol = 3, nrow = 3, page = i)
}
length_cc[[1]]
length_cc[[2]] + coord_cartesian(ylim = c(0, 25))
length_cc[[3]]


# Rebuilding with shifting lnalpha ----------------------------------------
#Build a vector of lnalpha values where regimes shift every 15-25 years
get_lnalpha_ts <- function(base, red){
  out <- 
    mapply(rep, each = runif(1000, 15, 25), MoreArgs = list(x = c(base, red)), SIMPLIFY = FALSE) %>% 
    do.call(c, .)
  out[1:input_sims]
}

lnalpha_regime <- mapply(get_lnalpha_ts, base = input$lnalpha, red = input$lb_pctMSY * input$beta, SIMPLIFY = FALSE)

#fixed escapement goals
sim_Seq_statusquo_list <- mapply(FUN = simSR_goal, 
                                 lnalpha = lnalpha_regime,  
                                 sigW = input$sigW, 
                                 phi = input$phi,
                                 lb_goal = input$lb_pctMSY,
                                 ub_goal = input$ub_pctMSY,
                                 MoreArgs = list(beta = beta, 
                                                 age0 = Chinook_age,
                                                 Sims0 = input_sims, 
                                                 sigN = 0.1,
                                                 sigF = 0,
                                                 Hfun = H_goal),
                                 SIMPLIFY = FALSE)
sim_Seq_statusquo_df <-
  lapply(sim_Seq_statusquo_list, as.data.frame) %>%
  lapply(function(x) rename(x, lnalpha_red = lnalpha)) %>%
  lapply(function(x){mutate(x, sim = row_number())}) %>%
  do.call("rbind", .) %>%
  mutate(lnalpha = rep(input$lnalpha, each = input_sims)) %>%
  filter(R != 0) %>% #remove rows where population perished
  mutate(resid = S - lb_goal) %>%  
  select(-starts_with("N_age")) %>%
  group_by(lnalpha, sigW, phi) %>%
  mutate(resid_roll = roll_mean(x = resid, 5, align = "right", fill = NA) / lb_goal)

#raise escapement goals when in Conservation or Management concern
sim_Seq_rebuild_list <- mapply(FUN = simSR_goal, 
                               lnalpha = lnalpha_regime,  
                               sigW = input$sigW, 
                               phi = input$phi,
                               lb_goal = input$lb_pctMSY,
                               ub_goal = input$ub_pctMSY,
                               lb_manage = input$ub_pctMSY,
                               ub_manage = input$ub_pctMSY,
                               MoreArgs = list(beta = beta, 
                                               age0 = Chinook_age,
                                               Sims0 = input_sims, 
                                               sigN = 0.1,
                                               sigF = 0,
                                               Hfun = H_soc),
                               SIMPLIFY = FALSE)
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

#Average S, N and R under each regime
sim_Seq_rebuild_df %>% group_by(lnalpha, sigW, phi, lnalpha_red) %>% summarise(S = mean(S), N = mean(N), R = mean(R)) %>% print( n = 100)

temp_sq <- sim_Seq_statusquo_df %>%
  mutate(scenario = "StatusQuo")
temp_rebuild <- sim_Seq_rebuild_df %>%
  mutate(scenario = "Rebuild")
sq_rebuild_regime_combined <- 
  rbind(temp_rebuild, temp_sq) %>%
  group_by(lnalpha, sigW, phi, scenario) %>%
  select(-F, -U, -N) %>% 
  mutate(change = case_when(SOC != lag(SOC) ~ TRUE, TRUE ~ FALSE),
         n_change = cumsum(change),
         deviation_lb = resid / lb_goal) %>%
  group_by(lnalpha, sigW, phi, scenario, n_change, SOC) %>%
  summarise(length = length(n_change),
            deviation_lb = mean(deviation_lb),
            miss = sum(resid < 0) / length,
            mean_S = mean(S))

sq_rebuild_regime_NSOC <- list()
for(i in 1:length(unique(input$lnalpha))){
  sq_rebuild_regime_NSOC[[i]] <-
    sq_rebuild_regime_combined %>% 
    group_by(lnalpha, sigW, phi, scenario, SOC) %>% 
    summarize(delta_SOC = length(n_change)) %>% 
    ggplot(aes(x = scenario, y  = delta_SOC, fill = SOC)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        ncol = 3, nrow = 3, page = i)
}
#sq_rebuild_regime_NSOC[[1]]
sq_rebuild_regime_NSOC[[2]]
#sq_rebuild_regime_NSOC[[3]]

pct_SOC <- list()
for(i in 1:length(unique(input$lnalpha))){
  pct_SOC[[i]] <-
    sq_rebuild_regime_combined %>% 
    group_by(lnalpha, sigW, phi, scenario, SOC) %>% 
    summarize(pct_SOC = sum(length)/input_sims) %>% 
    ggplot(aes(x = scenario, y  = pct_SOC, fill = SOC)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        ncol = 3, nrow = 3, page = i)
}
#length_SOC[[1]]
pct_SOC[[2]]
#length_SOC[[3]]

S_cc <- list()
for(i in 1:length(unique(input$lnalpha))){
  S_cc[[i]] <-
    sq_rebuild_regime_combined %>% 
    group_by(lnalpha, sigW, phi, scenario, SOC) %>% 
    filter(SOC %in% c("Conservation", "Management")) %>%
    ggplot(aes(x = scenario, y  = deviation_lb)) + 
    geom_boxplot() +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        ncol = 3, nrow = 3, page = i)
}
S_cc[[1]] + coord_cartesian(ylim = c(-0.5, 0.5))
S_cc[[2]] + coord_cartesian(ylim = c(-0.5, 0.5))
S_cc[[3]] + coord_cartesian(ylim = c(-0.5, 0.5))

length_cc <- list()
for(i in 1:length(unique(input$lnalpha))){
  length_cc[[i]] <-
    sq_rebuild_regime_combined %>% 
    group_by(lnalpha, sigW, phi, scenario, SOC) %>% 
    filter(SOC %in% c("Conservation", "Management")) %>%
    ggplot(aes(x = scenario, y  = length)) + 
    geom_boxplot() +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        ncol = 3, nrow = 3, page = i)
}
length_cc[[1]]
length_cc[[2]] + coord_cartesian(ylim = c(0, 25))
length_cc[[3]]




# #Shnute params
# #Less parameter correlation during estimation
# #Can. J. Fish. Aquat. Sci. 53: 12811293 (1996).
# #h: harvest rate
# #C: MSY
# 
# #to get h and C from alpha and beta
# h <- rep(NA, 101)
# alpha <- rep(NA, 100)
# c <- rep(NA, 100)
# 
# for(i in 1:100){
#   h[1] <- 0.5
#   alpha[i] <- exp(h[i])/(1 - h[i])
#   h[i + 1] <- h[i] + (1 - h[i])/(2 - h[i])*log(exp(1.5)/alpha[i]) #alpha = 1.5
# }
# h[100] #h
# h[100]^2 / ( 1 - h[100]) / 0.0001 #C, beta = 0.0001
# 
# #Plot equivalence
# ricker <- function(s) s*exp(1.5 - 0.0001*s)
# ricker_alt <- function(s, h = 0.64, c = 8484) s/(1 - h) * exp(h - h^2/(1 - h)*s/c)
# s <- 1:30000
# plot(s, ricker(s), type = 'l') #alpha, beta Ricker
# abline(0, 1)
# Smsy <- 1.5 / 0.0001 * (0.5 - 0.07 * 1.5) #Smsy for the ricker above
# Rmsy <- ricker(Smsy)
# C <- Rmsy - Smsy
# h <- C/Rmsy
# lines(ricker_alt(s, h, C), col = "red") #h, C Ricker
# #Find a value of h that provides negligible yielde for the same beta
# f <- function(h, beta = 0.0001, C) {h^2 / ( 1 - h) / beta - C}
# uniroot(f, c(0, 1), C = 100)
# lines(ricker_alt(s, 0.095, 100), col = "green")
# 
# #Curve changes with constant C and varied h
# test <-
#   data.frame(C = 8757,
#            h = seq(.4, .9, length.out = 4)) %>%
#   mutate(a = exp(h) / (1 - h),
#          b = h^2 / ((1 - h)*C))
# test
# plot(s, ricker_alt(s, test$h[1], test$C[1]))
# abline(0, 1)
# lines(s, ricker_alt(s, test$h[2], test$C[2]))
# lines(s, ricker_alt(s, test$h[3], test$C[3]))
# lines(s, ricker_alt(s, test$h[4], test$C[4]))
# 
# #Curve changes with constant h and varied C
# test2 <-
#   data.frame(C = c(100, 3000, 8757, 12000),
#              h = .57) %>%
#   mutate(a = exp(h) / (1 - h),
#          b = h^2 / ((1 - h)*C))
# test2
# 
# plot(s, ricker_alt(s, test2$h[4], test2$C[4]), col = "red")
# abline(0, 1)
# lines(s, ricker_alt(s, test2$h[1], test2$C[1]), col = "red")
# lines(s, ricker_alt(s, test2$h[2], test2$C[2]), col = "red")
# lines(s, ricker_alt(s, test2$h[3], test2$C[3]), col = "red")
# 
# #Neither of the above patterns seem realistic
# #Curve changes with constant beta and varied C
# #So can describe what I think is the most likely situation
# #Would be interesting to see how it works in estimation.
# test3 <-
#   data.frame(C = c(100, 3000, 8757, 12000),
#              h = c(uniroot(f, c(0, 1), C = 100)[[1]],
#                    uniroot(f, c(0, 1), C = 3000)[[1]],
#                    uniroot(f, c(0, 1), C = 8757)[[1]],
#                    uniroot(f, c(0, 1), C = 12000)[[1]])) %>%
#   mutate(a = exp(h) / (1 - h),
#          b = h^2 / ((1 - h)*C))
# test3
# 
# plot(s, ricker_alt(s, test3$h[4], test3$C[4]), col = "green")
# abline(0, 1)
# lines(s, ricker_alt(s, test3$h[1], test3$C[1]), col = "green")
# lines(s, ricker_alt(s, test3$h[2], test3$C[2]), col = "green")
# lines(s, ricker_alt(s, test3$h[3], test3$C[3]), col = "green")
