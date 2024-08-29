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
                       power = c(0.6, 0.7, 0.8)) #how to pick these. High for SOC sims. Should be lower for OPY sims.
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
sim_base <- 
  mapply(FUN = sim_Ricker,
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
         SIMPLIFY = FALSE) %>%
  do.call("rbind", .)

# * Plot time series --------------------------------------------------------
#writing function here since it has dubious use outside of this single application
plot_ts <- function(dat, lnalpha0){ 
  df <- 
    dat %>%
    select(lnalpha, sigW, phi, sim, lb_goal, ub_goal, S, N) %>%
    filter(lnalpha == lnalpha0) %>%
    pivot_longer(c(S, N), names_to = "stat", values_to = "value")
  
  limit <- max(quantile(df$value[df$stat == "N"], .95, na.rm = TRUE), max(df$ub_goal) * 1.25)
  
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

#plot_ts(sim_base, 1)
plot_ts(sim_base, 1.5)
#plot_ts(sim_base, 2)

# * Table of concern criteria -----------------------------------------------
sim_base %>%
  group_by(lnalpha, sigW, phi) %>%
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

#Need to redo. residual not calculated yet.
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
# % of average productivity ~ 25%
dat_pctred <- 
  input %>% 
  mutate(lnalpha_p = get_lnalpha_p(lnalpha, beta, sigW, phi),
         lnalpha_p_red = lb_pctMSY  * beta,
         pct_red = lnalpha_p_red / lnalpha_p)
mean(dat_pctred$pct_red)
median(dat_pctred$pct_red)
range(dat_pctred$pct_red)


# % of average productivity lower bound ~ .47%
input %>% 
  mutate(lnalpha_p_red = lb_pctMSY  * beta,
         lb_pctMSY_red = get_Smsy(lnalpha_p_red, beta, FALSE),
         pct_red = lb_pctMSY_red / lb_pctMSY) %>%
  print(n = 100)

#simulation results
# Simulation results
sim_Seq <- 
  mapply(FUN = sim_Ricker,
         lnalpha = input$lb_pctMSY*beta,  #reduced ln_alpha
         sigW = input$sigW,
         phi = input$phi,
         lb_goal = input$lb_pctMSY,
         ub_goal = input$ub_pctMSY,
         MoreArgs = list(beta = beta,
                         age0 = Chinook_age,
                         Sims0 = input_sims,
                         sigN = 0.2,
                         sigF = 0,
                         Hfun = H_goal),
         SIMPLIFY = FALSE) %>%
  lapply(function(x) rename(x, lnalpha_red = lnalpha)) %>%
  do.call("rbind", .) %>%
  left_join(input[, c("lnalpha", "lb_pctMSY", "ub_pctMSY")], by = c("lb_goal" = "lb_pctMSY", "ub_goal" = "ub_pctMSY"))

example_SOC <- sim_Seq[15:27, c(7, 9, 11:12, 19, 23:26)]
saveRDS(example_SOC, file = ".\\sims\\example_SOC.rds")

# * Plot time series --------------------------------------------------------
# because of the colors used teal below the lower bound represents times we fished below the lower bound w abundant fish... call the a yield concern
# while grey below the lower bound represent when the fish were not there to make the goal ... call that a management concern
# When N is far below the lb it is a conservation concern
plot_ts(sim_Seq, 1)
plot_ts(sim_Seq, 1.5)
plot_ts(sim_Seq, 2)



# * Table of concern criteria -----------------------------------------------
# estimation error for N introduces lb misses based on fishing and increased ub misses.
sim_Seq %>%
  group_by(lnalpha, sigW, phi) %>%
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

# * Compare to base -----------------------------------------------
temp_base <- sim_base %>%
  mutate(scenario = "base")
temp_Seq <- sim_Seq %>%
  mutate(scenario = "Seq") %>%
  select(-lnalpha_red)
baseSeq_combined <- 
  rbind(temp_base, temp_Seq) %>%
  group_by(lnalpha, sigW, phi, scenario) %>%
  select(-F, -U, -N) %>% 
  mutate(change = case_when(SOC != lag(SOC) ~ TRUE, TRUE ~ FALSE),
         n_change = cumsum(change),
         deviation_lb = (S - lb_goal) / lb_goal) %>%
  group_by(lnalpha, sigW, phi, scenario, n_change, SOC) %>%
  summarise(length = length(n_change),
            deviation_lb = mean(deviation_lb),
            miss = sum((S < lb_goal)) / length,
            mean_S = mean(S))

baseSeq_pct_SOC <- list()
for(i in 1:length(unique(input$lnalpha))){
  baseSeq_pct_SOC[[i]] <-
    baseSeq_combined %>% 
    group_by(lnalpha, sigW, phi, scenario, SOC) %>% 
    summarize(pct_SOC = sum(length)/input_sims) %>% 
    ggplot(aes(x = scenario, y  = pct_SOC, fill = SOC)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        ncol = 3, nrow = 3, page = i)
}
baseSeq_pct_SOC[[1]]
baseSeq_pct_SOC[[2]]
baseSeq_pct_SOC[[3]]

# * Average 5 year rolling residual as a diagnostic -------------------------
#Many concerns.
# plot_resid(sim_Seq_df, 1)
# plot_resid(sim_Seq_df, 1.5)
# plot_resid(sim_Seq_df, 2)


# Rebuild -----------------------------------------------------------------
# Can managing to a higher goal during management or conservation concerns reduce the duration of concern designations or 
# increase escapement during concern designations
sim_Seq_rebuild <- 
  mapply(FUN = sim_Ricker, 
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
                         sigN = 0.2,
                         sigF = 0,
                         Hfun = H_soc),
         SIMPLIFY = FALSE) %>%
  lapply(function(x) rename(x, lnalpha_red = lnalpha)) %>%
  do.call("rbind", .) %>%
  left_join(input[, c("lnalpha", "lb_pctMSY", "ub_pctMSY")], by = c("lb_goal" = "lb_pctMSY", "ub_goal" = "ub_pctMSY"))



# * Plot time series --------------------------------------------------------
# because of the colors used teal below the lower bound represents times we fished below the lower bound w abundant fish... call the a management concern
# teal above the upper bound represents underutilized yield... call that a yield concern
# while grey below the lower bound represent when the fish were not there to make the goal ... call that a conservation concern
# Note that some phi and sigma W combinations are not sustainable
plot_ts(sim_Seq_rebuild, 1)
plot_ts(sim_Seq_rebuild, 1.5)
plot_ts(sim_Seq_rebuild, 1.5) + coord_cartesian(xlim = c(600, 625))
plot_ts(sim_Seq_rebuild, 2)


# * Table of concern criteria -----------------------------------------------
sim_Seq_rebuild %>%
  group_by(lnalpha, sigW, phi) %>%
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

# * Compare to Seq  -----------------------------------------------
# * * pct of time in each SOC designation  -----------------------------------------------
temp_Seq <- sim_Seq %>%
  mutate(scenario = "Seq") %>%
  select(-lnalpha_red)
temp_Seq_rebuild <- sim_Seq_rebuild %>%
  mutate(scenario = "Seq_rebuild") %>%
  select(-lnalpha_red)
Seq_Seqrebuild_combined <- 
  rbind(temp_Seq, temp_Seq_rebuild) %>%
  group_by(lnalpha, sigW, phi, scenario) %>%
  select(-F, -U, -N) %>% 
  mutate(change = case_when(SOC != lag(SOC) ~ TRUE, TRUE ~ FALSE),
         n_change = cumsum(change),
         deviation_lb = (S - lb_goal) / lb_goal) %>%
  group_by(lnalpha, sigW, phi, scenario, n_change, SOC) %>%
  summarise(length = length(n_change),
            deviation_lb = mean(deviation_lb),
            miss = sum((S < lb_goal)) / length,
            mean_S = mean(S))

SeqSeqrebuild_pct_SOC <- list()
for(i in 1:length(unique(input$lnalpha))){
  SeqSeqrebuild_pct_SOC[[i]] <-
    Seq_Seqrebuild_combined %>% 
    group_by(lnalpha, sigW, phi, scenario, SOC) %>% 
    summarize(pct_SOC = sum(length)/input_sims) %>% 
    ggplot(aes(x = scenario, y  = pct_SOC, fill = SOC)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_x_discrete(name = "Management",
                     breaks = c("Seq", "Seq_rebuild"), 
                     labels = c("Standard", "Rebuilding")) +
    scale_y_continuous(name = "Percent of Simulations") +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        ncol = 3, nrow = 3, page = i)
}
SeqSeqrebuild_pct_SOC[[1]]
SeqSeqrebuild_pct_SOC[[2]]
SeqSeqrebuild_pct_SOC[[3]]

# * * mean abundance ------------------------------------------------------
#Summary statistics for status quo and Rebuilding fisheries under reduced productivity
#Plot of point estimates for each Ricker parameter combination
tab_rebuild <- sim_Seq_rebuild %>%
  filter(max(sim) == input_sims) %>%
  filter(SOC %in% c("Conservation", "Management")) %>%
  group_by(lnalpha, sigW, phi) %>% 
  summarise(S = mean(S), N = mean(N), R = mean(R)) %>%
  mutate(scenario = "Rebuild") %>%
  pivot_longer(cols = c(S, N, R), names_to = "stat", values_to = "value")
tab_statusquo <- sim_Seq %>%
  filter(max(sim) == input_sims) %>%
  filter(SOC %in% c("Conservation", "Management")) %>%
  group_by(lnalpha, sigW, phi) %>% 
  summarise(S = mean(S), N = mean(N), R = mean(R)) %>%
  mutate(scenario = "Status_Quo") %>%
  pivot_longer(cols = c(S, N, R), names_to = "stat", values_to = "value")
bind_rows(tab_rebuild, tab_statusquo) %>%
  pivot_wider(names_from = scenario, values_from = value) %>%
  ggplot(aes(x = Status_Quo, y = Rebuild, color = stat)) +
  geom_jitter(width = 50, height = 50) +
  geom_abline(slope = 1, intercept = 0)

#Compare different SOC designation statistics
  # mutate(change = case_when(SOC != lag(SOC) ~ TRUE, TRUE ~ FALSE),
  #        n_change = cumsum(change),
  #        deviation_lb = (S - lb_goal) / lb_goal) %>%
  # group_by(lnalpha, sigW, phi, scenario, n_change, SOC) %>%
  # summarise(length = length(n_change),
  #           deviation_lb = mean(deviation_lb),
  #           miss = sum(S < lb_goal) / length,
  #           mean_S = mean(S))

# Regimes - shifting lnalpha ----------------------------------------
#Build a vector of lnalpha values where regimes shift every 15-25 years
get_lnalpha_ts <- function(base, red){
  out <- 
    mapply(rep, each = runif(1000, 15, 25), MoreArgs = list(x = c(base, red)), SIMPLIFY = FALSE) %>% 
    do.call(c, .)
  out[1:input_sims]
}

lnalpha_regime <- mapply(get_lnalpha_ts, base = input$lnalpha, red = input$lb_pctMSY * input$beta, SIMPLIFY = FALSE)

#fixed escapement goals
sim_regime_statusquo <- 
  mapply(FUN = sim_Ricker, 
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
         SIMPLIFY = FALSE) %>%
  lapply(function(x) rename(x, lnalpha_red = lnalpha)) %>%
  do.call("rbind", .) %>%
  left_join(input[, c("lnalpha", "lb_pctMSY", "ub_pctMSY")], by = c("lb_goal" = "lb_pctMSY", "ub_goal" = "ub_pctMSY"))

#raise escapement goals when in Conservation or Management concern
sim_regime_rebuild <- 
  mapply(FUN = sim_Ricker,
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
         SIMPLIFY = FALSE)  %>%
  lapply(function(x) rename(x, lnalpha_red = lnalpha)) %>%
  do.call("rbind", .) %>%
  left_join(input[, c("lnalpha", "lb_pctMSY", "ub_pctMSY")], by = c("lb_goal" = "lb_pctMSY", "ub_goal" = "ub_pctMSY"))


# * Compare to strategies  -----------------------------------------------
# * * pct of time in each SOC designation  -----------------------------------------------
temp_regime <- sim_regime_statusquo %>%
  mutate(scenario = "regime_statusquo") %>%
  select(-lnalpha_red)
temp_regime_rebuild <- sim_regime_rebuild %>%
  mutate(scenario = "regime_rebuild") %>%
  select(-lnalpha_red)
regime_combined <- 
  rbind(temp_regime, temp_regime_rebuild) %>%
  group_by(lnalpha, sigW, phi, scenario) %>%
  select(-F, -U, -N) %>% 
  mutate(change = case_when(SOC != lag(SOC) ~ TRUE, TRUE ~ FALSE),
         n_change = cumsum(change),
         deviation_lb = (S - lb_goal) / lb_goal) %>%
  group_by(lnalpha, sigW, phi, scenario, n_change, SOC) %>%
  summarise(length = length(n_change),
            deviation_lb = mean(deviation_lb),
            miss = sum((S < lb_goal)) / length,
            mean_S = mean(S))

regime_pct_SOC <- list()
for(i in 1:length(unique(input$lnalpha))){
  regime_pct_SOC[[i]] <-
    regime_combined %>% 
    group_by(lnalpha, sigW, phi, scenario, SOC) %>% 
    summarize(pct_SOC = sum(length)/input_sims) %>% 
    ggplot(aes(x = scenario, y  = pct_SOC, fill = SOC)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        ncol = 3, nrow = 3, page = i)
}
regime_pct_SOC[[1]]
regime_pct_SOC[[2]]
regime_pct_SOC[[3]]

# * * mean abundance ------------------------------------------------------
#Summary statistics for status quo and Rebuilding fisheries under reduced productivity.

#Plot of point estimates for each Ricker parameter combination
tab_regime_rebuild <- sim_regime_rebuild %>%
  filter(max(sim) == input_sims) %>%
  filter(SOC %in% c("Conservation", "Management")) %>%
  group_by(lnalpha, sigW, phi) %>% 
  summarise(S = mean(S), N = mean(N), R = mean(R)) %>%
  mutate(scenario = "Rebuild") %>%
  pivot_longer(cols = c(S, N, R), names_to = "stat", values_to = "value")
tab_regime_statusquo <- sim_regime_statusquo %>%
  filter(max(sim) == input_sims) %>%
  filter(SOC %in% c("Conservation", "Management")) %>%
  group_by(lnalpha, sigW, phi) %>% 
  summarise(S = mean(S), N = mean(N), R = mean(R)) %>%
  mutate(scenario = "Status_Quo") %>%
  pivot_longer(cols = c(S, N, R), names_to = "stat", values_to = "value")
bind_rows(tab_regime_rebuild, tab_regime_statusquo) %>%
  pivot_wider(names_from = scenario, values_from = value) %>%
  ggplot(aes(x = Status_Quo, y = Rebuild, color = stat)) +
  geom_jitter(width = 0, height = 0) +
  geom_abline(slope = 1, intercept = 0)


# Depensation ----------------------------------------
#example plot
#parameter set
temp_par <- gamma_par(1 / beta, Ricker(1.5, 0.0001, 1 / beta), c(1, 1.1, 1.6))

#depensation points
dat_dep <- 
  as.data.frame(temp_par) %>%
  filter(c != 1) %>%
  rowwise() %>%
  mutate(S = uniroot(function(S) {-a * S^(c-2) * (b * S - c + 1)*exp(-b * S)}, c(1, 10000), tol = 0.00001)$root,
         R = SRgamma(a, b, c, S = S),
         logRS = log(SRgamma(a, b, c, S = S) / S),
         gamma = as.character(c)) %>%
  ungroup() %>%
  pivot_longer(c(R, logRS), names_to = "Parameter", values_to = "Value")

#critical values
dat_crit <- 
  as.data.frame(temp_par) %>%
  filter(c > 1.2) %>%
  rowwise() %>%
  mutate(S = uniroot(function(S) {SRgamma(a, b, c, S) - S}, c(5, 2500), tol = 0.00001)$root,
         R = SRgamma(a, b, c, S = S),
         logRS = log(SRgamma(a, b, c, S = S) / S),
         gamma = as.character(c)) %>%
  ungroup() %>%
  pivot_longer(c(R, logRS), names_to = "Parameter", values_to = "Value")

#create plot
data.frame(S = rep(seq(0, 20000, length.out = 100), times = length(temp_par[[3]])),
           alpha = rep(temp_par[[1]], each = 100),
           beta = rep(temp_par[[2]], each = 100),
           gamma = rep(temp_par[[3]], each = 100)) %>%
  mutate(R = SRgamma(alpha =  alpha, beta = beta, gamma = gamma, S),
         logRS = log(R / S)) %>%
  pivot_longer(c(R, logRS), names_to = "Parameter", values_to = "Value") %>%
  ggplot(aes(x = S, y = Value, color = as.character(gamma))) +
  geom_line() +
  geom_point(data = dat_dep) +
  geom_segment(data = dat_dep, aes(xend = S, y = -Inf, yend = Value), linetype = 4) +
  geom_segment(data = dat_crit, aes(xend = S, y = -Inf, yend = Value), linetype = 2) +  
  geom_abline(aes(slope = x, intercept = y), data.frame(x = 1, y = 0, Parameter = "R")) +
  geom_hline(aes(yintercept = y), data.frame(y = 0, Parameter = "logRS")) +
  scale_x_continuous(breaks = seq(0, 20000, 2500)) +
  scale_color_discrete(name = "Gamma", 
                       labels = c("1 (Ricker)", "1.1", "1.6")) +
  theme_bw() +
  facet_grid(Parameter ~ ., scales = "free_y")

#Create input parameters
input <- 
  input %>%
  rowwise() %>%
  mutate(Rmax = Ricker(lnalpha = get_lnalpha_p(lnalpha, beta, sigW, phi), beta = beta, S = 1 / beta),
         Rmax_Seq = Ricker(lnalpha = lb_pctMSY * beta, beta = beta, S = 1 / beta),
         Seq = get_lnalpha_p(lnalpha, beta, sigW, phi) / beta,
         c = 1.1,
         a = gamma_par(1 / beta, Rmax, c)[[1]],
         b = gamma_par(1 / beta, Rmax, c)[[2]],
         a_Seq = gamma_par(1 / beta, Rmax_Seq, c)[[1]],
         b_Seq = gamma_par(1 / beta, Rmax_Seq, c)[[2]]) %>%
  rowwise() %>%
  #mutate(b_match = uniroot(function(b) {a * Seq^c * exp(-b * Seq) - Seq}, interval = c(b/1.1, b/.9), tol = 0.00001)$root,
         #crit = uniroot(function(S) a * S^c * exp(-b_match * S) - S, c(5, 1000), tol = 0.00001)$root,
         #dep = uniroot(function(S) {-a * S^(c-2) * (b_match * S - c + 1)*exp(-b_match * S)}, c(0, 5000), tol = 0.00001)$root,
         #lnRS_gamma_crit = log(SR_gamma(a = a, b = b_match, c = c, S = crit) / crit),
         #lnRS_gamma_dep = log(SR_gamma(a = a, b = b_match, c = c, S = dep) / dep)) %>%
  ungroup()


gamma_plots <- list()
for(i in 1:length(unique(input$lnalpha))){
  gamma_plots[[i]] <- 
    input %>%
    slice(rep(1:n(), each = 251)) %>% 
    mutate(S = rep(seq(0, 5 * 1 / max(beta), by = 200), times = nrow(input)),
           R_Ricker = Ricker(lnalpha = get_lnalpha_p(lnalpha, beta, sigW, phi), beta = beta, S = S),
           RSeq_Ricker = Ricker(lnalpha = lb_pctMSY * beta, beta = beta, S = S),
           R_gamma = SRgamma(alpha = a, beta = b, gamma = c, S = S),
           RSeq_gamma = SRgamma(alpha = a_Seq, beta = b_Seq, gamma = c, S = S)) %>%
    ggplot(aes(x = S, y = R_Ricker)) +
    geom_rect(aes(xmin = lb_pctMSY, xmax = ub_pctMSY, ymin = -Inf, ymax = Inf), fill = "grey95") +
    geom_line() +
    geom_line(aes(y = RSeq_Ricker), color = "black", linetype = 2) +
    geom_line(aes(y = R_gamma), color = scales::hue_pal()(3)[[3]]) +
    geom_line(aes(y = RSeq_gamma), color = scales::hue_pal()(3)[[3]], linetype = 2) +
    geom_abline() +
    scale_y_continuous(name = "R") +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        #scales = "free_y",
                        ncol = 3, nrow = 3, page = i)
} 

#gamma_plots[[1]] + coord_cartesian(xlim = c(0, 20000), ylim = c(0, 20000))
gamma_plots[[2]] + coord_cartesian(xlim = c(0, 30000), ylim = c(0, 30000))
#gamma_plots[[3]] + coord_cartesian(xlim = c(0, 50000), ylim = c(0, 50000))




lnRS_plots <- list()
for(i in 1:length(unique(input$lnalpha))){
  lnRS_plots[[i]] <- 
    input %>%
    slice(rep(1:n(), each = 201)) %>% 
    mutate(S = rep(seq(0, 6000, by = 6000*1/200), times = nrow(input)),
           R_Ricker = Ricker(lnalpha = get_lnalpha_p(lnalpha, beta, sigW, phi), beta = beta, S = S),
           RSeq_Ricker = Ricker(lnalpha = lb_pctMSY * beta, beta = beta, S = S),
           lnRS_Ricker = log(R_Ricker / S),
           lnRSSeq_Ricker = log(RSeq_Ricker / S),
           R_gamma = SRgamma(alpha = a, beta = b, gamma = c, S = S),
           RSeq_gamma = SRgamma(alpha = a_Seq, beta = b_Seq, gamma = c, S = S),
           lnRS_gamma = log(R_gamma / S),
           lnRSSeq_gamma = log(RSeq_gamma / S)) %>%
    rowwise() %>%
    mutate(dep = uniroot(function(S) {-a * S^(c-2) * (b * S - c + 1)*exp(-b * S)}, c(0, 5000), tol = 0.00001)$root,
           lnRS_gamma_dep = log(SRgamma(a, b, c, S = dep) / dep),
           dep_Seq = uniroot(function(S) {-a_Seq * S^(c-2) * (b_Seq * S - c + 1)*exp(-b_Seq * S)}, c(0, 5000), tol = 0.00001)$root,
           lnRS_gamma_dep_Seq = log(SRgamma(a_Seq, b_Seq, c, S = dep_Seq) / dep_Seq)) %>%
    ungroup() %>%
    ggplot(aes(x = S, y = lnRS_Ricker)) +
    geom_rect(aes(xmin = lb_pctMSY, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "grey95") +
    geom_line() +
    geom_line(aes(y = lnRSSeq_Ricker), color = "black", linetype = 2) +
    geom_line(aes(y = lnRS_gamma), color = scales::hue_pal()(3)[[3]]) +
    geom_line(aes(y = lnRSSeq_gamma), color = scales::hue_pal()(3)[[3]], linetype = 2) +
    geom_point(aes(x = dep, y = lnRS_gamma_dep), color = scales::hue_pal()(3)[[3]]) +
    geom_point(aes(x = dep_Seq, y = lnRS_gamma_dep_Seq), color = scales::hue_pal()(3)[[3]]) +
    geom_hline(yintercept = 0) +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        #scales = "free_y",
                        ncol = 3, nrow = 3, page = i)
} 

lnRS_plots[[1]]
gamma_plots[[1]] + coord_cartesian(xlim = c(0, 500), ylim = c( 0, 500))
lnRS_plots[[2]]
lnRS_plots[[3]]

# * Historic Productivity ---------------------------------------------------
sim_gamma <- 
  mapply(FUN = sim_SRgamma, 
         alpha = input$a,
         beta = input$b,
         gamma = input$c,
         sigW = input$sigW, 
         phi = input$phi,
         lb_goal = input$lb_pctMSY,
         ub_goal = input$ub_pctMSY,
         MoreArgs = list(age0 = Chinook_age,
                         Sims0 = input_sims, 
                         sigN = 0.2,
                         sigF = 0,
                         Hfun = H_goal),
         SIMPLIFY = FALSE) %>%
  do.call("rbind", .) %>%
  left_join(input[, c("lnalpha", "lb_pctMSY", "ub_pctMSY")], by = c("lb_goal" = "lb_pctMSY", "ub_goal" = "ub_pctMSY"))

#plot_ts(sim_base, 1)
plot_ts(sim_gamma, 1.5) +coord_cartesian(ylim = c( 0, 30000))
#plot_ts(sim_base, 2)

# * * Compare to sim_base from Ricker  -----------------------------------------------
# * * * pct of time in each SOC designation  -----------------------------------------------
temp_gamma <- sim_gamma %>%
  select(-alpha, -beta, -gamma) %>%
  mutate(scenario = "gamma")
temp_base2 <- temp_base %>%
  select(-lnalpha.y)
basegamma_combined <- 
  rbind(temp_base2, temp_gamma) %>%
  group_by(lnalpha, sigW, phi, scenario) %>%
  select(-F, -U, -N) %>% 
  mutate(change = case_when(SOC != lag(SOC) ~ TRUE, TRUE ~ FALSE),
         n_change = cumsum(change),
         deviation_lb = (S - lb_goal) / lb_goal) %>%
  group_by(lnalpha, sigW, phi, scenario, n_change, SOC) %>%
  summarise(length = length(n_change),
            deviation_lb = mean(deviation_lb),
            miss = sum((S < lb_goal)) / length,
            mean_S = mean(S))

basegamma_pct_SOC <- list()
for(i in 1:length(unique(input$lnalpha))){
  basegamma_pct_SOC[[i]] <-
    basegamma_combined %>% 
    group_by(lnalpha, sigW, phi, scenario, SOC) %>% 
    summarize(pct_SOC = sum(length)/input_sims) %>% 
    ggplot(aes(x = scenario, y  = pct_SOC, fill = SOC)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        ncol = 3, nrow = 3, page = i)
}
basegamma_pct_SOC[[1]]
basegamma_pct_SOC[[2]]
basegamma_pct_SOC[[3]]

# no abundance comparisons because they focus on abundance during 
# management and conservation concerns.

# * Low Productivity --------------------------------------------------------
sim_Seq_gamma <- 
  mapply(FUN = sim_SRgamma, 
         alpha = input$a_Seq,
         beta = input$b_Seq,
         gamma = input$c,
         sigW = input$sigW, 
         phi = input$phi,
         lb_goal = input$lb_pctMSY,
         ub_goal = input$ub_pctMSY,
         MoreArgs = list(age0 = Chinook_age,
                         Sims0 = input_sims, 
                         sigN = 0.2,
                         sigF = 0,
                         Hfun = H_goal),
         SIMPLIFY = FALSE) %>%
  do.call("rbind", .) %>%
  left_join(input[, c("lnalpha", "lb_pctMSY", "ub_pctMSY")], by = c("lb_goal" = "lb_pctMSY", "ub_goal" = "ub_pctMSY"))

#plot_ts(sim_Seq_gamma, 1)
plot_ts(sim_Seq_gamma, 1.5)
#plot_ts(sim_Seq_gamma, 2)


# * * Compare to sim_Seq from Ricker  -----------------------------------------------
# * * * pct of time in each SOC designation  -----------------------------------------------
temp_Seq_gamma <- sim_Seq_gamma %>%
  select(-alpha, -beta, -gamma) %>%
  mutate(scenario = "Seq_gamma")
temp_Seq2 <- temp_Seq %>%
  select(-lnalpha.y)
Seq_combined <- 
  rbind(temp_Seq2, temp_Seq_gamma) %>%
  group_by(lnalpha, sigW, phi, scenario) %>%
  select(-F, -U, -N) %>% 
  mutate(change = case_when(SOC != lag(SOC) ~ TRUE, TRUE ~ FALSE),
         n_change = cumsum(change),
         deviation_lb = (S - lb_goal) / lb_goal) %>%
  group_by(lnalpha, sigW, phi, scenario, n_change, SOC) %>%
  summarise(length = length(n_change),
            deviation_lb = mean(deviation_lb),
            miss = sum((S < lb_goal)) / length,
            mean_S = mean(S))

Seq_pct_SOC <- list()
for(i in 1:length(unique(input$lnalpha))){
  Seq_pct_SOC[[i]] <-
    Seq_combined %>% 
    group_by(lnalpha, sigW, phi, scenario, SOC) %>% 
    summarize(pct_SOC = sum(length)/input_sims) %>% 
    ggplot(aes(x = scenario, y  = pct_SOC, fill = SOC)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    facet_grid_paginate(paste0("\u03C6: ", phi) ~ paste0("ln(\u03B1): ", lnalpha) + paste0("\u03C3: ", sigW), 
                        ncol = 3, nrow = 3, page = i)
}
Seq_pct_SOC[[1]]
Seq_pct_SOC[[2]]
Seq_pct_SOC[[3]]

# * * * mean abundance ------------------------------------------------------
#Summary statistics for status quo and Rebuilding fisheries under reduced productivity.

#Plot of point estimates for each Ricker parameter combination
tab_Seq <- sim_Seq %>%
  filter(max(sim) == input_sims) %>%
  filter(SOC %in% c("Conservation", "Management")) %>%
  group_by(lnalpha, sigW, phi) %>% 
  summarise(S = mean(S), N = mean(N), R = mean(R)) %>%
  mutate(scenario = "Seq") %>%
  pivot_longer(cols = c(S, N, R), names_to = "stat", values_to = "value")
tab_Seq_gamma <- sim_Seq_gamma %>%
  filter(max(sim) == input_sims) %>%
  filter(SOC %in% c("Conservation", "Management")) %>%
  group_by(lnalpha, sigW, phi) %>% 
  summarise(S = mean(S), N = mean(N), R = mean(R)) %>%
  mutate(scenario = "Seq_gamma") %>%
  pivot_longer(cols = c(S, N, R), names_to = "stat", values_to = "value")
bind_rows(tab_Seq, tab_Seq_gamma) %>%
  pivot_wider(names_from = scenario, values_from = value) %>%
  ggplot(aes(x = Seq, y = Seq_gamma, color = stat)) +
  geom_jitter(width = 0, height = 0) +
  geom_abline(slope = 1, intercept = 0)


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
