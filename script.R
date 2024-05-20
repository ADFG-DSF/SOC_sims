# Simulate SOC listing criteria under a variety of productivity regimes

# Author: Adam Reimer
# Version: 2023-09-15

# Packages
packs <- c("tidyverse", "RcppRoll", "knitr", "ggforce")
lapply(packs, require, character.only = TRUE)

# source functions
function_files <- list.files(path=".\\functions")
lapply(function_files, function(x) source(paste0(".\\functions\\", x)))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# * Inputs -----------------------------------------------------------
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
  mutate(beta = beta, lna_p = lnalpha + (sigW * sigW / 2 / (1 - phi * phi)),
         lb_p = lna_p/beta*(0.5 - 0.07 * lna_p) * 0.8,
         ub_p = lna_p/beta*(0.5 - 0.07 * lna_p) * 1.6) %>%
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

# Simulation base case -------------------------------------------------------------------
#First simulation sigN = 0
#Not very realistic
# Simulation results
# sim_base_list <- mapply(FUN = simSR_goal,
#                         lnalpha = input$lnalpha,
#                         sigW = input$sigW,
#                         phi = input$phi,
#                         lb_sim = input$lb_p,
#                         ub_sim = input$ub_p,
#                         MoreArgs = list(beta = beta,
#                                         age0 = Chinook_age,
#                                         Sims = input_sims,
#                                         sigN = 0.2,
#                                         sigF = 0,
#                                         Hfun = H_goal),
#                         SIMPLIFY = FALSE)
# saveRDS(sim_base_list, file = ".\\sim_base_list.R")
sim_base_list <- readRDS(file = ".\\sim_base_list.R")

# Create dataframe --------------------------------------------------------
sim_base_df <-
  lapply(sim_base_list, as.data.frame) %>%
  lapply(function(x){mutate(x, sim = row_number())}) %>%
  do.call("rbind", .) %>%
  filter(R != 0) %>% #remove rows where population perished
  mutate(resid = S - lb) %>%  
  select(-starts_with("N_age")) %>%
  group_by(lnalpha, sigW, phi) %>%
  mutate(resid_roll = roll_mean(x = resid, 5, align = "right", fill = NA) / lb)

# * Plot time series --------------------------------------------------------
#writing function here since it has dubious use outside of this single application
plot_ts <- function(dat, lnalpha0){ 
  df <- 
    dat %>%
    select(lnalpha, sigW, phi, sim, lb, ub, S, N) %>%
    filter(lnalpha == lnalpha0) %>%
    pivot_longer(c(S, N), names_to = "stat", values_to = "value")
  
  limit <- quantile(df$value[df$stat == "N"], .95, na.rm = TRUE)
  
      ggplot(df, aes(x = sim, y = value, color = stat)) +
        geom_line(alpha = 0.4) +
        geom_hline(aes(yintercept = lb), linetype = 2) +
        geom_hline(aes(yintercept = ub), linetype = 2) +
        scale_y_continuous(limits = c(0, limit)) +
        scale_x_continuous(name = "Simulation year") +
        facet_grid(paste0("\u03C6: ", phi) ~ paste0("\u03C3: ", sigW)) +
        ggtitle(label = bquote(log(alpha): .(lnalpha0)))
}
# because of the colors used teal below the lower bound represents times we fished below the lower bound w abundant fish... call the a management concern
# teal above the upper bound represents underutilized yield... call that a yield concern
# while grey below the lower bound represent when the fish were not there to make the goal ... call that a conservation concern
# Note that some phi and sigma W combinations are not sustainable

plot_ts(sim_base_df, 1)
plot_ts(sim_base_df, 1.5)
plot_ts(sim_base_df, 2)



# Table of concern criteria -----------------------------------------------
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
                      "P(N < lower bound)",
                      "P(S < lb, N > lb)",
                      "P(S > upper bound)",
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


# Average 5 year rolling residual as a diagnostic -------------------------
plot_resid <- function(dat, lnalpha0){
  dat %>%
    filter(lnalpha == lnalpha0, !is.na(SOC)) %>%
    mutate(resid_pct = resid / lb) %>%
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


# # Reduced alpha: lb = Seq  -------------------------------------
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

# Simulation params -------------------------------------------------------
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
Ricker_plots_Seq

# * static reduced alpha -------------------------------------------------------------
static_Seq <- simSR_goal(lb_p*beta, beta, 0, 0,
                         age0 = c("1" = 1),
                         Sims0 = 1,
                         Hfun = H_target,
                         target = lb_p,
                         lb_sim = lb_p, ub_sim = ub_p)
static_Seq$U

#some harvest when random variability is included
sim_Seq <- simSR_goal(lb*beta, beta, 0.5, 0,
                      age0 = c("1" = 1),
                      Sims0 = 100,
                      Hfun = H_target,
                      target = lb_p,
                      lb_sim = lb_p, ub_sim = ub_p)
hist(sim_Seq$U)
mean(sim_Seq$U == 0) # mostly closed fisheries
mean(sim_Seq$U == 0 & sim_Seq$S < lb) # Often miss lower bound wo fishing
mean(sim_Seq$U > 0)
mean(sim_Seq$N_age * sim_Seq$U) # average yield under reduced productivity
mean(sim_Seq$N_age * sim_Seq$U) / (Ricker(lnalpha, beta, lb_p) - lb_p) # ~ 8% of MSY under average productivity conditions


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
#                           lb_sim = input$lb_p, #same bounds
#                           ub_sim = input$ub_p,
#                           MoreArgs = list(beta = beta,
#                                           age0 = Chinook_age,
#                                           Sims = input_sims,
#                                           sigN = 0.1,
#                                           sigF = 0,
#                                            Hfun = H_goal),
#                           SIMPLIFY = FALSE)
# saveRDS(sim_Seq_N1_list, file = ".\\sim_Seq_N1_list.R")
sim_Seq_N1_list <- readRDS(file = ".\\sim_Seq_N1_list.R")


# Create dataframe --------------------------------------------------------
sim_Seq_N1_df <-
  lapply(sim_Seq_N1_list, as.data.frame) %>%
  lapply(function(x) rename(x, lnalpha_red = lnalpha)) %>%
  lapply(function(x){mutate(x, sim = row_number())}) %>%
  do.call("rbind", .) %>%
  mutate(lnalpha = rep(input$lnalpha, each = input_sims)) %>%
  filter(R != 0) %>% #remove rows where population perished
  mutate(resid = S - lb) %>%  
  select(-starts_with("N_age")) %>%
  group_by(lnalpha, sigW, phi) %>%
  mutate(resid_roll = roll_mean(x = resid, 5, align = "right", fill = NA) / lb)



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
                      "P(N < lower bound)",
                      "P(S < lb, N > lb)",
                      "P(S > upper bound)",
                      "Conservation concern",
                      "Management concern",
                      "Yield concern",
                      "Residual"))

# Criteria occurrence ####
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





# could be used to simulate a rebuild
sim_Seq_rebuild_list <- mapply(FUN = simSR_goal, 
                          lnalpha = input$lb_p*beta,  #reduced ln_alpha 
                          sigW = input$sigW, 
                          phi = input$phi, 
                          MoreArgs = list(beta = beta, 
                                          age0 = Chinook_age,
                                          lb_sim = 1 / beta,
                                          ub_sim = Inf,
                                          Sims = 1500, 
                                          sigN = 0.1,
                                          sigF = 0,
                                          Hfun = H_goal),
                          SIMPLIFY = FALSE)
# Create dataframe --------------------------------------------------------
sim_Seq_rebuild_df <-
  lapply(sim_Seq_rebuild_list, as.data.frame) %>%
  lapply(function(x) rename(x, lnalpha_red = lnalpha)) %>%
  lapply(function(x){mutate(x, sim = row_number())}) %>%
  do.call("rbind", .) %>%
  mutate(lnalpha = rep(input$lnalpha, each = input_sims)) %>%
  filter(R != 0) %>% #remove rows where population perished
  mutate(resid = S - lb) %>%  
  select(-starts_with("N_age")) %>%
  group_by(lnalpha, sigW, phi) %>%
  mutate(resid_roll = roll_mean(x = resid, 5, align = "right", fill = NA) / lb)



# * Plot time series --------------------------------------------------------
# because of the colors used teal below the lower bound represents times we fished below the lower bound w abundant fish... call the a management concern
# teal above the upper bound represents underutilized yield... call that a yield concern
# while grey below the lower bound represent when the fish were not there to make the goal ... call that a conservation concern
# Note that some phi and sigma W combinations are not sustainable
plot_ts(sim_Seq_rebuild_df, 0.5)
plot_ts(sim_Seq_rebuild_df, 1)
plot_ts(sim_Seq_rebuild_df, 1.5)
plot_ts(sim_Seq_rebuild_df, 2)
data.frame(input, 
           dif_S = sapply(sim_Seq_rebuild_list, function(x) {mean(x$S, na.rm = TRUE)}) - sapply(sim_Seq_N1_list, function(x) {mean(x$S, na.rm = TRUE)}), 
           dif_N = sapply(sim_Seq_rebuild_list, function(x) {median(x$N, na.rm = TRUE)}) - sapply(sim_Seq_N1_list, function(x) {median(x$N, na.rm = TRUE)}),
           Harvest = sapply(sim_Seq_rebuild_list, function(x) {mean(x$U != 0)}))
standard <- sim_Seq_N1_df %>% filter(lnalpha == 1, sigW == .5, phi == .3) %>% select(S)
rebuild <- sim_Seq_rebuild_df %>% filter(lnalpha == 1, sigW == .5, phi == .3) %>% select(S)
plot(standard$S, rebuild$S)
abline(a = 0, b = 1)

sim_Seq_rebuild_df

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
                      "P(N < lower bound)",
                      "P(S < lb, N > lb)",
                      "P(S > upper bound)",
                      "Conservation concern",
                      "Management concern",
                      "Yield concern",
                      "Residual"))

# # 
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
# #Find a value of h that provides negligible yielded for the same beta
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
