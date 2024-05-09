# Simulate SOC listing criteria under a variety of productivity regimes

# Author: Adam Reimer
# Version: 2023-09-15

# Packages
packs <- c("tidyverse", "RcppRoll", "knitr")
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
vec_lnalpha <- seq(.5, 2, length.out = 4)
vec_sigW <- seq(0.25, 1, length.out = 4)
vec_phi <- seq(0, .9, length.out = 4)
beta <- 0.0001
df_power <- data.frame(lnalpha = vec_lnalpha,
                       power = c(0.25, 0.5, 0.6, 0.8))
input <- 
  expand.grid(lnalpha = vec_lnalpha, sigW = vec_sigW, phi = vec_phi) %>%
  mutate(lna_p = lnalpha + (sigW * sigW / 2 / (1 - phi * phi)),
         lb_p = lna_p/beta*(0.5 - 0.07 * lna_p) * 0.8,
         ub_p = lna_p/beta*(0.5 - 0.07 * lna_p) * 1.6) %>%
  left_join(df_power, by = "lnalpha")

Chinook_age <- c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02)

# Simulation base case -------------------------------------------------------------------
#First simulation sigN = 0
#Not very realistic
# Simulation results
# sim_base_grid <- mapply(FUN = simSR_goal, 
#                         lnalpha = input$lnalpha, 
#                         sigW = input$sigW, 
#                         phi = input$phi, 
#                         lb_sim = input$lb_p,
#                         ub_sim = input$ub_p,
#                         MoreArgs = list(beta = beta, 
#                                         age0 = Chinook_age, 
#                                         Sims = 1500, 
#                                         sigN = 0.2,
#                                         sigF = 0,
#                                         Hfun = H_goal),
#                         SIMPLIFY = FALSE)
# saveRDS(sim_base_grid, file = ".\\sim_base_grid.R")
sim_base_grid <- readRDS(file = ".\\sim_base_grid.R")

# * Plot time series --------------------------------------------------------
#writing function here since it has dubious use outside of this single application
plot_ts <- function(dat, lnalpha0){ 
  df <- lapply(dat, as.data.frame) %>%
          lapply(function(x){mutate(x, sim = row_number())}) %>%
          do.call("rbind", .) %>%
          rowwise() %>%
          mutate(N = N_age.1 + N_age.2 + N_age.3 + N_age.4 + N_age.5) %>%
          filter(lnalpha %in% lnalpha0) %>%
          pivot_longer(c(S, N), names_to = "stat", values_to = "value")
  
  limit <- quantile(df$value[df$stat == "N"], .95)
  
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

plot_ts(sim_base_grid, 0.5)
plot_ts(sim_base_grid, 1)
plot_ts(sim_base_grid, 1.5)
plot_ts(sim_base_grid, 2)


# Create dataframe --------------------------------------------------------
sim_base_df <-
  lapply(sim_base_grid, as.data.frame) %>%
    lapply(function(x){mutate(x, sim = row_number())}) %>%
    do.call("rbind", .) %>%
    filter(R != 0) %>%
    mutate(N = rowSums(pick(starts_with("N_age."))),
           cc = ifelse(U == 0 & S <= lb, TRUE, FALSE), #conservation concern defined as missing the goal when the fish were not there
           mc = ifelse(U > 0 & S <= lb, TRUE, FALSE), #management concern = missing the goal due to fishing
           yc = ifelse(S > ub, TRUE, FALSE), #Yield concern = going over the top end
           miss = ifelse(S <= lb, TRUE, FALSE),
           resid = S - lb) %>%  
    select(-starts_with("N_age")) %>%
    group_by(lnalpha, sigW, phi) %>%
    mutate(miss4.5 = ifelse(roll_sum(x = miss, 5, align = "right", fill = NA) >= 4, TRUE, FALSE),
           miss5.5 = ifelse(roll_sum(x = miss, 5, align = "right", fill = NA) >= 5, TRUE, FALSE),
           cc4.5 = ifelse(roll_sum(x = cc, 5, align = "right", fill = NA) >= 5, TRUE, FALSE),
           mc4.5 = ifelse(roll_sum(x = mc, 5, align = "right", fill = NA) >= 5, TRUE, FALSE),
           yc4.5 = ifelse(roll_sum(x = yc, 5, align = "right", fill = NA) >= 5, TRUE, FALSE),
           resid_MD = roll_mean(x = resid, 5, align = "right", fill = NA),
           resid_MP = resid_MD / lb,
           dev = -2 * log(plnorm(R, log(S*exp(lnalpha - beta*S)), sigW)),
           SOC0 = ifelse(cc4.5 == TRUE, "No fish", ifelse(miss4.5 == TRUE, "w Harvest", "No SOC")),
           SOC = factor(SOC0, levels = c("No fish", "w Harvest", "No SOC")))


# Table of concern criteria -----------------------------------------------
sim_base_df %>%
  summarise(N = length(sim),
            below_lb = sum(miss, na.rm = TRUE),
            below_lb_fish = sum(cc, na.rm = TRUE),
            below_lb_harvest = sum(mc, na.rm = TRUE),
            above_ub = sum(yc, na.rm = TRUE),
            miss4.5 = sum(miss4.5, na.rm = TRUE),
            cc4.5 = sum(cc4.5, na.rm = TRUE),
            mc4.5 = sum(mc4.5, na.rm = TRUE),
            yc4.5 = sum(yc4.5, na.rm = TRUE),
            resid = sum(resid_MP < 0, na.rm = TRUE)) %>%
  arrange(phi, sigW, lnalpha) %>%
  kable()

# Criteria occurrence
sim_base_df %>%
  summarise_at(.vars = vars(ends_with("4.5")), ~ mean(., na.rm = TRUE)) %>%
  pivot_longer(ends_with("4.5"), names_to = "Criteria", values_to = "Probability") %>%
  ggplot(aes(x = sigW, y = Probability, color = Criteria)) +
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
plot_resid(sim_base_df, 0.5)
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

# example
S<- 1:30000
plot(S, S*exp(lnalpha - beta*S))
abline(a = 0, b = 1)
abline(v = lb_p, lt = 2)
lines(S, S*exp(lb_p*beta - beta*S), col = "red")
# 33% reduction in productivity in this case
lb_p*beta / lnalpha

plot(S, log(S*exp(lnalpha - beta*S) / S))
abline(v = lb_p, lt = 2)
lines(S, log(S*exp(lb_p*beta - beta*S) / S), col = "red")

# * static reduced alpha -------------------------------------------------------------
static_Seq <- simSR_goal(lb_p*beta, beta, 0, 0,
                         age0 = c("1" = 1),
                         Sims0 = 1,
                         Hfun = H_target,
                         target = lb_p)
static_Seq$U

#some harvest when random variability is included
sim_Seq <- simSR_goal(lb*beta, beta, 0.5, 0,
                      age0 = c("1" = 1),
                      Sims0 = 100,
                      Hfun = H_target,
                      target = lb_p)
hist(sim_Seq$U)
mean(sim_Seq$U == 0) # mostly closed fisheries
mean(sim_Seq$U == 0 & sim_Seq$S < lb) # Often miss lower bound wo fishing
mean(sim_Seq$U > 0)
mean(sim_Seq$N_age * sim_Seq$U) # average yield under reduced productivity
##mean(sim_Seq$N_age * sim_Seq$U) / ?  # ~ 10% of MSY under average productivity conditions FIX

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
sim_Seq_N1_grid <- mapply(FUN = simSR_goal, 
                          lnalpha = input$lb_p*beta,  #reduced ln_alpha 
                          sigW = input$sigW, 
                          phi = input$phi, 
                          lb_sim = input$lb_p, #same bounds
                          ub_sim = input$ub_p,
                          MoreArgs = list(beta = beta, 
                                          age0 = Chinook_age, 
                                          Sims = 1500, 
                                          sigN = 0.1,
                                          sigF = 0,
                                           Hfun = H_goal),
                          SIMPLIFY = FALSE)
#saveRDS(sim_Seq_N1_grid, file = ".\\sim_Seq_N1_grid.R")
sim_Seq_N1_grid <- readRDS(file = ".\\sim_Seq_N1_grid.R")


# * Plot time series --------------------------------------------------------
# because of the colors used teal below the lower bound represents times we fished below the lower bound w abundant fish... call the a management concern
# teal above the upper bound represents underutilized yield... call that a yield concern
# while grey below the lower bound represent when the fish were not there to make the goal ... call that a conservation concern
# Note that some phi and sigma W combinations are not sustainable
plot_ts_red <- function(dat, lnalpha0){ 
  df <- lapply(dat, as.data.frame) %>%
    lapply(function(x){mutate(x, sim = row_number())}) %>%
    mapply(function(x, y) {cbind(x, lnalpha_group = y)}, ., input$lnalpha, SIMPLIFY = FALSE) %>%
    do.call("rbind", .) %>%
    rowwise() %>%
    mutate(N = N_age.1 + N_age.2 + N_age.3 + N_age.4 + N_age.5) %>%
    filter(lnalpha_group %in% lnalpha0) %>%
    pivot_longer(c(S, N), names_to = "stat", values_to = "value")
  
  limit <- quantile(df$value[df$stat == "N"], .95)
  
  ggplot(df, aes(x = sim, y = value, color = stat)) +
    geom_line(alpha = 0.4) +
    geom_hline(aes(yintercept = lb), linetype = 2) +
    geom_hline(aes(yintercept = ub), linetype = 2) +
    scale_y_continuous(limits = c(0, limit)) +
    scale_x_continuous(name = "Simulation year") +
    facet_grid(paste0("\u03C6: ", phi) ~ paste0("\u03C3: ", sigW)) +
    ggtitle(label = bquote(log(alpha): .(lnalpha0)))
}

plot_ts_red(sim_Seq_N1_grid, 0.5)
plot_ts_red(sim_Seq_N1_grid, 1)
plot_ts_red(sim_Seq_N1_grid, 1.5)
plot_ts_red(sim_Seq_N1_grid, 2)


# * Create dataframe --------------------------------------------------------
sim_Seq_N1_df <-
  lapply(sim_Seq_N1_grid, as.data.frame) %>%
  lapply(function(x){mutate(x, sim = row_number())}) %>%
  mapply(function(x, y) {cbind(x, lnalpha_group = y)}, ., input$lnalpha, SIMPLIFY = FALSE) %>%
  do.call("rbind", .) %>%
  filter(R != 0) %>%
  mutate(N = rowSums(pick(starts_with("N_age."))),
         cc = ifelse(U == 0 & S <= lb, TRUE, FALSE), #conservation concern defined as missing the goal when the fish were not there
         mc = ifelse(U > 0 & S <= lb, TRUE, FALSE), #management concern = missing the goal due to fishing
         yc = ifelse(S > ub, TRUE, FALSE), #Yield concern = going over the top end
         miss = ifelse(S <= lb, TRUE, FALSE),
         resid = S - lb) %>%  
  select(-starts_with("N_age")) %>%
  group_by(lnalpha, sigW, phi) %>%
  mutate(miss4.5 = ifelse(roll_sum(x = miss, 5, align = "right", fill = NA) >= 4, TRUE, FALSE),
         miss5.5 = ifelse(roll_sum(x = miss, 5, align = "right", fill = NA) >= 5, TRUE, FALSE),
         cc4.5 = ifelse(roll_sum(x = cc, 5, align = "right", fill = NA) >= 5, TRUE, FALSE),
         mc4.5 = ifelse(roll_sum(x = mc, 5, align = "right", fill = NA) >= 5, TRUE, FALSE),
         yc4.5 = ifelse(roll_sum(x = yc, 5, align = "right", fill = NA) >= 5, TRUE, FALSE),
         resid_MD = roll_mean(x = resid, 5, align = "right", fill = NA),
         resid_MP = resid_MD / lb,
         dev = -2 * log(plnorm(R, log(S*exp(lnalpha - beta*S)), sigW)),
         SOC0 = ifelse(cc4.5 == TRUE, "No fish", ifelse(miss4.5 == TRUE, "w Harvest", "No SOC")),
         SOC = factor(SOC0, levels = c("No fish", "w Harvest", "No SOC")))


# * Table of concern criteria -----------------------------------------------
# estimation error for N introduces lb misses based on fishing and increased ub misses.
sim_Seq_N1_df %>%
  summarise(N = length(sim),
            below_lb = sum(miss, na.rm = TRUE),
            below_lb_fish = sum(cc, na.rm = TRUE),
            below_lb_harvest = sum(mc, na.rm = TRUE),
            above_ub = sum(yc, na.rm = TRUE),
            miss4.5 = sum(miss4.5, na.rm = TRUE),
            cc4.5 = sum(cc4.5, na.rm = TRUE),
            mc4.5 = sum(mc4.5, na.rm = TRUE),
            yc4.5 = sum(yc4.5, na.rm = TRUE),
            resid = sum(resid_MP < 0, na.rm = TRUE)) %>%
  kable()

# Criteria occurrence #### needs lnalpha_group
sim_Seq_N1_df %>%
  group_by(lnalpha_group, lnalpha, sigW, phi) %>%
  summarise_at(.vars = vars(ends_with("4.5")), ~ mean(., na.rm = TRUE)) %>%
  pivot_longer(ends_with("4.5"), names_to = "Criteria", values_to = "Probability") %>%
  ggplot(aes(x = sigW, y = Probability, color = Criteria)) +
  geom_line() +
  facet_grid(phi ~ lnalpha_group, labeller = label_bquote(rows = phi: .(phi), cols = log(alpha): .(lnalpha)))

# * Average 5 year rolling residual as a diagnostic -------------------------
# Average 5 year rolling residual as a diagnostic -------------------------
plot_resid_red <- function(dat, lnalpha0){
  dat %>%
    filter(lnalpha_group == lnalpha0, !is.na(SOC)) %>%
    mutate(resid_pct = resid / lb) %>%
    ggplot(aes(x = resid_pct, fill = SOC)) +
    geom_histogram() +
    scale_x_continuous(limits = c(-1, 1)) + 
    geom_vline(aes(xintercept = 0)) +
    facet_grid(paste0("\u03C6: ", phi) ~ paste0("\u03C3: ", sigW),
               scales = "free_y") +
    ggtitle(label = bquote(log(alpha): .(lnalpha0))) #.(lnalpha) pulls from the parent frame instead of the function's environment
}

#Few concerns. Notice that in some cases "miss" diagnostics are not residual outliers.
plot_resid_red(sim_Seq_N1_df, 0.5)
plot_resid_red(sim_Seq_N1_df, 1)
plot_resid_red(sim_Seq_N1_df, 1.5) plot_resid(sim_base_df, 1.5)
plot_resid_red(sim_Seq_N1_df, 2)





# could be used to simulate a rebuild
sim_Seq_rebuild_grid <- mapply(FUN = simSR_goal, 
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
plot_ts_red(sim_Seq_rebuild_grid, 0.5)
plot_ts_red(sim_Seq_rebuild_grid, 1)
plot_ts_red(sim_Seq_rebuild_grid, 1.5)
plot_ts_red(sim_Seq_rebuild_grid, 2)


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
