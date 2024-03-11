# Simulate SOC listing criteria under a variety of productivity regimes

# Author: Adam Reimer
# Version: 2023-09-15

# Packages
packs <- c("tidyverse")
lapply(packs, require, character.only = TRUE)

# source functions
function_files <- list.files(path=".\\functions")
lapply(function_files, function(x) source(paste0(".\\functions\\", x)))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# age-at-maturity 3-7 years -----------------------------------------------
#Kenai Chinook posterior for finding reasonable population parameter estimates
# post <- readRDS(file = "S:\\RTS\\Reimer\\KenaiSRA\\posts\\post_run_test3.rds")
# lnalpha.y <- post$q50$lnalpha.vec[,2]
# median(lnalpha.y[(length(lnalpha.y)-4):length(lnalpha.y)])
# median(lnalpha.y[(length(lnalpha.y)-9):length(lnalpha.y)])

# Base Case -----------------------------------------------
# * static base case -------------------------------------------------------------
# ** constant SR params -----------------------------------------------------
lnalpha <- 1.5
beta <- 0.0001
sigW <- 0.5

static_base <- simSR_goal(lnalpha, beta, 0, 0,
                       age0 = c("1" = 1),
                       Sims0 = 1,
                       Hfun = H_target,
                       target = .8)
(MSY_base <- static_base$N_age * static_base$U) #MSY for base case


# * simulation base case -------------------------------------------------------------------
# Is there concern when a stock falls below the escapement goal wo a change in SR parameters?
# How often does that occur?
# ** inputs -----------------------------------------------------------
# Grid of SR parameter values
vec_lnalpha <- seq(.5, 2, length.out = 4)
vec_sigW <- seq(0.25, 1, length.out = 4)
vec_phi <- seq(.15, .9, length.out = 4)
input <- 
  expand.grid(lna = vec_lnalpha, s = vec_sigW, p = vec_phi) %>%
  mutate(lna_p = lna + (s * s / 2 / (1 - p * p)),
         lb_p = lna_p/beta*(0.5 - 0.07 * lna_p) * 0.8)

arg_age <- c('3' = 0.1, '4' = 0.2, '5' = 0.3, '6' = 0.38, '7' = 0.02)

# Simulation results
sim_base_grid <- mapply(FUN = simSR_goal, 
                        lnalpha = input$lna, 
                        sigW = input$s, 
                        phi = input$p, 
                        target = input$lb_p,
                        MoreArgs = list(beta = beta, age0 = arg_age, Sims = 1500, Hfun = H_target),
                        SIMPLIFY = FALSE)
#annual restrictions
restrict_base_grid <- lapply(sim_base_grid, function(x) x$U == 0)
#SOC restrictions. Missed lower goal in 4 of 5 consecutive years.
soc_base_grid <- lapply(restrict_base_grid, function(x){
                        temp <- rep(NA, 4)
                        for(i in 5:length(x)) if(sum(x[(i - 4):i]) >= 4){temp[i] <- TRUE} 
                        else{temp[i] <- FALSE}
                        temp})
#SOC probability for each parameter combination in the grid.
soc_base <- sapply(soc_base_grid, mean, na.rm = TRUE)

# # glm: restrictions ~ lna_y
# misses <- unlist(soc_base_grid[49:64])
# lna_y <-
#   lapply(sim_base_grid[49:64], `[`, "lnalpha.y") %>%
#     unlist(recursive = FALSE) %>%
#     lapply(`[`, 5:1500) %>%
#     unlist()
# 
# mod <- glm(misses ~ lna_y, family = "binomial")
# preds <- predict(mod, data.frame(lna_y = seq(0, 2, length.out = 20)), type = "response")
# plot(seq(0, 2, length.out = 20), preds)

#Probability of a SOC recommendation for each grid combination
#SOC rare w constant SR params although some king stocks may get there.
#High variability and autocorrelation are needed, particularly when productivity is high.
cbind(input, soc = soc_base) %>%
  ggplot(aes(x = s, y = p, fill = soc)) +
  geom_tile() +
  facet_wrap(. ~ lna, nrow = 2, ncol = 2, labeller = label_bquote(log(alpha): .(lna))) +
  scale_fill_gradient2(low = "blue", 
                       high = "orange", 
                       mid = "black",
                       midpoint = 0.5,
                       limits = c(0, 1)) +
  scale_x_continuous(name = expression("\u03C3"), breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(name = expression("\u03C6"), breaks = seq(0, 1, by = 0.25))

#Distribution of ln_alpha_y for years w and wo a SOC recommendation.
lapply(sim_base_grid, as.data.frame) %>%
  mapply(function(x, y) {cbind(x, miss = y)}, ., soc_base_grid, SIMPLIFY = FALSE) %>%
  do.call("rbind", .) %>%
  filter(!is.na(miss)) %>%
  ggplot(aes(x = miss, y = lnalpha.y)) +
    #geom_boxplot(outlier.shape = NA) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_hline(aes(yintercept = lnalpha), linetype = 2) +
    scale_y_continuous(name = paste0("Annual log(", expression("\u03B1"), ")"), limits = c(-4, 6)) +
    scale_x_discrete(name = "Was the Run Below the Lower Bound of the Goal in 4 of the last 5 years?") +
    facet_wrap(. ~ lnalpha, nrow = 4, ncol = 4, labeller = label_bquote("mean log(alpha)": .(lnalpha)))





# Reduced alpha: lb = Seq  -------------------------------------
#What if productivity changes but our goal remains based on historic productivity?
#Note S_eq = lnalpha/beta so lnalpha = goal_lb*beta is a productivity change that would reduce yield to 0 at the lower bound.
# lower bound goal which is 0.8% of Smsy
lb <- lnalpha/beta*(0.5 - 0.07 * lnalpha) * 0.8
#Assume goal set using corrected lnalpha
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
mean(sim_Seq$U == 0 & sim_Seq$S < lb) # Often miss lower bound
mean(sim_Seq$U > 0)
mean(sim_Seq$U > 0 & round(sim_Seq$S) == lb)
mean(sim_Seq$N_age * sim_Seq$U) # average yield
mean(sim_Seq$N_age * sim_Seq$U) / MSY_base # ~ 10% of MSY under average productivity conditions

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
#reduced productivity goal (not used in simulation)
input$lb_p*beta/beta * (0.5 - 0.07 * input$lb_p*beta) *.8
# % of average productivity goal ~ 38%
input$lb_p*beta/beta * (0.5 - 0.07 * input$lb_p*beta) *.8 / input$lb_p

#simulation results
sim_Seq_grid <- mapply(FUN = simSR_goal, 
                       lnalpha = input$lb_p*beta,  #reduced ln_alpha
                       sigW = input$s, 
                       phi = input$p, 
                       target0 = input$lb_p, #same lb
                       MoreArgs = list(beta = beta, age0 = arg_age, Sims0 = 1500, Hfun = H_target),
                       SIMPLIFY = FALSE)
#annual restrictions
restrict_Seq_grid <- lapply(sim_Seq_grid, function(x) x$U == 0)
#SOC restrictions. Missed lower goal in 4 of 5 consecutive years.
soc_Seq_grid <- lapply(restrict_Seq_grid, function(x){
  temp <- rep(NA, 4)
  for(i in 5:length(x)) if(sum(x[(i - 4):i]) >= 4){temp[i] <- TRUE} 
  else{temp[i] <- FALSE}
  temp})
#SOC probability for each parameter combination in the grid.
soc_Seq <- sapply(soc_Seq_grid, mean, na.rm = TRUE)



#Probability of a SOC recommendation for each grid combination
#bc the goal is not too high from current productivity it's pretty easy to hit SOC w low variability particularly under low productivity stocks.
#adds a second way to get there...in addition to high phi and sigma.
cbind(input, soc = soc_Seq) %>%
  ggplot(aes(x = s, y = p, fill = soc)) +
  geom_tile() +
  facet_wrap(. ~ lna, nrow = 2, ncol = 2, labeller = label_bquote(log(alpha): .(lna))) +
  scale_fill_gradient2(low = "blue", 
                       high = "orange", 
                       mid = "black",
                       midpoint = 0.5,
                       limits = c(0, 1)) +
  scale_x_continuous(name = expression("\u03C3"), breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(name = expression("\u03C6"), breaks = seq(0, 1, by = 0.25))

#Distribution of ln_alpha_y for years w and wo a SOC recommendation.
#Productivity low with approximately 50% of runs failing to meet the goal.
lapply(sim_Seq_grid, as.data.frame) %>%
  mapply(function(x, y) {cbind(x, miss = y)}, ., soc_Seq_grid, SIMPLIFY = FALSE) %>%
  mapply(function(x, y) {cbind(x, lnalpha_group = y)}, ., input$lna, SIMPLIFY = FALSE) %>%
  do.call("rbind", .) %>%
  filter(!is.na(miss)) %>%
  ggplot(aes(x = miss, y = lnalpha.y)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_hline(aes(yintercept = lnalpha_group), linetype = 2) +
  scale_y_continuous(name = paste0("Annual log(", expression("\u03B1"), ")"), limits = c(-4, 6)) +
  scale_x_discrete(name = "Was the Run Below the Lower Bound of the Goal in 4 of the last 5 years?") +
  facet_wrap(. ~ lnalpha_group, nrow = 4, ncol = 4, labeller = label_bquote("mean log(alpha)": .(lnalpha_group)))




# # Pinks: age-at-maturity 2 years  -----------------------------------------------
# # * Constant SR params -----------------------------------------------------
# #Is there concern when a stock falls below the escapement goal wo a change in SR parameters?
# # ** grid simulation -------------------------------------------------------------------
# #SOC rare w constant SR params.
# # *** inputs -----------------------------------------------------------
# #Grid of SR parameter values
# vec_lnalpha_pinks <- seq(1.5, 3, length.out = 10)
# vec_sigW_pinks <- seq(0, 1.2, length.out = 10)
# vec_phi_pinks <- c(-0.5, -0, 0.5, 0.75)
# input_pinks <- 
#   expand.grid(lna = vec_lnalpha_pinks, s = vec_sigW_pinks, p = vec_phi_pinks) %>%
#   mutate(lna_p = lna + (s * s / 2 / (1 - p * p)),
#          lb_p = lna_p/beta*(0.5 - 0.07 * lna_p) * 0.8)
# 
# sim_basepinks_grid <- mapply(FUN = simSR_goal, 
#                              lnalpha = input_pinks$lna, 
#                              sigW = input_pinks$s, 
#                              phi = input_pinks$p, 
#                              target0 = input_pinks$lb_p,
#                              MoreArgs = list(beta = beta, age0 = c('1' = 0.2, '2' = 0.98), Sims0 = 1500, Hfun = H_target),
#                              SIMPLIFY = FALSE)
# restrict_basepinks_grid <- lapply(sim_basepinks_grid, function(x) x$U == 0)
# soc_basepinks_grid <- sapply(restrict_basepinks_grid, function(x){
#   temp <- NA
#   for(i in 5:length(x)) if(sum(x[(i - 4):i]) >= 4){temp[i - 4] <- TRUE} 
#   else{temp[i - 4] <- FALSE}
#   mean(temp)}
# )
# cbind(input_pinks, soc = soc_basepinks_grid) %>%
#   ggplot(aes(x = lna, y = s, fill = soc)) +
#   geom_tile() +
#   facet_wrap(. ~ p, nrow = 2, ncol = 2, labeller = label_bquote(phi: .(p))) +
#   scale_fill_gradient2(low = "blue", 
#                        high = "orange", 
#                        mid = "black",
#                        midpoint = 0.5,
#                        limits = c(0, 1)) +
#   scale_x_continuous(name = "Log(\u03B1)", breaks = seq(1.5, 3, by = 0.5)) +
#   scale_y_continuous(name = expression("\u03C3"), breaks = seq(0, 1.2, by = 0.1))

# # *delta: lower bound = Seq  -------------------------------------
# #Real concern is accompanied by changing SR parameters.
# #Note S_eq = lnalpha/beta so lnalpha = goal_lb*beta is a productivity change that would reduce yield to 0.
# # ** grid simulation -------------------------------------------------------------------
# sim_Seqpinks_grid <- mapply(FUN = simSR_goal, 
#                             lnalpha = input_pinks$lb_p*beta, 
#                             sigW = input_pinks$s, 
#                             phi = input_pinks$p, 
#                             target0 = input_pinks$lb_p, 
#                             MoreArgs = list(beta = beta, age0 = arg_age, Sims0 = 1000, Hfun = H_target),
#                             SIMPLIFY = FALSE)
# restrict_Seqpinks_grid <- lapply(sim_Seqpinks_grid, function(x) x$U == 0)
# soc_Seqpinks_grid <- sapply(restrict_Seqpinks_grid, function(x){
#   temp <- NA
#   for(i in 5:length(x)) if(sum(x[(i - 4):i]) >= 4){temp[i - 4] <- TRUE} 
#   else{temp[i - 4] <- FALSE}
#   mean(temp)}
# )
# cbind(input_pinks, soc = soc_Seqpinks_grid) %>%
#   ggplot(aes(x = lna, y = s, fill = soc)) +
#   geom_tile() +
#   facet_wrap(. ~ p, nrow = 2, ncol = 2, labeller = label_bquote(phi: .(p))) +
#   scale_fill_gradient2(low = "blue", 
#                        high = "orange", 
#                        mid = "black",
#                        midpoint = 0.5,
#                        limits = c(0, 1)) +
#   scale_x_continuous(name = "Log(\u03B1)", breaks = seq(1.5, 3, by = 0.5)) +
#   scale_y_continuous(name = expression("\u03C3"), breaks = seq(0, 1.2, by = 0.1))
# 
# 
# soc_Seqpinks_grid9yrs <- sapply(restrict_Seqpinks_grid, function(x){
#   temp <- NA
#   for(i in 9:length(x)) if(sum(x[(i - 8):i]) >= 8){temp[i - 8] <- TRUE} 
#   else{temp[i - 8] <- FALSE}
#   mean(temp)}
# )
# cbind(input_pinks, soc = soc_Seqpinks_grid9yrs) %>%
#   ggplot(aes(x = lna, y = s, fill = soc)) +
#   geom_tile() +
#   facet_wrap(. ~ p, nrow = 2, ncol = 2, labeller = label_bquote(phi: .(p))) +
#   scale_fill_gradient2(low = "blue", 
#                        high = "orange", 
#                        mid = "black",
#                        midpoint = 0.5,
#                        limits = c(0, 1)) +
#   scale_x_continuous(name = "Log(\u03B1)", breaks = seq(1.5, 3, by = 0.5)) +
#   scale_y_continuous(name = expression("\u03C3"), breaks = seq(0, 1.2, by = 0.1))
# 
# 
# 
# 
# 
# #Shnute params
# 
# h <- rep(NA, 101)
# alpha <- rep(NA, 100)
# c <- rep(NA, 100)
# 
# for(i in 1:100){
#   h[1] <- 0.5
#   alpha[i] <- exp(h[i])/(1 - h[i])
#   h[i + 1] <- h[i] + (1 - h[i])/(2 - h[i])*log(exp(1.5)/alpha[i])
# }
# h
# h[100]^2 / ( 1 - h[100]) / 0.0001
# 
# ricker <- function(s) s*exp(1.5 - 0.0001*s)
# ricker_alt <- function(s, h = 0.64, c = 8484) s/(1 - h) * exp(h - h^2/(1 - h)*s/c)
# s <- 1:30000
# plot(s, ricker(s), type = 'l')
# abline(0, 1)
# lines(ricker_alt(s, 0.595, 8757), col = "red")
# lines(ricker_alt(s, 0.1, 100), col = "green")
# 
# 
# test <- 
#   data.frame(C = 8757, 
#            h = seq(.4, .9, length.out = 4)) %>%
#   mutate(a = exp(h) / (1 - h),
#          b = h^2 / ((1 - h)*C))
# test
# lines(s, ricker_alt(s, test$h[1], test$C[1]))
# lines(s, ricker_alt(s, test$h[2], test$C[2]))
# lines(s, ricker_alt(s, test$h[3], test$C[3]))
# lines(s, ricker_alt(s, test$h[4], test$C[4]))
# lines(s, ricker_alt(s, test$h[5], test$C[5]))
# lines(s, ricker_alt(s, test$h[6], test$C[6]))
# lines(s, ricker_alt(s, test$h[7], test$C[7]))
# lines(s, ricker_alt(s, test$h[8], test$C[8]))
# 
# 1.5/0.0001*(0.5-0.07*1.5)
# ricker(5925) - 5925
# (ricker(5925) - 5925) / ricker(5925)
# 
# 
# 
# 
# 
# # Smsy <- 1.5/0.0001*(0.5 - 0.07 * 1.5)
# # MSY <- Smsy * exp(1.5 - 0.0001 * Smsy) - Smsy
# # lna2 <- log((MSY / 2 + Smsy) / Smsy) + 0.0001 * Smsy
# # simSR_goal(lna2, 0.0001, 0, 0,
# #            age0 = c("1" = 1),
# #            Sims0 = 25,
# #            Hfun = H_target,
# #            target = Smsy)
# # 
# # simSR_goal(1.5, 0.0001, 0, 0,
# #            age0 = c("1" = 1),
# #            Sims0 = 5,
# #            Hfun = H_target,
# #            target = Smsy)