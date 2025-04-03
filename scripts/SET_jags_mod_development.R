# SET threshold

# Author: Adam Reimer
# Version: 2025-02-13

# Packages
packs <- c("tidyverse", "jagsUI", "mgcv", "ggpubr", "flextable")
lapply(packs, require, character.only = TRUE)

# source functions
function_files <- list.files(path=".\\functions")
lapply(function_files, function(x) source(paste0(".\\functions\\", x)))
# Used this code to experiment with the jags code.
# Mostly used it to troubleshoot one dataset at a time.
# Sloppy but I left it that way on purpose to see some of the model development.
--------------------------------------------------------------------------------------
# Simulation scenarios
# Note this file has all of the scenarios and 50 replicate simulated data sets.
# also has estimation results from an old model but i'm try to improve upon that.
rep_scenarios_gp <- readRDS(file = ".\\rep_scenarios_gp2.rds")


# TEST DATA SET -----------------------------------------------------------
test_data <- 
  rep_scenarios_gp %>% 
  filter(lnalpha == 2, sigma == 0.5, pct_lb == 0.9, rep == 25, gamma == 1.3)


# Simple Model ------------------------------------------------------------
# Simple priors - first thoughts wo much thought
# naming convention gamma_RS_change: gamma spawner-recruit, on logRS scale, changepoint model 
# i.e (pick a point the time series where the productivity regime changes from high to low)
out_simple <- 
  jags(data = test_data$data_jags[[1]],
       parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d"),
       model.file = ".\\scripts\\gamma_RS_change.txt",
       n.chains = 3,
       n.iter = 5e4,
       n.burnin = 1e4,
       n.thin = 480, 
       parallel = TRUE)

# previous results (and with trial and error) indicate we estimate gamma poorly. 
# originally I thought this was ok because estimation of the entire curve seems ok 
# but mostly when productivity is very low.
# Once I realized that I started this code to figure it out.
summary(out_simple)
plot(out_simple$sims.list$lnalpha[,1], out_simple$sims.list$beta[,1])
plot(out_simple$sims.list$lnalpha[,1], out_simple$sims.list$gamma[,1])
plot(out_simple$sims.list$gamma[,1], out_simple$sims.list$beta[,1])

plot(out_simple$sims.list$lnalpha[,1], out_simple$sims.list$lnalpha[,2])
plot(out_simple$sims.list$lnalpha[,2], out_simple$sims.list$beta[,2])
# This is the one that is gross. and why I refined the model. 
# Why are we estimating 2 things wtih this kind of relationship?
plot(out_simple$sims.list$gamma[,2], out_simple$sims.list$lnalpha[,2])
# I also like how spread out the estiamtes are. for at least some datasets
hist(out_simple$sims.list$gamma[, 2])
hist((out_simple$sims.list$gamma[, 2] - 1) / out_simple$sims.list$beta[, 2])


# SET Model ------------------------------------------------------------
# My first impulse was to remove gamma from the model after I realized I could estimate 
# the SET directly
# naming convention gamma_RS_change_set: gamma spawner-recruit, on logRS scale, changepoint model,
# estimate the set directly
out_set <- 
  jags(data = test_data$data_jags[[1]],
       parameters.to.save = c("lnalpha", "beta", "set", "sigma", "y_d", "delta"),
       model.file = ".\\scripts\\gamma_RS_change_set.txt",
       n.chains = 3,
       n.iter = 5e4,
       n.burnin = 1e4,
       n.thin = 480, 
       parallel = TRUE)

summary(out_set)
plot(out_set$sims.list$lnalpha[,1], out_set$sims.list$beta[,1])
plot(out_set$sims.list$lnalpha[,1], out_set$sims.list$set[,1])
plot(out_set$sims.list$set[,1], out_set$sims.list$beta[,1])

plot(out_set$sims.list$lnalpha[,1], out_set$sims.list$lnalpha[,2])
plot(out_set$sims.list$lnalpha[,2], out_set$sims.list$beta[,2])
plot(out_set$sims.list$set[,2], out_set$sims.list$beta[,2])
# not a ton better
plot(out_set$sims.list$set[,2], out_set$sims.list$lnalpha[,2])
# and we see the same pattern if we back calculate gamma
plot(out_set$sims.list$set[,2] * out_set$sims.list$beta[,2] + 1, out_set$sims.list$lnalpha[,2])

# Scale Model ------------------------------------------------------------
# Turns out there is an algebraic relationship between lnalpha and gamma in the gamma 
# spawner-recruit model provided Rmax and Smax are fixed
# lnalpha = log(Rmax) - gamma * [log(Smax) - 1]
# This model is inspired by using that relationship to reduce the estimation burden on lnalpha
# Instead of estimating both the scale and the slope of the relationship we will tell the model the 
# slope and only estimate the scale
# naming convention gamma_RS_change_scale: gamma spawner-recruit, on logRS scale, changepoint model,
# only estimate the scale... use the lnalpha ~ gamma relationship

# * theoretical relationship ------------------------------------------------
# build a dataset
temp_data = 
  expand.grid(
    gamma = seq(1, 1.6, 0.1), 
    lnalpha = c(0.375, 0.75, 1, 1.25, 1.5, 1.75)) %>%
  mutate(beta = 0.0001,
         Rmax = exp(lnalpha) / beta / exp(1)) %>%
  rowwise() %>%
  mutate(a = gamma_par(1 / beta, exp(lnalpha)/beta/exp(1), gamma)[[1]],
         b = gamma_par(1 / beta, exp(lnalpha)/beta/exp(1), gamma)[[2]])
#plot it with a regression line
temp_data %>%
  ggplot(aes(x = gamma, y = log(a), color = as.character(lnalpha))) +
  geom_point() +
  stat_smooth(method = lm)
#note that the relationship is exact
lm(log(a) ~ gamma, temp[temp$lnalpha == 0.75, ]) %>% summary()
# If you do the algebra in intercept is the log of Rmax
log(unique(temp_data[temp_data$lnalpha == 0.75, ]$Rmax))
# and the slope is the -[log(Smax) - 1]
-(log(1/unique(temp_data[temp_data$lnalpha == 0.75, ]$beta))-1)


# * run model -------------------------------------------------------------
out_scale <- 
  jags(data = test_data$data_jags[[1]],
       parameters.to.save = c("lnalpha", "beta", "gamma", "sigma", "y_d", "scale"),
       model.file = ".\\scripts\\gamma_RS_change_scale.txt",
       n.chains = 3,
       n.iter = 5e4,
       n.burnin = 1e4,
       n.thin = 480, 
       parallel = TRUE)

# It seems like we recover the true gamma better with this model
# will need to prove it by running all of the replicate datasets.
test_data
summary(out_scale)
plot(out_scale$sims.list$lnalpha[,1], out_scale$sims.list$beta[,1])
plot(out_scale$sims.list$lnalpha[,1], out_scale$sims.list$gamma[,1])
plot(out_scale$sims.list$gamma[,1], out_scale$sims.list$beta[,1])

plot(out_scale$sims.list$lnalpha[,1], out_scale$sims.list$lnalpha[,2])
plot(out_scale$sims.list$lnalpha[,2], out_scale$sims.list$beta[,2])

# Still gross but this time the slope is a direct calculation
plot(out_scale$sims.list$gamma[,2], out_scale$sims.list$lnalpha[,2])
#theoretical line
abline(log(test_data$Rmax_red), -(log(test_data$Smax) - 1), col = "green")
#model implied line
abline(out_scale$q50$scale, -(log(1 / out_scale$q50$beta[1]) - 1))
#actual estimated line
(coef <- lm(out_scale$sims.list$lnalpha[,2] ~ out_scale$sims.list$gamma[,2]))
abline(coef$coefficients[1], coef$coefficients[2], col = 'red')
# Marginal improvements in spread
hist(out_scale$sims.list$gamma[, 2])
hist((out_scale$sims.list$gamma[, 2] - 1) / out_scale$sims.list$beta[, 2])

# perfect relationship between lnalpha and gamma
# posterior Rmax
Rmax <- 
  exp(out_scale$sims.list$lnalpha[,2]) * 
  (out_scale$sims.list$gamma[,2] / out_scale$sims.list$beta[,2])^out_scale$sims.list$gamma[,2] *
  exp(-out_scale$sims.list$gamma[,2])
hist(Rmax)
# Posterior Smax
Smax <- out_scale$sims.list$gamma[,2] / out_scale$sims.list$beta[,2]
hist(Smax)
# lnalpha from linear model exact match for posterior estiamtes
lnalpha <- log(Rmax) - out_scale$sims.list$gamma[,2]*(log(Smax) -1)
plot(lnalpha, out_scale$sims.list$lnalpha[, 2])
