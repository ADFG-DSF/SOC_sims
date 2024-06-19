# Simulation of the goal setting process

# Author: Adam Reimer
# Version: 2024-06-5

# Packages
packs <- c("tidyverse", "ggforce", "broom", "modelr")
lapply(packs, require, character.only = TRUE)

# source functions
function_files <- list.files(path=".\\functions")
lapply(function_files, function(x) source(paste0(".\\functions\\", x)))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#Pull a simulation to get a bunch of SR pairs
sim_base_list <- readRDS(file = ".\\sim_base_list.R")
input[5, ]
#lnalpha = 1.5, sigW = 0.5, phi = 0.0
sim_base_list <- sim_base_list[[5]]

# get SR pairs for a Ricker regression
data <- 
  data.frame(R = sim_base_list$R,
             S = sim_base_list$S) %>%
  mutate(lnRS = ifelse(is.finite(log(R/S)) & !is.na(log(R/S)), log(R/S), NA),
         sim = row_number()) %>%
  filter(!is.na(lnRS)) %>%
  select(-R)

#bootstrap residuals for each 25 year group. 
years <- 25
boot_reps <- 500
mods <- 
  data %>%
  mutate(group = (sim - 1) %/% years) %>% #number 25 year groups
  group_by(group) %>%
  nest() %>%
  mutate(mod = map(data, ~ lm(lnRS ~ S, data = .)), #regular regression 
         fit = map(mod, ~ .$fitted.values),
         resid = map(mod, ~ .$residuals),
         lnalpha = map_dbl(mod, ~ .$coefficients[1]),
         beta = map_dbl(mod, ~ -.$coefficients[2]),
         sigma = map_dbl(mod, ~ summary(.)$sigma),
         MSY = get_MSY(lnalpha, beta, TRUE, sigma = sigma, phi = 0)) %>%
  slice(rep(1:n(), each = boot_reps)) %>% #duplicate each group
  mutate(boot = rep(1:boot_reps)) %>% #number boots reps in each group
  group_by(group, boot) %>%
  mutate(resamp = map(resid, ~ sample(., replace = TRUE)),
         lnRS_boot = map2(fit, resamp, function(x, y) x + y),
         data_boot = map2(data, lnRS_boot, function(x, y) 
           data.frame(data, lnRS_boot) %>% setNames(c("S", "lnRS", "sim", "lnRS_boot"))),
         mod_boot = map(data_boot, ~ lm(lnRS_boot ~ S, data = .)),
         lnalpha_boot = map_dbl(mod_boot, ~ .$coefficients[1]),
         beta_boot = map_dbl(mod_boot, ~ -.$coefficients[2]),
         sigma_boot = map_dbl(mod_boot, ~ summary(.)$sigma))

#Pull boot params as it's own df for speed.
params <-
  mods %>%
  mutate(Smsy_boot = get_Smsy(lnalpha_boot, beta_boot, TRUE, sigma = sigma_boot, phi = 0)) %>% #by group and boot
  ungroup() %>%
  mutate(thirds_lnalpha = cut(lnalpha_boot, quantile(lnalpha_boot, c(0, 0.33, 0.67, 1))), #overall
         thirds_beta = cut(beta_boot, quantile(beta_boot, c(0, 0.33, 0.67, 1)))) %>%
  filter(beta_boot > 0, Smsy_boot > 0) %>% #add p.value ?
  rowwise() %>%
  mutate(lb = optimise(get_bounds, #by group and boot
                       1:Smsy_boot,
                       lnalpha = lnalpha_boot,
                       beta = beta_boot,
                       pct_MSY = 0.8,
                       correct = TRUE,
                       sigma = sigma_boot,
                       phi = 0)$minimum,
         ub = optimise(get_bounds,
                       Smsy_boot:(Smsy_boot*3),
                       lnalpha = lnalpha_boot,
                       beta = beta_boot,
                       pct_MSY = 0.8,
                       correct = TRUE,
                       sigma = sigma_boot,
                       phi = 0)$minimum) %>%
  select(group, MSY, boot, lnalpha_boot, beta_boot, sigma_boot, lb, ub, thirds_lnalpha, thirds_beta)
#lnalpha and beta are related and it effects estimation of our upper bound.
ggplot(params, aes(x = lnalpha_boot, y = beta_boot)) +
  geom_point()
ggplot(params, aes(x = ub)) +
  geom_histogram() +
  scale_x_continuous(limits = c( 0, 20000)) +
  facet_grid(thirds_lnalpha ~ thirds_beta)

#Histogram of the mean bootstrap parameter estimates for each group
params_group <- 
  params %>%
  group_by(group) %>%
  summarize(lnalpha_boot = median(lnalpha_boot),
            beta_boot = median(beta_boot),
            sigma_boot = median(sigma_boot))
hist(params_group$lnalpha_boot)
hist(params_group$beta_boot)
hist(params_group$sigma_boot)


#Optimum Yield Profile
#max S for PYP plot
Smax <- median(params$lnalpha_boot)/ median(params$beta_boot)

#Create OYP data set
OYP <- 
  params %>% 
  group_by(group, boot) %>%
  mutate(MSY_boot = get_MSY(lnalpha_boot, beta_boot, TRUE, sigma = sigma_boot, phi = 0)) %>%
  slice(rep(1:n(), times = 100)) %>%
  mutate(Sstar = seq(1, Smax, length.out = 100),
         SY = get_SY(lnalpha_boot, beta_boot, Sstar, TRUE, sigma = sigma_boot, phi = 0),
         OYP = SY > MSY_boot * 0.8) %>%
  group_by(group, Sstar) %>%
  summarise(OYP = mean(OYP))

#Function to plot true OYP
SY <- function(S){
  out <- S >= input$lb[5] & S <= input$ub[5]
  as.numeric(out)}
#OYP curve if each 25 year block was use3d to assess OYP probabilities.
OYP_temporal <- mods %>%
  select(lnalpha, beta, sigma) %>%
  mutate(MSY_boot = get_MSY(lnalpha, beta, TRUE, sigma = sigma, phi = 0)) %>%
  slice(rep(1:n(), times = 100)) %>%
  mutate(Sstar = seq(1, Smax, length.out = 100),
         SY = get_SY(lnalpha, beta, Sstar, TRUE, sigma = sigma, phi = 0),
         OYP = SY > MSY_boot * 0.8) %>%
  group_by(Sstar) %>%
  summarise(OYP = mean(OYP))
#OYP curves you might get (one for each 25 year group)
ggplot(OYP, aes(x= Sstar, y = OYP, color = as.character(group))) +
  geom_line() +
  geom_line(data = OYP_temporal, aes(x = Sstar, y = OYP), linewidth = 0.75, inherit.aes = FALSE) +
  geom_function(fun = SY, color = "black", linetype = 2, linewidth = 0.75) +
  guides(color = "none")


#Notes for later
#Need a new estimation procedure...
#only underestimates in both parameters predict yield near upper end of fixed value OYP
temp <- OYP[OYP$OYP == TRUE & OYP$Sstar > 7000, ] %>% 
  group_by(group) %>% 
  summarise(lna = mean(lnalpha_boot), beta = mean(beta_boot))
hist(temp$lna)
hist(temp$beta)

temp2 <- OYP[OYP$OYP == FALSE & OYP$Sstar > 7000, ] %>% 
  group_by(group) %>% 
  summarise(lna = mean(lnalpha_boot), beta = mean(beta_boot))
hist(temp2$lna)
hist(temp2$beta)







#try bayes
#bootstrap residuals for each 25 year group. 
mods_bayes <- 
  data %>%
  mutate(group = (sim - 1) %/% years) %>% #number 25 year groups
  group_by(group) %>%
  filter(group %in% 1:20) %>%
  nest() %>%
  mutate(mod = map(data, ~ stan_glm(lnRS ~ S, data = ., prior = normal(-0.0001, 0.0001))), #regular regression 
         fit = map(mod, ~ .$fitted.values),
         resid = map(mod, ~ .$residuals),
         lnalpha = map_dbl(mod, ~ .$coefficients[1]),
         beta = map_dbl(mod, ~ -.$coefficients[2]),
         sigma = map_dbl(mod, ~ sigma(.)),
         sample = map(mod, ~ as.matrix(.)[sample(1:4000, 500), ]),
         MSY = get_MSY(lnalpha, beta, TRUE, sigma = sigma, phi = 0))

params_bayes <- 
  mods_bayes %>% 
  select(group, sample) %>%
  unnest(cols = c(sample)) %>%
  mutate(Smsy = get_Smsy(sample[, 1], -sample[, 2], TRUE, sigma = sample[, 3], phi = 0))

#Optimum Yield Profile
#max S for PYP plot
Smax_bayes <- median(params$sample[, 1])/ -median(params$sample[, 2])

#Create OYP data set
OYP_bayes <- 
  params_bayes %>% 
  mutate(MSY_boot = get_MSY(sample[, 1], -sample[, 2], TRUE, sigma = sample[, 3], phi = 0)) %>%
  slice(rep(1:n(), each = 100)) %>%
  mutate(Sstar = rep(seq(1, Smax, length.out = 100), times = 500),
         SY = get_SY(sample[, 1], -sample[, 2], Sstar, TRUE, sigma = sample[, 3], phi = 0),
         OYP = SY > MSY_boot * 0.8) %>%
  group_by(group, Sstar) %>%
  summarise(OYP = mean(OYP))

#OYP curves you might get (one for each 25 year group)
ggplot(OYP_bayes, aes(x= Sstar, y = OYP, color = as.character(group))) +
  geom_line() +
  geom_line(data = OYP_temporal, aes(x = Sstar, y = OYP), linewidth = 0.75, inherit.aes = FALSE) +
  geom_function(fun = SY, color = "black", linetype = 2, linewidth = 0.75) +
  guides(color = "none")

#Notes for later
#I still think bayes estimation with a moderately informative lognormal prior on Smax may improve estimation but...
#It also looks like we are simulating too high of fishing power for many fisheries
plot(data$S[51:75], data$R[51:75], xlim = c(0, 35000), ylim = c(0, 35000))
abline(a =0 , b = 1)
