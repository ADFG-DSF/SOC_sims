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
  mutate(Smsy_boot = get_Smsy(lnalpha_boot, beta_boot, TRUE, sigma = sigma_boot, phi = 0)) %>%
  filter(beta_boot > 0, Smsy_boot > 0) %>% #add p.value ?
  select(group, MSY, boot, lnalpha_boot, beta_boot, sigma_boot)


#Histogram of the mean bootstrap parameter estimates for each group
params_group <- 
  params %>%
  group_by(group) %>%
  summarize(lnalpha_boot = mean(lnalpha_boot),
            beta_boot = mean(beta_boot),
            sigma_boot = mean(sigma_boot))
hist(params_group$lnalpha_boot)
hist(params_group$beta_boot)
hist(params_group$sigma_boot)
ggplot(params_group, aes(x = lnalpha_boot, y = beta_boot, color = sigma_boot)) +
  geom_point()

#Optimum Yield Profile
#max S for PYP plot
Smax <- mean(params$lnalpha_boot)/ mean(params$beta_boot) * 2

#Create OYP data set
OYP <- 
  params %>% 
  slice(rep(1:n(), times = 100)) %>%
  mutate(Sstar = seq(1, Smax, length.out = 100),
         MSY_boot = get_MSY(lnalpha_boot, beta_boot, TRUE, sigma = sigma_boot, phi = 0),
         SY = get_SY(lnalpha_boot, beta_boot, Sstar, TRUE, sigma = sigma_boot, phi = 0),
         OYP = SY > MSY_boot * 0.8) %>%
  group_by(group, Sstar) %>%
  summarise(OYP = mean(OYP))

#Function to plot true OYP
SY <- function(S){
  out <- (Ricker(lnalpha = input$lna_p[5], beta = input$beta[5], S) - S) > 
    get_MSY(input$lna_p[5], input$beta[5], FALSE) * 0.80
  as.numeric(out)}
#OYP curves you might get (one for each 25 year group)
ggplot(OYP, aes(x= Sstar, y = OYP, color = as.character(group))) +
  geom_line() +
  geom_function(fun = SY, color = "black") +
  guides(color = "none")

