#Simulate salmon population dynamics with age-at-maturity and management strategy inputs.
#Arguments:
#   a, b, c: gamma SR function parameters
#   sigW: process error standard deviation on the log scale
#   phi: autocorrelation parameter
#   age0: a named vector where the names are numeric total age formatted as a characters(e,g, '3'), and the vector values are proportional age-at-maturity
#   Sims0: number of simulations to run
#   sigN: estimation error
#   sigF: management error
#   Hfun: a function which provides an annual harvest rate based on the management strategy simulated by the function. 
sim_SRgamma <- function(alpha, beta, gamma, sigW, phi, age0, Sims0, sigN = 0, sigF = 0, Hfun = "none", ...){
  #Check arguments
  stopifnot(Hfun %in% c("none", "rate", "goal", "soc"))
  if(Hfun == "rate"){
    if(!("harvest_rate" %in% names(list(...)))) stop("Provide a harvest rate (harvest_rate).")
    if(list(...)["harvest_rate"] <= 0 | list(...)["harvest_rate"] >= 1) stop("Provide a harvest between 0 and 1.")
  }
  if(Hfun %in% c("goal", "soc")){
    if(!("power" %in% names(list(...)))) stop("Provide a maximum fishing power (power).")
    if(list(...)["power"] <= 0 | list(...)["power"] >= 1) stop("Provide a maximum fishing power between 0 and 1.")
    if(!("lb_goal" %in% names(list(...)))) stop("Provide an escapemnt goal lower bound (lb_goal).")
    if(!("ub_goal" %in% names(list(...)))) stop("Provide an escapemnt goal upper bound (ub_goal).")
  } 
  if(Hfun == "soc"){
    if(!("misses_in" %in% names(list(...)))) stop("Provide the number of years the annual run size must be below the lower bound of the escapement goal to list a SOC (misses_in).")
    if(!("makes_out" %in% names(list(...)))) stop("Provide the number of years the annual run size must be above the lower bound of the escapement goal to de-list a SOC (makes_out).")
    if(!("window" %in% names(list(...)))) stop("Provide the number of years considered during SOC designations (window).")
    if(!("lb_soc" %in% names(list(...)))) stop("Provide an escapemnt goal lower bound for use when listed as a SOC (lb_soc).")
    if(!("ub_soc" %in% names(list(...)))) stop("Provide an escapemnt goal upper bound for use when listed as a SOC (ub_soc).")
  }
  
  #Harvest functions
  H_none <- function(){
    U <- 0
    return(U)
  }

  H_rate <- function(harvest_rate, ...){
    U <- harvest_rate
    return(U)
  }

  H_goal <- function(N_sim, lb_goal, ub_goal, power, ...){
    U_lb <- U_ub <- U <- NA
    
    U_lb <- if(is.na(N_sim) == FALSE & N_sim > lb_goal){(N_sim - lb_goal) / N_sim} else{0} # Find U to fish to lb of goal
    U_ub <- if(is.na(N_sim) == FALSE & N_sim > ub_goal){(N_sim - ub_goal) / N_sim} else{0} # Find U to fish to ub of goal
    U <- min(runif(1, U_ub, U_lb), power)
    
    return(U)
  }
  
  #SOC functions
  flag_soc <- function(N_sim, lb_goal = NA, ...){  
    N_lt_lb <- if(is.na(N_sim) | is.na(lb_goal)){NA} else(if(N_sim < lb_goal){TRUE} else(FALSE))
    return(N_lt_lb)
  }
  
  get_SOC <- function(SOC_sim, flag, misses_in, makes_out, window, ...){
    #stopifnot(!is.na(SOC_sim))
      from_SOC <- function(x){
        if(sum(x, na.rm = TRUE) > (window - makes_out)){"SOC"} else("No_concern")
      }
      from_noconcern <- function(x){
        if(sum(x, na.rm = TRUE) >= misses_in){"SOC"} else("No_concern")
      }
    out <- switch(SOC_sim,
                  "SOC" = from_SOC(flag),
                  "No_concern" = from_noconcern(flag),
                  "error")
    
    return(out)
  }
  
  # age-at-maturity input
  A <- if(is.vector(age0)){length(age0)} else {dim(age0)[2]}
  a.min <- as.numeric(names(age0)[1])
  a.max <- a.min + A - 1
  
  #Set simulation length and create empty objects
  Sims = Sims0 + 100
  
  # Year increment (1 for most stocks, a.min for stocks that mature at one age).
  step <- if(A == 1){a.min}else{1}
  
  # Create empty objects and output vector ids
  R <- E1R <- E2R <- redresid <- rep(0, Sims + step) #Notice: R it length(S)+step
  R_out <- (length(R) - Sims0 + 1):length(R)
  N_age <- matrix(NA, Sims + a.max, A)
  N_age_out <- (dim(N_age)[1] - a.max - Sims0 + 1):(dim(N_age)[1] - a.max)
  S <- N <- N_hat <- epsF <- U <- Ft <- H <- rep(NA, Sims)
  S_out <- (length(S) - Sims0 + 1):length(S)
  if(Hfun == "soc"){
    miss_goal <- lb_manage <- ub_manage <- rep(NA, Sims)
    SOC <- rep("No_concern", Sims)
    }
  
  # Create vectors for time varying inputs
  alpha_vec <- if(length(alpha) == 1){rep(alpha, Sims)}else(c(rep(alpha[1], 100), alpha))
  age <- if(is.vector(age0)){matrix(age0, Sims, A, byrow = TRUE)} else {age0}
  
  # Initial values for incomplete broods
  R0 <- alpha[[1]] * (gamma / beta)^gamma * exp(-gamma) * 0.75 #fraction of Rmax
  for (c in 1:(A + a.min)) { 
    R[c] <- if(c%%step == 0){R0}else{0}
    for (a in 1:A) {
      N_age[c + a.max - (a - 1), (A + 1 - a)] <- age[c, (A + 1 - a)] * R[c]
    }
  }
  
  # Simulate population dynamics
  for (c in (A + a.min):Sims) {
    
    # Calculate values for this brood year
    # Annual Run size
    N[c] <- sum(N_age[c, ])
    N_hat[c] <- N[c]*rlnorm(1, sdlog = sigN)
    
    #SOC status
    if(Hfun == "soc"){
      miss_goal[c] <- flag_soc(N_sim = N[c], ...)
      SOC[c] <- get_SOC(SOC_sim = SOC[c], flag = miss_goal[(c - list(...)$window + 1):(c)], ...)
      lb_manage[c] <- if(!is.na(SOC[c]) & SOC[c] == "SOC"){list(...)$lb_soc} else(list(...)$lb_goal)
      ub_manage[c] <- if(!is.na(SOC[c]) & SOC[c] == "SOC"){list(...)$ub_soc} else(list(...)$ub_goal)
    }
    
    # Calculate harvest rate and instantaneous mortality
    U[c] <- 
      switch(Hfun,
             none = H_none(),
             rate = H_rate(...),
             goal = H_goal(N_sim = N_hat[c], ...),
             soc = H_goal(N_sim = N_hat[c], 
                          lb_goal = lb_manage[c], 
                          ub_goal = ub_manage[c], 
                          power = list(...)$power)
      )
    epsF[c] <- rnorm(1, 0, sigF)
    Ft[c] <- -log(1 - U[c])*exp(epsF[c])
    
    # Spawning abundance
    S[c] <- N[c] * exp(-Ft[c])
    H[c] <- N[c] - S[c]
    
    # Recruitment (step years ahead of the S that produced it)
    E1R[c + step] <- SRgamma(alpha_vec[c], beta, gamma, S[c])
    E2R[c + step] <- E1R[c + step]*exp(phi * redresid[c])
    R[c + step] <- if(is.na(S[c]) | S[c] <= 1){0} else{E2R[c + step]*rlnorm(1, 0, sigW)}
    redresid[c + step] <- log(R[c + step] / E1R[c + step])
    
    # Distribute recruits to future runs
    for (a in 1:A) {
      N_age[c + a.max - (a - 1), (A + 1 - a)] <- age[c, (A + 1 - a)] * R[c + step]
    }
  }
  
  # Return object
  # return(list(R = R[1:10], N = N[1:10], N.age = N_age[1:10, ]))
  out <- data.frame(sim = 1:Sims0, # Fill dataframe with parameter values used in simulation
                    alpha = alpha_vec[(length(alpha_vec) - Sims0):(length(alpha_vec) - 1)],
                    beta = rep(beta, Sims0),
                    gamma = rep(gamma, Sims0),
                    sigW = rep(sigW, Sims0),
                    phi = rep(phi, Sims0),
                    sigN = rep(sigN, Sims0),
                    sigF = rep(sigF, Sims0),
                    S = S[S_out],
                    F = Ft[S_out],
                    U = U[S_out],
                    R = R[R_out], #R is longer than the others!
                    N = N[S_out],
                    N_hat = N_hat[S_out],
                    N_age = N_age[N_age_out, ]) 
  
  if(Hfun %in% c("goal")){
    out <- 
      out %>%
      mutate(lb_goal = rep(list(...)$lb_goal, Sims0),
             ub_goal = rep(list(...)$ub_goal, Sims0)) %>%
      relocate(lb_goal, ub_goal, .after = sigF)
  }
  
  if(Hfun == "soc"){
    out <- 
      out %>%
      mutate(SOC = SOC[S_out],
             lb_goal = lb_manage[S_out],
             ub_goal = ub_manage[S_out]) %>%
      relocate(lb_goal, ub_goal, .after = sigF)
  }
  
  return(out %>% filter(R != 0))
}
sim_SRgamma(exp(1.5), 0.0001, 1, 0.5, 0, c('2' = 1), Sims0 = 100,
            Hfun = "goal", power = 0.3, 
            lb_goal = 3000, ub_goal = 10000)
sim_SRgamma(exp(1.5), 0.0001, 1, 0.5, 0, c('2' = 0.05, '3' = 0.90, '4' = 0.05), Sims0 = 100,
            Hfun = "soc", power = 0.3, 
            lb_goal = 3000, ub_goal = 10000, 
            lb_soc = 4000, ub_soc = 9000,
            misses_in = 4, makes_out = 4, window = 5)



#rep_scenarios_gp0 <- readRDS(file = ".\\rep_scenarios_seqlb_scale_SOC.rds")
rep_scenarios_gp0 %>% 
  filter(rep == 1, lnalpha_1 == 1.5, pct_MSY == 0.9, gamma == 1.3, sigma == 0.5, pct_lb == 1) %>%
  rowwise() %>%
  mutate(Smsy_gamma =
           optimize(
             function(a, b, g, x){
               SRgamma(alpha = a, beta = b, gamma = g, S = x) - x
             },
             interval = c(0, 4 * 1 / beta),
             maximum = TRUE,
             a = a_1,
             b = b,
             g = gamma)$maximum,
         MSY_gamma = SRgamma(alpha = a_1, beta = b, gamma = gamma, S = Smsy_gamma) - Smsy_gamma) %>%
  mutate(lb_gamma =
           optimize(
             function(a, b, g, pct_MSY, MSY_gamma, x){
               ((SRgamma(alpha = a, beta = b, gamma = g, S = x) - x) - pct_MSY * MSY_gamma)^2
             },
             interval = c(Smsy_gamma / 10, Smsy_gamma),
             a = a_1,
             b = b,
             g = gamma,
             pct_MSY = pct_MSY,
             MSY = MSY_gamma)$minimum,
         ub_gamma =
           optimize(
             function(a, b, g, pct_MSY, MSY_gamma, x){
               ((SRgamma(alpha = a, beta = b, gamma = g, S = x) - x) - pct_MSY * MSY_gamma)^2
             },
             interval = c(Smsy_gamma, Smsy_gamma * 3),
             a = a_1,
             b = b,
             g = gamma,
             pct_MSY = 0.7,
             MSY = MSY_gamma)$minimum) %>%
  select(lb_gamma, ub_gamma)


test_goal <- sim_SRgamma(c(rep(0.382, 100), rep(0.125, 900)), 0.00013, 1.3, 0.5, 0, c('3' = 0.25, '4' = 0.5, '5' = 0.25), Sims0 = 900,
            Hfun = "goal", power = 0.3, 
            lb_goal = 3840, ub_goal = 10347)

test_soc <- sim_SRgamma(c(rep(0.382, 100), rep(0.125, 900)), 0.00013, 1.3, 0.5, 0, c('3' = 0.25, '4' = 0.5, '5' = 0.25), Sims0 = 900,
                    Hfun = "soc", power = 0.3, 
                    lb_goal = 3840, ub_goal = 10347, 
                    lb_soc = 4586, ub_soc = 10607,
                    misses_in = 1, makes_out = 4, window = 5)
plot(test_soc$sim, test_soc$S, type = "l", col = "red")
lines(test_goal$sim, test_goal$S)
abline(h = 3840)
abline(h = 4586, col = "red")

 