Ricker <- function(lnalpha, beta, S)  {S*exp(lnalpha - beta*S)}

#Simulate salmon population dynamics with age-at-maturity and management strategy inputs.
#Arguments:
#   lnalpha, beta, sigW, phi: Ricker SR parameters
#   age0: a named vector where the names are numeric total age formatted as a characters(e,g, '3'), and the vector values are proportional age-at-maturity
#   Sims0: number of simulations to run
#   sigN: estimation error
#   sigF: management error
#   Hfun: a function which provides an annual harvest rate based on the management strategy simulated by the function. 
simSR_goal <- function(lnalpha, beta, sigW, phi, age0, Sims0, sigN = 0, sigF = 0, Hfun, ...){
  Sims = Sims0 + 100
  
  #format age-at-maturity input and draw indexing values
  A <- length(age0)
  a.min <- as.numeric(names(age0)[1])
  a.max <- a.min + A - 1
  age <- if(is.vector(age0)){matrix(age0, Sims + A + a.min + 1, A, byrow = TRUE)} else {age0}
  # Ricker parameters
  Seq <- lnalpha / beta
  
  # initial values
  R <- E1R <- E2R <- whiteresid <- lnalpha.y <- rep(NA, Sims + A + a.min + 1)
  redresid <- rep(0, Sims + A + a.min + 1)
  N_age <- matrix(NA, Sims + A + a.min + a.max + 1, A)
  for (c in 1:(A + a.min)) { 
    R[c] <- Seq
    for (a in 1:A) {
      N_age[c + a.max - (a - 1), (A + 1 - a)] <- age[c, (A + 1 - a)] * R[c]
    }
  }
  
  S <- N <- N_hat <- epsF <- U <- Ft <- H <- Rhat <- lnRhatShat <- fittedR <- cc <- mc <- yc <- SOC <- lb_manage <- ub_manage <- lb_goal <- ub_goal <- rep(NA, Sims + A + a.min + 1)
  
  # recursive portion...
  for (c in (A + a.min):(Sims + A + a.min)) {
    
    epsF[c] <- rnorm(1, 0, sigF)
    N[c] <- sum(N_age[c, ])
    N_hat[c] <- N[c]*rlnorm(1, sdlog = sigN)
    temp_U <- Hfun(N_sim = N_hat[c], SOC_sim = SOC[c - 1], ...)
    U[c] <- temp_U[[1]]
    lb_manage[c] <- temp_U[[2]]
    ub_manage[c] <- temp_U[[3]]
    Ft[c] <- -log(1 - U[c])*exp(epsF[c])
    S[c] <- N[c] * exp(-Ft[c])
    temp_cc <- get_cc(N_sim = N[c], ...)
    cc[c] <- temp_cc[[1]]
    lb_goal[c] <- temp_cc[[2]]
    mc[c] <- get_mc(S_sim = S[c], ...)[[1]]
    temp_yc <- get_yc(S_sim = S[c], ...)
    yc[c] <- temp_yc[[1]]
    ub_goal[c] <- temp_yc[[2]]
    SOC[c] <- if(c >= (A + a.min + 5)){get_SOC(cc_sim = cc[(c-4):c], mc_sim = mc[(c-4):c], yc_sim = yc[(c-4):c])} else(NA)
    H[c] <- N[c] - S[c]
    
    E1R[c + 1] <- S[c]*exp(lnalpha - beta*S[c])
    E2R[c + 1] <- E1R[c + 1]*exp(phi * redresid[c])
    R[c + 1] <- if(S[c] <= 100){0} else{E2R[c + 1]*rlnorm(1, 0, sigW)}
    redresid[c + 1] <- log(R[c + 1] / E1R[c + 1])
    lnalpha.y[c + 1] <- lnalpha + redresid[c + 1]
    
    for (a in 1:A) {
      N_age[(c + 1) + a.max - (a - 1), (A + 1 - a)] <- age[c, (A + 1 - a)] * R[c + 1]
    }
  }
  
  return(list(lnalpha = rep(lnalpha, Sims0),
              sigW = rep(sigW, Sims0),
              phi = rep(phi, Sims0),
              S = S[(length(S) - Sims0):(length(S) - 1)],
              F = Ft[(length(Ft) - Sims0):(length(Ft) - 1)],
              U = U[(length(U) - Sims0):(length(U) - 1)],
              R = R[(length(R) + 1 - Sims0):length(R)],
              N = N[(length(N) - Sims0):(length(N) - 1)],
              N_hat = N_hat[(length(N_hat) - Sims0):(length(N_hat) - 1)],
              lnalpha.y = lnalpha.y[(length(lnalpha.y) + 1 - Sims0):length(lnalpha.y)],
              N_age = N_age[(dim(N_age)[1] + 1 - A - a.min - Sims0):(dim(N_age)[1] - A - a.min), ],
              lb_goal = lb_goal[(length(lb_goal) - Sims0):(length(lb_goal) - 1)],
              ub_goal = ub_goal[(length(ub_goal) - Sims0):(length(ub_goal) - 1)],
              lb_manage = lb_manage[(length(lb_manage) - Sims0):(length(lb_manage) - 1)],
              ub_manage = ub_manage[(length(ub_manage) - Sims0):(length(ub_manage) - 1)],
              cc = cc[(length(cc) - Sims0):(length(cc) - 1)],
              mc = mc[(length(mc) - Sims0):(length(mc) - 1)],
              yc = yc[(length(yc) - Sims0):(length(yc) - 1)],
              SOC = SOC[(length(SOC) - Sims0):(length(SOC) - 1)]
  ))
}

#Function to determine if a "conservation concern" occurred for each year of a simulation where conservation concern is defined as a 
#run size that failed to meet the goal wo harvest.
#Arguments:
#   N_sim: Run size from the simulation.
#   lb_manage: Lower bound of the escapement goal based on historical SR relationship.
get_cc <- function(N_sim, lb_goal, ...){  
  cc <- if(N_sim < lb_goal){TRUE} else(FALSE)
  return(list(cc, lb = lb_goal))
}

#Function to determine if a "management concern" occurred for each year of a simulation where management concern is defined as a 
#escapement that failed to meet the goal.
#Arguments:
#   S_sim: Run size from the simulation.
#   lb_sim: Lower bound of the escapement goal based on historical SR relationship.
get_mc <- function(S_sim, lb_goal, ...){  
  mc <- if(S_sim < lb_goal){TRUE} else(FALSE)
  return(list(mc, lb = lb_goal))
}

#Function to determine if a "yield concern" occurred for each year of a simulation where management concern is defined as a 
#escapement that exceeded the goal.
#Arguments:
#   S_sim: Run size from the simulation.
#   ub_sim: Lower bound of the escapement goal based on historical SR relationship.
get_yc <- function(S_sim, ub_goal, ...){  
  yc <- if(S_sim > ub_goal){TRUE} else(FALSE)
  return(list(yc, ub = ub_goal))
}

#Function to determine if SOC status for each year of a simulation.
#Arguments:
#   N_sim: Run size from the simulation.
#   lb_sim: Lower bound of the escapemtng goal from the simulation.
get_SOC <- function(cc_sim, mc_sim, yc_sim){  
  SOC <- 
    if(sum(cc_sim, na.rm = TRUE) >= 4){"Conservation"} else(
      if(sum(mc_sim, na.rm = TRUE) >= 4){"Management"} else(
        if(sum(yc_sim, na.rm = TRUE) >= 4){"Yield"} else("No concern")))
  return(SOC)
}

#Function to find U that achieves S at some percentage of bias corrected Smsy
#Arguments:
#   N: run size
#   target0: If a number between 0 and 1 then the percent of Smsy ADF&G will manage for. Defaults to 80% which was Bernard era guidance on escapement goal lower bounds.
#            If a integer greater than 1 then the actual escapement ADF&G will manage for.
#   power: maximum harvest rate for the fishery.
H_target <- function(N_sim, target0 = 0.8, power = .85, ...){
  lna <- with(parent.frame(), lnalpha)
  b <- with(parent.frame(), beta)
  sig <- with(parent.frame(), sigW)
  phi <- with(parent.frame(), phi)
  lnap <- lna + 0.5 * sig * sig / (1 - phi * phi)
  Smsy <- lna / b * (0.5 - 0.07 * lna)
  Smsyp <- lnap / b * (0.5 - 0.07 * lnap)
  
  if(target0 >= 0 & target0 <= 1){target <- Smsyp * target0} else{target = target0}
  
  U <- c(if(N_sim > target){min((N_sim - (target)) / N_sim, power)} else{0})
  
  return(list(U, target = target, power = power))
}

# #Function to find U that archives S somewhere within the goal
#Arguments:
#   N: run size
#   lb_manage: Lower bound of the actual escapement ADF&G will manage for.
#   ub_manage: Upper bound of the actual escapement ADF&G will manage for.
#   power: maximum harvest rate for the fishery.
H_goal <- function(N_sim, lb_manage, ub_manage, power = .85, ...){
  U_lb <- U_ub <- U <- NA

  U_lb <- if(N_sim > lb_manage){(N_sim - lb_manage) / N_sim} else{0} # Find U to fish to lb of goal
  U_ub <- if(N_sim > ub_manage){(N_sim - ub_manage) / N_sim} else{0} # Find U to fish to ub of goal
  U <- min(runif(1, U_ub, U_lb), power)

  return(list(U, lb = lb_manage, ub = ub_manage, power = power))
}

# #Function to find U that archives S somewhere within the goal
#Arguments:
#   N: run size
#   lb_manage: Lower bound of the actual escapement ADF&G will manage for.
#   ub_manage: Upper bound of the actual escapement ADF&G will manage for.
#   power: maximum harvest rate for the fishery.
H_soc <- function(N_sim, SOC_sim, lb_goal, ub_goal, lb_manage, ub_manage, power = .85, ...){
  U_lb <- U_ub <- U <- NA
  
  U_lb <- if(SOC_sim %in% "Conservation"){
      if(N_sim > lb_manage){(N_sim - lb_manage) / N_sim} else{0}
    }
    else{
      if(N_sim > lb_goal){(N_sim - lb_goal) / N_sim} else{0} # Find U to fish to lb of goal
    }
  U_ub <- if(SOC_sim %in% "Conservation"){
      if(N_sim > ub_manage){(N_sim - ub_manage) / N_sim} else{0}
    } 
    else{
      if(N_sim > ub_goal){(N_sim - ub_goal) / N_sim} else{0} # Find U to fish to ub of goal
    }
  U <- min(runif(1, U_ub, U_lb), power)
  
  return(list(U, lb = lb_manage, ub = ub_manage, power = power))
}

