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
  
  S <- N <- N_hat <- epsF <- U <- lb <- ub <- Ft <- H <- Rhat <- lnRhatShat <- fittedR <- cc <- mc <- yc <- SOC <- rep(NA, Sims + A + a.min + 1)
  
  # recursive portion...
  for (c in (A + a.min):(Sims + A + a.min)) {
    
    epsF[c] <- rnorm(1, 0, sigF)
    N[c] <- sum(N_age[c, ])
    N_hat[c] <- N[c]*rlnorm(1, sdlog = sigN)
    temp <- Hfun(N_sim = N_hat[c], ...)
    U[c] <- temp[[1]]
    lb[c] <- temp[[2]]
    ub[c] <- temp[[3]]
    Ft[c] <- -log(1 - U[c])*exp(epsF[c])
    S[c] <- N[c] * exp(-Ft[c])
    cc[c] <- get_cc(N_sim = N[c], ...)
    mc[c] <- get_mc(S_sim = S[c], ...)
    yc[c] <- get_yc(S_sim = S[c], ...)
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
              lb = lb[(length(lb) - Sims0):(length(lb) - 1)],
              ub = ub[(length(ub) - Sims0):(length(ub) - 1)],
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
#   lb_sim: Lower bound of the escapemtng goal from the simulation.
get_cc <- function(N_sim, lb_sim, ...){  
  cc <- if(N_sim < lb_sim){TRUE} else(FALSE)
  return(cc)
}

#Function to determine if a "management concern" occurred for each year of a simulation where management concern is defined as a 
#escapement that failed to meet the goal.
#Arguments:
#   S_sim: Run size from the simulation.
#   lb_sim: Lower bound of the escapemtng goal from the simulation.
get_mc <- function(S_sim, lb_sim, ...){  
  mc <- if(S_sim < lb_sim){TRUE} else(FALSE)
  return(mc)
}

#Function to determine if a "yield concern" occurred for each year of a simulation where management concern is defined as a 
#escapement that exceeded the goal.
#Arguments:
#   S_sim: Run size from the simulation.
#   ub_sim: Lower bound of the escapemtng goal from the simulation.
get_yc <- function(S_sim, ub_sim, ...){  
  yc <- if(S_sim > ub_sim){TRUE} else(FALSE)
  return(yc)
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

lna_p = 1 + (0.5 * 0.5 / 2 / (1 - 0.5 * 0.5))
lb_p = lna_p/0.0001*(0.5 - 0.07 * lna_p) * 0.8
ub_p = lna_p/0.0001*(0.5 - 0.07 * lna_p) * 1.6
temp <- simSR_goal(.7, 0.0001, 0.5, 0.5, age0 = Chinook_age, 5000, sigN = 0.4, sigF = 0, Hfun = H_goal, lb_sim = 3904, ub_sim = 7808, power = 0.4)
temp2 <- data.frame(temp$S, temp$N, temp$U, temp$lb, temp$ub, temp$cc, temp$mc, temp$yc, temp$SOC)
apply(temp$N_age, 1, sum)[1:10]
temp$N[1:10]
temp$S[1:10]
temp$R[1:10]

table(temp$cc, temp$mc, temp$yc)
table(temp$cc, temp$SOC)
table(temp$mc, temp$SOC)
table(temp$yc, temp$SOC)
table(temp$cc, temp$mc, temp$yc, temp$SOC)

