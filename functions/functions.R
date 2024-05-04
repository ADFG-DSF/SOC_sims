Ricker <- function(lnalpha, beta, S)  {S*exp(lnalpha - beta*S)}

#Simulate salmon population dynamics with age-at-maturity and management strategy inputs.
#Arguments:
#   lnalpha, beta, sigW, phi: Ricker SR parameters
#   age0: a named vector where the names are numeric total age formatted as a characters(e,g, '3'), and the vector values are proportional age-at-maturity
#   Sims: number of simulations to run
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
  
  S <- N <- N_hat <- epsF <- U <- lb <- ub <- Ft <- H <- Rhat <- lnRhatShat <- fittedR <- rep(NA, Sims + A + a.min + 1)
  
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
              S = S[(length(S) + 1 - 1 - Sims0):(length(S) - 1)],
              F = Ft[(length(Ft) + 1 - 1 - Sims0):(length(Ft) - 1)],
              U = U[(length(U) + 1 - 1 - Sims0):(length(U) - 1)],
              R = R[(length(R) + 1 - Sims0):length(R)],
              N = N[(length(N) + 1 - Sims0):length(N)],
              N_hat = N_hat[(length(N_hat) + 1 - Sims0):length(N_hat)],
              lnalpha.y = lnalpha.y[(length(lnalpha.y) + 1 - Sims0):length(lnalpha.y)],
              N_age = N_age[(dim(N_age)[1] + 1 - A - a.min - Sims0):(dim(N_age)[1] - A - a.min) + 1, ],
              lb = lb[(length(lb) + 1 - Sims0):length(lb)],
              ub = ub[(length(ub) + 1 - Sims0):length(ub)]
  ))
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
#   goal: If a number between 0 and 1 then the percent of Smsy ADF&G will manage for. Defaults to 80% which was Bernard era guidance on escapement goal lower bounds.
#            If a integer greater than 1 then the actual escapement ADF&G will manage for.
#   power: maximum harvest rate for the fishery.
H_goal <- function(N_sim, lb_sim, ub_sim, power = .85, ...){
  # stopifnot(length(goal) == 2)
  # stopifnot(goal[1] > 1 & goal[2] < (lnalpha + 0.5 * sigW * sigW) / beta & goal[2] > goal[1])
  # lna <- with(parent.frame(), lnalpha)
  # b <- with(parent.frame(), beta)
  # sig <- with(parent.frame(), sigW)
  # phi <- with(parent.frame(), phi)
  # lnap <- lna + 0.5 * sig * sig / (1 - phi * phi)
  # Smsyp <- lnap / b * (0.5 - 0.07 * lnap)
  # goal <- c(0.8 * Smsyp, 1.6 * Smsyp)
  U_lb <- U_ub <- U <- NA

  U_lb <- if(N_sim > lb_sim){(N_sim - lb_sim) / N_sim} else{0} # Find U to fish to lb of goal
  U_ub <- if(N_sim > ub_sim){(N_sim - ub_sim) / N_sim} else{0} # Find U to fish to ub of goal
  U <- min(runif(1, U_ub, U_lb), power)

  return(list(U, lb = lb_sim, ub = ub_sim, power = power))
}
