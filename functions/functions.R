#Simulate salmon population dynamics with age-at-maturity and management strategy inputs.
#Arguments:
#   lnalpha, beta, sigW, phi: Ricker SR parameters
#   age0: a named vector where the names are numeric total age formatted as a characters(e,g, '3'), and the vector values are proportional age-at-maturity
#   Sims: number of simulations to run
#   sigF: management error
#   Hfun: a function which provides an annual harvest rate based on the management strategy simulated by the function. 
simSR_goal <- function(lnalpha, beta, sigW, phi, age0, Sims0, sigF = 0, Hfun, ...){
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
  
  S <- N <- epsF <- U <- Ft <- H <- Rhat <- lnRhatShat <- fittedR <- rep(NA, Sims + A + a.min + 1)
  
  # recursive portion...
  for (c in (A + a.min):(Sims + A + a.min)) {
    
    epsF[c] <- rnorm(1, 0, sigF)
    N[c] <- sum(N_age[c, ])
    temp <- Hfun(N = N[c], ...)#H_target(N[c])
    U[c] <- temp[[1]]
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
              U = U[(length(U) + 1 - 1 - Sims0):(length(U) - 1)],
              R = R[(length(R) + 1 - Sims0):length(R)],
              lnalpha.y = lnalpha.y[(length(lnalpha.y) + 1 - Sims0):length(lnalpha.y)],
              N_age = N_age[(dim(N_age)[1] + 1 - A - a.min - Sims0):(dim(N_age)[1] - A - a.min), ]
  ))
}


#Function to find U that achieves S at some percentage of bias corrected Smsy
#Arguments:
#   N: run size
#   target0: If a number between 0 and 1 then the percent of Smsy ADF&G will manage for. Defaults to 80% which was Bernard era guidance on escapement goal lower bounds.
#            If a integer greater than 1 then the actual escapement ADF&G will manage for.
#   power: maximum harvest rate for the fishery.
H_target <- function(N, target0 = 0.8, power = .85, ...){
  lna <- with(parent.frame(), lnalpha)
  b <- with(parent.frame(), beta)
  sig <- with(parent.frame(), sigW)
  phi <- with(parent.frame(), phi)
  lnap <- lna + 0.5 * sig * sig / (1 - phi * phi)
  Smsy <- lna / b * (0.5 - 0.07 * lna)
  Smsyp <- lnap / b * (0.5 - 0.07 * lnap)
  
  if(target0 >= 0 & target0 <= 1){target <- Smsyp * target0} else{target = target0}
  
  U <- c(if(N > target){min((N - (target)) / N, power)} else{0})
  
  return(list(U, target = target, power = power))
}


plot_recruit<- function(dat){
  limit <- quantile(dat$R, probs = .9)
  
  ggplot(dat, aes(x = sim, y = R)) +
    geom_line(alpha = 0.2) +
    geom_hline(aes(yintercept = lb_p), linetype = 2) +
    scale_y_continuous(limits = c(0, limit)) +
    scale_x_continuous(name = "Simulation year") +
    facet_grid(paste0("\u03C6: ", phi) ~ paste0("\u03C3: ", sigW),)
}
