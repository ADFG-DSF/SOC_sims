# # Functions ---------------------------------------------------------------
# Ricker <- function(S, lnalpha, beta) {
#   S*exp(lnalpha - beta*S)
# }
# 
# #Simulate a single stock
# # goal: a vector of length 1 (if providing a target harvest rate) or 2 (if providing an escapement goal range).
# # sigS: measurement error
# # phi: residual autocorrelation
# # power: fishing power
# # sigF: implementation error
# simulateSR_goal <- function(lnalpha, beta, sigW, N, goal, sigS = 0, phi = 0, power = 1, sigF = 0) {
#   stopifnot(length(goal) %in% 1:2)
#   if(length(goal) == 1){stopifnot(goal >= 0 & goal < 1)}
#   if(length(goal) == 2){stopifnot(goal[1] > 1 & goal[2] < (lnalpha + 0.5 * sigW * sigW) / beta & goal[2] > goal[1])}
#   stopifnot(power >= 0 & power <=1) # fishing power is 0 for no harvest and 1 when the fishery can fully control escapement
#   
#   # ----- initial values ----- #
#   # initial value for S: Seq minus some harvest
#   S <- 900
#   
#   # initial value for observed S
#   Shat <- S*rlnorm(1, sdlog=sigS)
#   
#   # initializing all other values
#   redresid <- 0
#   E1R <- E2R <- whiteresid <- epsF <- U1 <- U2 <- U <- Ft <- H <- Rhat <- lnRhatShat <- fittedR <- NA
#   R <- 1 #This value in never used but needed to evaluate U
#   
#   # recursive portion...
#   for(i in 2:(N + 50 +1)) { #added 50 (will discard later)
#     E1R[i] <- S[i-1]*exp(lnalpha - beta*S[i-1])
#     E2R[i] <- E1R[i]*exp(phi*redresid[i-1])
#     R[i] <- E2R[i]*rlnorm(1,0,sigW)
#     redresid[i] <- log(R[i]/E1R[i])
#     whiteresid[i] <- log(R[i]/E2R[i])
#     epsF[i] <- rnorm(1,0,sigF)
#     if(length(goal) == 2){
#       U1[i] <- if(R[i] > goal[1]){(R[i] - goal[1]) / R[i]} else{0} # Find U to fish to lb of goal
#       U2[i] <- if(R[i] > goal[2]){(R[i] - goal[2]) / R[i]} else{0} # Find U to fish to ub of goal
#     }
#     U[i] <- if(length(goal) == 1){goal} else{min(runif(1, U2[i], U1[i]), power)} # draw a harvest rate
#     Ft[i] <- -log(1 - U[i])*exp(epsF[i])
#     S[i] <- R[i]*exp(-Ft[i])
#     Shat[i] <- S[i]*rlnorm(1, sdlog=sigS)
#     H[i] <- R[i]-S[i]
#     Rhat[i] <- Shat[i]+H[i]
#     lnRhatShat[i] <- log(Rhat[i]/Shat[i])
#   }
#   
#   return(list(S=Shat[1:N + 50],
#               R=Rhat[2:(N+1) + 50],
#               Strue=S[1:N + 50],
#               Rtrue=R[2:(N+1) + 50],
#               yield=R[2:(N+1) + 50] - S[1:N + 50],
#               goal = goal,
#               lnalpha = lnalpha, 
#               beta = beta,
#               sigW = sigW))
# }
# 
# #Simulate 2 stocks w same dynamics but fished to different goals.
# simSR_goal <- function(lnalpha, beta, sigW, phi, N, sigF = 0, Hfun, ...){
#   # ----- initial values ----- #
#   # initial value for S: Seq minus some harvest
#   S <- c(900, rep(NA, N + 50))
#   
#   # initializing all other values
#   redresid <- 0
#   E1R <- E2R <- whiteresid <- epsF <- U <- Ft <- H <- Rhat <- lnRhatShat <- fittedR <-
#     rep(NA, N + 51)
#   R <- c(1, rep(NA, N + 50))
#   
#   # recursive portion...
#   for(i in 2:(N + 50 + 1)) { #added 50 (will discard later)
#     E1R[i] <- S[i-1]*exp(lnalpha - beta*S[i-1]) 
#     E2R[i] <- E1R[i]*exp(phi*redresid[i-1])
#     R[i] <- E2R[i]*rlnorm(1,0,sigW)
#     redresid[i] <- log(R[i]/E1R[i])
#     whiteresid[i] <- log(R[i]/E2R[i])
#     epsF[i] <- rnorm(1,0,sigF)
#     temp <- Hfun(R = R[i], ...)
#     U[i] <- temp[[1]]
#     Ft[i] <- -log(1 - U[i])*exp(epsF[i])
#     S[i] <- R[i]*exp(-Ft[i])
#     H[i] <- R[i]-S[i]
#   }
#   return(data.frame(sim = 1:N,
#                     S = S[1:N + 50],
#                     R = R[2:(N+1) + 50],
#                     Y = R[2:(N+1) + 50] - S[1:N + 50],
#                     U = U[2:(N+1) +50],
#                     power = temp[["power"]],
#                     lnalpha = lnalpha,
#                     beta = beta,
#                     sigW = sigW))
# }
# 
# #Simulate 2 stocks w same dynamics but fished to different goals.
# simSR_goal <- function(lnalpha, beta, sigW, phi, age0, N, sigF = 0, Hfun, ...){
#   age <- matrix(age0, N, length(age0), byrow = TRUE)
#   A <- dim(age)[2]
#   a.min <- names(age)[1]
#   a.max <- a.min + A
#   lnalpha.p  <- lnalpha + 0.5 * sigW * sigW / (1-phi*phi) 
#   Seq <- lnalpha.p / beta
#   
#   # ----- initial values ----- #
#   for (c in 1:a.max) { 
#     R[c] <- Ricker(lnalpha.p, beta, Seq) #Ricker mean at Seq
#   }
# 
#   # initializing all other values
#   redresid <- 0
#   E1R <- E2R <- whiteresid <- epsF <- U <- Ft <- H <- Rhat <- lnRhatShat <- fittedR <- rep(NA, N + 100 +a.min + 1)
#   
#   # recursive portion...
#   for (c in (A + a.min):(N + A - 1 + 100)) { #added 100 (will discard later)
#     E1R[c] <- S[c - 1]*exp(lnalpha - beta*S[c - 1]) 
#     E2R[c] <- E1R[c]*exp(phi*redresid[c - 1])
#     R[c] <- E2R[c]*rlnorm(1, 0, sigW)
#     redresid[c] <- log(R[c] / E1R[c])
#     whiteresid[c] <- log(R[c] / E2R[c])
#   }
#   
#   for (a in 1:A) {
#     for (c in a:(N + (a - 1))) {
#       N[c - (a - 1), (A + 1 - a)] <- p[c, (A + 1 - a)] * R[c]
#     }
#   }
#   
#   for (y in 1:N) {
#     epsF[y] <- rnorm(1, 0, sigF)
#     temp <- Hfun(N = N[y], ...)
#     U[y] <- temp[[1]]
#     Ft[y] <- -log(1 - U[y])*exp(epsF[y])
#     S[y] <- sum(N[y, ]) * exp(-Ft[y])
#     H[y] <- R[y] - S[y]
#   }
#   
#   return(data.frame(sim = 1:N,
#                     S = S,
#                     R = R,
#                     N = N))
#   
#   # return(data.frame(sim = 1:N,
#   #                   S = S[1:N + 50],
#   #                   R = R[2:(N+1) + 50],
#   #                   Y = R[2:(N+1) + 50] - S[1:N + 50],
#   #                   U = U[2:(N+1) +50],
#   #                   power = temp[["power"]],
#   #                   lnalpha = lnalpha,
#   #                   beta = beta,
#   #                   sigW = sigW))
# }
# 
# 
# #Function to find U that archives S somewhere within the goal
# H_goal <- function(R, goal, power = .85, ...){
#   # stopifnot(length(goal) == 4)
#   # stopifnot(goal[1, 1] > 1 & goal[1, 2] < (lnalpha + 0.5 * sigW * sigW) / beta & goal[1, 2] > goal[1, 1])
#   # stopifnot(goal[2, 1] > 1 & goal[2, 2] < (lnalpha + 0.5 * sigW * sigW) / beta& goal[2, 2] > goal[2, 1])
#   U_lb <- U_ub <- U <- c(NA, NA)
#   
#   for(j in 1:2){
#     U_lb[j] <- if(R[j] > goal[j, 1]){(R[j] - goal[j, 1]) / R[j]} else{0} # Find U to fish to lb of goal
#     U_ub[j] <- if(R[j] > goal[j, 2]){(R[j] - goal[j, 2]) / R[j]} else{0} # Find U to fish to ub of goal
#     U[j] <- min(runif(1, U_ub[j], U_lb[j]), power)
#   }
#   return(list(U, goal = goal, power = power))
# }
# 
# #Function to find U that archives Smsy for each goal
# H_Smsy <- function(R, power = 0.85){
#   lna <- with(parent.frame(), lnalpha)
#   b <- with(parent.frame(), beta)
#   sig <- with(parent.frame(), sigW)
#   lnap <- lna + 0.5 * sig * sig
#   Smsy <- lna / b * (0.5 - 0.07 * lna)
#   Smsyp <- lnap / b * (0.5 - 0.07 * lnap)
#   
#   U <- c(if(R > Smsyp){min((R - Smsyp) / R, power)} else{0})
#   
#   return(list(U, power = power))
# }
# 
# #Function to find U that archives S a some percentage of the goal range
# H_target <- function(R, goal, target, power = .85, ...){
#   S_target <- U <- c(NA, NA)
#   S_target <- c(goal[1,1] + target * (goal[1,2] - goal[1,1]), 
#                 goal[2,1] + target * (goal[2,2] - goal[2,1]))
#   
#   U <- c(if(R[1] > S_target[1]){min((R[1] - S_target[1]) / R[1], power)} else{0},
#          if(R[2] > S_target[2]){min((R[2] - S_target[2]) / R[2], power)} else{0})
#   
#   return(list(U, goal = goal, power = power))
# }
# 
# #Simulate 2 stocks w same dynamics but fished to different goals.
# simSR_goals <- function(lnalpha, beta, sigW, N, sigF = 0, Hfun, ...){
#   # ----- initial values ----- #
#   # initial value for S: Seq minus some harvest
#   S <- matrix(c(900, rep(NA, N + 50), 900, rep(NA, N + 50)), N + 51, 2)
#   
#   # initializing all other values
#   redresid <- 0
#   E1R <- E2R <- redresid <- whiteresid <- epsF <- U <- Ft <- H <- Rhat <- lnRhatShat <- fittedR <-
#     matrix(c(rep(NA, N + 51), rep(NA, N + 51)), N + 51, 2)
#   R <-  matrix(c(1, rep(NA, N + 50), 1, rep(NA, N + 50)), N + 51, 2)
#   
#   # recursive portion...
#   for(i in 2:(N + 50 + 1)) { #added 50 (will discard later)
#     E1R[i, ] <- S[(i-1), ]*exp(lnalpha - beta*S[(i-1), ]) 
#     R[i, ] <- E1R[i, ]*rlnorm(1,0,sigW)
#     epsF[i] <- rnorm(1,0,sigF)
#     temp <- Hfun(R = R[i, ], ...)
#     U[i,] <- temp[[1]]
#     Ft[i, ] <- -log(1 - U[i, ])*exp(epsF[i])
#     S[i, ] <- R[i, ]*exp(-Ft[i, ])
#     H[i, ] <- R[i, ]-S[i, ]
#   }
#   return(data.frame(sim = rep(1:N, times = 2),
#                     goal = rep(c("mode", "mean"), each = N),
#                     S = c(S[1:N + 50, 1], S[1:N + 50, 2]),
#                     R = c(R[2:(N+1) + 50, 1], R[2:(N+1) + 50, 2]),
#                     Y = c(R[2:(N+1) + 50, 1] - S[1:N + 50, 1], R[2:(N+1) + 50, 2] - S[1:N + 50, 2]),
#                     U = c(U[2:(N+1) +50, 1], U[2:(N+1) +50, 2]),
#                     power = temp[["power"]],
#                     lnalpha = lnalpha,
#                     beta = beta,
#                     sigW = sigW))
# }
# 
# #table of realized yield near Smsy... Not sure I need this
# tab_Smsy <- function(list){
#   lnalpha_c <- list$lnalpha + .5 * list$sigW * list$sigW
#   data.frame(median = median(list$yield[list$S < list$lnalpha / list$beta * (0.5 - 0.07 * list$lnalpha) + 50 &
#                                           list$S > list$lnalpha / list$beta * (0.5 - 0.07 * list$lnalpha) - 50]),
#              mean = mean(list$yield[list$S < lnalpha_c / list$beta * (0.5 - 0.07 * lnalpha_c) + 10 & 
#                                       list$S > lnalpha_c / list$beta * (0.5 - 0.07 * lnalpha_c) - 10]))
# }