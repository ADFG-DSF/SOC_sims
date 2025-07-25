sim_tv <-   sim_SRgamma(c(rep(exp(2), 150), rep(exp(1.5), 150), rep(exp(1), 150), rep(exp(0.46), 150), rep(exp(0.38), 150), rep(exp(0.28), 150)), 
                        beta = 0.0001,
                        gamma = 1,
                        sigW = 0.25,
                        phi = 0,
                        age0 = c('3' = 0.1, '4' = 0.38, '5' = 0.3, '6' = 0.2, '7' = 0.02),
                        Sims0 = 900,
                        Hfun = H_goal,
                        lb_goal = 3893,
                        ub_goal = 10529,
                        power = 0.6,
                        sigF = 0,
                        sigN = 0) %>%
  mutate(lnRS = log(R/S),
         lnalpha = round(log(alpha), 2),
         group = ifelse(lnalpha %in% c(2, 0.46), "2", ifelse(lnalpha %in% c(1.5, 0.38), "1.5", "1"))) #plot the data and the data generating model

colors <- scales::hue_pal()(6)  

plot_simtv <-   
  sim_tv %>%
  ggplot(aes(x = S, y = lnRS, color = as.character(lnalpha))) +
  geom_point() +
  stat_smooth(data = sim_tv[sim_tv$lnalpha %in% c(2, 0.46), ], span = 0.9, 
              aes(x = S, y = lnRS),
              se = FALSE,
              inherit.aes = FALSE,
              color = "black") +
  stat_smooth(data = sim_tv[sim_tv$lnalpha %in% c(1.5, 0.38), ], span = 0.9, 
              aes(x = S, y = lnRS),
              se = FALSE,
              inherit.aes = FALSE,
              color = "black") +
  stat_smooth(data = sim_tv[sim_tv$lnalpha %in% c(1, 0.28), ], span = 0.9, 
              aes(x = S, y = lnRS),
              se = FALSE,
              inherit.aes = FALSE,
              color = "black") +
  theme_bw() +
  annotate("rect", xmin = 3893, xmax = 10529, ymin = -Inf, ymax = Inf, fill = 'gray', alpha = .4) +
  facet_grid(group ~ .)

plot_simtv

