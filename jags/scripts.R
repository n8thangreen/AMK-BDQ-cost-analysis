
#' ---
#' title: "cost analyis: scripts.R"
#'
#' author: "N Green"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     keep_md: TRUE
#' ---

# using IPD created in raw-data/data-clean.R and raw-data/costs.R
# a Bayesian MCMC is run in jags in parallel over different scenarios.
# See patient_costs_methods.pdf for details of scenarios and how cost
# are calculated.
#
# Plots of mean posterior distributions and cost difference density strips
# are created.
# Posterior predicitive plots are also carried-out for model checking.


library(R2jags)
library(R2WinBUGS)
library(purrr)

source(here::here("jags", "Utils.R"))

# data_input: list by scenario of IPD costs of treatments
load(here::here("data", "data.RData"))

# extract from list means for bdq and amik and create single dataframe
dat <-
  data_input %>%
  map("mean") %>%
  lapply(function(x) x[ ,c("cost_inj", "cost_pill")]) %>%
  do.call(data.frame, args = .)

# transform to log scale for MCMC sampler
ldat <- map_df(dat, function(x) log(as.numeric(x) + 0.001))

stat_name <- "mean" #"max" #"min"
ldat$centre <- data_input[[1]][[stat_name]]$centre_patientid

# copnvert from difftime to numeric
col_difftime <- lapply(dat, function(x) class(x) == "difftime") %>% unlist()
dat[, col_difftime] <- map_df(dat[, col_difftime], as.numeric, unit = "days")

mean_costs <- colMeans(dat)

n_centre <- 4

# all scenarios in single list
# in order to map across
dataJags <-
  list(
    # observed durations
    rho0_obs =
      list(s = n_centre,
           n = nrow(ldat),
           c1 = ldat$rho0.obs.cost_inj,         # patient cost amik
           c2 = ldat$rho0.obs.cost_pill,        # patient cost bdq
           child = ldat$centre),
    rho0.1_obs =
      list(s = n_centre,
           n = nrow(ldat),
           c1 = ldat$rho0.1.obs.cost_inj,       # patient cost amik
           c2 = ldat$rho0.1.obs.cost_pill,      # patient cost bdq
           child = ldat$centre),
    rho0.33_obs =
      list(s = n_centre,
           n = nrow(ldat),
           c1 = ldat$rho0.33.obs.cost_inj,      # patient cost amik
           c2 = ldat$rho0.33.obs.cost_pill,     # patient cost bdq
           child = ldat$centre),
    # 6 months both treatment
    rho0_6mo_days =
      list(s = n_centre,
           n = nrow(ldat),
           c1 = ldat$rho0.6mo_days.cost_inj,    # patient cost amik
           c2 = ldat$rho0.6mo_days.cost_pill,   # patient cost bdq
           child = ldat$centre),
    rho0.1_6mo_days =
      list(s = n_centre,
           n = nrow(ldat),
           c1 = ldat$rho0.1.6mo_days.cost_inj,  # patient cost amik
           c2 = ldat$rho0.1.6mo_days.cost_pill, # patient cost bdq
           child = ldat$centre),
    rho0.33_6mo_days =
      list(s = n_centre,
           n = nrow(ldat),
           c1 = ldat$rho0.33.6mo_days.cost_inj, # patient cost amik
           c2 = ldat$rho0.33.6mo_days.cost_pill,# patient cost bdq
           child = ldat$centre),
    # 6 months bdq vs 8 months amik
    rho0_6mbdq_8mamik =
      list(s = n_centre,
           n = nrow(ldat),
           c1 = ldat$rho0.8mo_days.cost_inj,    # patient cost amik
           c2 = ldat$rho0.6mo_days.cost_pill,   # patient cost bdq
           child = ldat$centre),
    rho0.1_6mbdq_8mamik =
      list(s = n_centre,
           n = nrow(ldat),
           c1 = ldat$rho0.1.8mo_days.cost_inj,  # patient cost amik
           c2 = ldat$rho0.1.6mo_days.cost_pill, # patient cost bdq
           child = ldat$centre),
    rho0.33_6mbdq_8mamik =
      list(s = n_centre,
           n = nrow(ldat),
           c1 = ldat$rho0.33.8mo_days.cost_inj, # patient cost amik
           c2 = ldat$rho0.33.6mo_days.cost_pill,# patient cost bdq
           child = ldat$centre)#,
    ## 8 months both treatment
    # rho0_8mo_days =
    #   list(s = n_centre,
    #        n = nrow(ldat),
    #        c1 = ldat$rho0.8mo_days.cost_inj,
    #        c2 = ldat$rho0.8mo_days.cost_pill,
    #        child = ldat$centre),
    # rho0.1_8mo_days =
    #   list(s = n_centre,
    #        n = nrow(ldat),
    #        c1 = ldat$rho0.1.8mo_days.cost_inj,
    #        c2 = ldat$rho0.1.8mo_days.cost_pill,
    #        child = ldat$centre),
    # rho0.33_8mo_days =
    #   list(s = n_centre,
    #        n = nrow(ldat),
    #        c1 = ldat$rho0.33.8mo_days.cost_inj,
    #        c2 = ldat$rho0.33.8mo_days.cost_pill,
    #        child = ldat$centre)
  )

# to monitor
params <-
  c("mu.c", "sigma.c",
    "mu.centre", "sigma.centre",
    # "mu.c_nat", "sigma.c_nat",
    # "mu.centre_nat", "sigma.centre_nat",
    "c1.rep", "c2.rep",
    "c1.ppost", "c2.ppost"
    # "delta_c",
    # "thresh"
  )

# starting values
inits <- function(){
  list(
    logsigma.c = c(5, 5), #runif(2, 500, 1500),
    logsigma.centre = c(5, 5), #runif(2, 500, 1500),
    mu.centre = c(10, 10)
  )
}

#test
# n_iter <- 1000
# n_burnin <- 10
# n_thin <- 1

n_iter <- 1e6
n_burnin <- 1e3
n_thin <- 1e2 #floor((n_iter - n_burnin)/500)

# jags_partial <-
#   partial(jags,
#           # inits = inits,
#           parameters.to.save = params,
#           model.file = "jags2.txt",
#           n.chains = 2,
#           n.iter = n_iter,
#           n.burnin = n_burnin,
#           n.thin = n_thin,
#           DIC = TRUE,
#           working.directory = here::here("jags"),
#           progress.bar = "text")


##############
## run MCMC ##
##############

library(parallel)

no_cores <- detectCores() - 1

# initiate cluster
cl <- makeCluster(no_cores)#, outfile = "temp_logfile.txt")

# clusterEvalQ(cl, library(data.tree))
# clusterEvalQ(cl, source("sample_subset_pop_dectree.R"))

set.seed(12345)

ptm <- proc.time()

out <- parLapplyLB(cl,
                   dataJags,
                   fun = jags,
                   inits = list(inits(), inits()),
                   parameters.to.save = params,
                   model.file = "jags2.txt",
                   n.chains = 2,
                   n.iter = n_iter,
                   n.burnin = n_burnin,
                   n.thin = n_thin,
                   DIC = TRUE,
                   working.directory = here::here("jags"),
                   progress.bar = "text")

ptm_elapsed <- (proc.time() - ptm)/60
ptm_elapsed

stopCluster(cl)

#serial
# out <- map(dataJags, jags_partial)

# lapply(out, print,
#        # digits = 3,
#        intervals = c(0.025, 0.975))

BUGSoutput <- map(out, "BUGSoutput")

save(BUGSoutput, file = here::here("jags", "data", "BUGSoutput.RData"))


# tranform from log scale back to natural
inv_lnorm_transform <- function(mu, sigma)
  exp(mu + 0.5*(sigma)^2)


#########
# plots #
#########


library(reshape2)
library(ggplot2)
library(mcmcplots)

summary_table <- list()
delta_c_all <- list()
m.centre_long <- NULL

x11()
par(mfcol = c(3,2))

for (j in seq_along(BUGSoutput)) {

  attach.bugs(BUGSoutput[[j]])

  ## health economic analysis

  m.centre <- matrix(NA, nrow = n.sims, ncol = 2)
  m.c <- array(NA, dim = c(n.sims, n_centre, 2))

  for (t in 1:2) {
    m.centre[ ,t] <- inv_lnorm_transform(mu.centre[ ,t],
                                         sigma.centre[t])

    for (i in seq_len(n_centre)) {
      m.c[ ,i,t] <- inv_lnorm_transform(mu.c[ ,i,t],
                                        sigma.c[t])
    }
  }

  delta_c_all[[j]] <- m.centre[ ,2] - m.centre[ ,1]

  # drop outliers
  # m.centre[ ,1][m.centre[ ,1] > 100000] <- NA
  # m.centre[ ,2][m.centre[ ,1] > 100000] <- NA

  summary_table[[j]] <-
    rbind(
      do.call(data.frame,
              list(mean = apply(m.centre, 2, mean, na.rm = TRUE),
                   sd = apply(m.centre, 2, sd, na.rm = TRUE),
                   median = apply(m.centre, 2, median, na.rm = TRUE),
                   min = apply(m.centre, 2, min, na.rm = TRUE),
                   max = apply(m.centre, 2, max, na.rm = TRUE),
                   n = apply(m.centre, 2, length))),
      do.call(data.frame,
              list(mean = apply(m.c[,,1], 2, mean, na.rm = TRUE),
                   sd = apply(m.c[,,1], 2, sd, na.rm = TRUE),
                   median = apply(m.c[,,1], 2, median, na.rm = TRUE),
                   min = apply(m.c[,,1], 2, min, na.rm = TRUE),
                   max = apply(m.c[,,1], 2, max, na.rm = TRUE),
                   n = apply(m.c[,,1], 2, length))),
      do.call(data.frame,
              list(mean = apply(m.c[,,2], 2, mean, na.rm = TRUE),
                   sd = apply(m.c[,,2], 2, sd, na.rm = TRUE),
                   median = apply(m.c[,,2], 2, median, na.rm = TRUE),
                   min = apply(m.c[,,2], 2, min, na.rm = TRUE),
                   max = apply(m.c[,,2], 2, max, na.rm = TRUE),
                   n = apply(m.c[,,2], 2, length)))
      # do.call(data.frame,
      #         list(mean = mean(BUGSoutput[[j]]$sims.list[["thresh"]]),
      #              sd = sd(BUGSoutput[[j]]$sims.list[["thresh"]]),
      #              median = median(BUGSoutput[[j]]$sims.list[["thresh"]]),
      #              min = min(BUGSoutput[[j]]$sims.list[["thresh"]]),
      #              max = max(BUGSoutput[[j]]$sims.list[["thresh"]]),
      #              n = length(BUGSoutput[[j]]$sims.list[["thresh"]]))),
      # do.call(data.frame,
      #         list(mean = mean(BUGSoutput[[j]]$sims.list[["delta_c"]]),
      #              sd = sd(BUGSoutput[[j]]$sims.list[["delta_c"]]),
      #              median = median(BUGSoutput[[j]]$sims.list[["delta_c"]]),
      #              min = min(BUGSoutput[[j]]$sims.list[["delta_c"]]),
      #              max = max(BUGSoutput[[j]]$sims.list[["delta_c"]]),
      #              n = length(BUGSoutput[[j]]$sims.list[["delta_c"]])))
    )

  summary_table[[j]] <- data.frame("name" = c("m.centre1", "m.centre2",
                                              "m.c1.1", "m.c1.2", "m.c1.3", "m.c1.4",
                                              "m.c2.1", "m.c2.2", "m.c2.3", "m.c2.4"
                                              #"thresh", "delta_c"
  ),
  summary_table[[j]])

  m.centre_long <- rbind(m.centre_long,
                         cbind(j, m.centre))

  # between centres ---

  hist(m.centre[ ,1],
       main = names(dataJags)[j],
       xlab = "Cost (£)",
       # xlim = c(10000, 50000),
       xlim = c(10000, 50000),
       # ylim = c(0,250),
       ylim = c(0, 0.0003),
       breaks = 100,
       freq = FALSE)
  hist(m.centre[ ,2],
       add = TRUE,
       col = "red",
       breaks = 100,
       freq = FALSE)
  abline(v = mean_costs[2*j - 1])
  abline(v = mean_costs[2*j], col = "red")

  # raw data
  # hist(dat[ ,2*j - 1],
  #      add = TRUE,
  #      col = "green",
  #      breaks = 100,
  #      freq = FALSE)
  # hist(dat[ ,2*j],
  #      add = TRUE,
  #      col = "blue",
  #      breaks = 100,
  #      freq = FALSE)


  mc_long1 <- melt(m.c[,,1])
  mc_long1$Var2 <- as.factor(mc_long1$Var2)
  mc_long2 <- melt(m.c[,,2])
  mc_long2$Var2 <- as.factor(mc_long2$Var2)


  # within centre ---

  ggplot(mc_long1, aes(x = value,
                       # fill = Var2,
                       color = Var2)) +
    theme_bw() +
    geom_density(alpha = 0.6) #+
  # geom_histogram(
  #   # fill = "white",
  #   alpha = 0.5,
  #   position = "identity")

}

names(summary_table) <- names(dataJags)

save(summary_table, file = here::here("jags", "data", "summary_table.RData"))
save(out, file = here::here("jags", "data", "jags_output.Rdata"))


###########
# ggplots #
###########

colnames(m.centre_long)[2] <- "amik"
colnames(m.centre_long)[3] <- "bdq"

##TODO##
## fudge: repeat sims that should have the same value anyway
# m.center_temp <- m.centre_long
# j: 1) rho0_obs
#    2) rho0.1_obs
#    3) rho0.33_obs
#    4) rho0_6mo_days
#    5) rho0.1_6mo_days
#    6) rho0.33_6mo_days
#    7) rho0_6mbdq_8mamik
#    8) rho0.1_6mbdq_8mamik
#    9) rho0.33_6mbdq_8mamik

m.centre_long[m.centre_long[ ,"j"] == 2, "amik"] <- m.centre_long[m.centre_long[ ,"j"] == 1, "amik"]
m.centre_long[m.centre_long[ ,"j"] == 3, "amik"] <- m.centre_long[m.centre_long[ ,"j"] == 1, "amik"]
m.centre_long[m.centre_long[ ,"j"] == 5, "amik"] <- m.centre_long[m.centre_long[ ,"j"] == 4, "amik"]
m.centre_long[m.centre_long[ ,"j"] == 6, "amik"] <- m.centre_long[m.centre_long[ ,"j"] == 4, "amik"]
m.centre_long[m.centre_long[ ,"j"] == 4, "bdq"] <-  m.centre_long[m.centre_long[ ,"j"] == 1, "bdq"]
m.centre_long[m.centre_long[ ,"j"] == 2, "bdq"] <-  m.centre_long[m.centre_long[ ,"j"] == 5, "bdq"]
m.centre_long[m.centre_long[ ,"j"] == 3, "bdq"] <-  m.centre_long[m.centre_long[ ,"j"] == 6, "bdq"]

m.centre_long[m.centre_long[ ,"j"] == 7, "bdq"] <-  m.centre_long[m.centre_long[ ,"j"] == 4, "bdq"]
m.centre_long[m.centre_long[ ,"j"] == 8, "bdq"] <-  m.centre_long[m.centre_long[ ,"j"] == 5, "bdq"]
m.centre_long[m.centre_long[ ,"j"] == 9, "bdq"] <-  m.centre_long[m.centre_long[ ,"j"] == 6, "bdq"]
m.centre_long[m.centre_long[ ,"j"] == 7, "amik"] <-  m.centre_long[m.centre_long[ ,"j"] == 8, "amik"]
m.centre_long[m.centre_long[ ,"j"] == 9, "amik"] <-  m.centre_long[m.centre_long[ ,"j"] == 8, "amik"]

# trim outliers
MAX_COST <- 50000
m.centre_long <- m.centre_long[m.centre_long[, "amik"] < MAX_COST, ]
m.centre_long <- m.centre_long[m.centre_long[, "bdq"] < MAX_COST, ]



# from [j |amik| bda| delta] to [j |variable| value]
# where variable is {amik,bdq,delta}
# then append Tx an rho
m.centre_long2 <- reshape2::melt(data.frame(m.centre_long), id.vars = "j")
m.centre_long2 <- merge(m.centre_long2,
                        data.frame(j = 1:9,
                                   Tx = rep(c("obs", "6 months", "8 months"), each = 3),
                                   rho = c(0, 0.1, 0.33)),
                        by = "j")

hist_dat <- m.centre_long2[m.centre_long2$Tx != "8 months", ]



# histogram
p <-
  ggplot(hist_dat, aes(x = value, fill = variable)) +
  geom_histogram(position = "identity", binwidth = 50, alpha = 0.7) +
  # geom_freqpoly(position = "identity", binwidth = 500, alpha = 0.7) +
  facet_grid(rho ~ Tx, labeller = labeller(Tx = as_labeller(c("6 months" = "6 months of injectables",
                                                              # "8 months" = "8 months",
                                                              "obs" =   "Observed duration of injectables")),
                                           rho = as_labeller(c("0" =    "No change in duration",
                                                               "0.1" =  "10% reduction in admission",
                                                               "0.33" = "33% reduction in admission")))) +
  xlim(0, 50000) +
  ylim(0, 500) +
  # ylim(0, 1000) +
  theme_bw() +
  xlab("Posterior mean cost per patient (£)") +
  labs(fill = "Treatment") +
  scale_fill_discrete(labels = c("Amikacin", "Bedaquiline")) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        text = element_text(size = 20),
        plot.margin = margin(2,2,2,2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11))
p

p_bw <- p + scale_color_grey() + scale_fill_grey(labels = c("Amikacin", "Bedaquiline"))
p_bw

# export
ggsave(p,
       filename = here::here("plots", "cost-histograms.tiff"),
       width = 25, height = 25, units = "cm", dpi = 300)
ggsave(p_bw,
       filename = here::here("plots", "cost-histograms_bw.tiff"),
       width = 25, height = 25, units = "cm", dpi = 300)

save(m.centre_long2, file = "data/m_centre_long2.RData")


# density strip plot -----------------------------------------------------

library(denstrip)

m.centre_long <- cbind(m.centre_long,
                       delta = m.centre_long[ ,"bdq"] - m.centre_long[ ,"amik"])

delta_c_all <- split(x = m.centre_long[ ,"delta"],
                     f = m.centre_long[ , "j"])

# remove outliers so histogram bins are not too big
# otherwise looks blocky
delta_c_all[[1]] <- delta_c_all[[1]][delta_c_all[[1]] < 40000]
delta_c_all[[2]] <- delta_c_all[[2]][delta_c_all[[2]] < 40000]
delta_c_all[[3]] <- delta_c_all[[3]][delta_c_all[[3]] < 40000]
delta_c_all[[4]] <- delta_c_all[[4]][delta_c_all[[4]] < 40000]
delta_c_all[[5]] <- delta_c_all[[5]][delta_c_all[[5]] < 40000]
delta_c_all[[6]] <- delta_c_all[[6]][delta_c_all[[6]] < 40000]
delta_c_all[[7]] <- delta_c_all[[7]][delta_c_all[[7]] < 40000]
delta_c_all[[8]] <- delta_c_all[[8]][delta_c_all[[8]] < 40000]
delta_c_all[[9]] <- delta_c_all[[9]][delta_c_all[[9]] < 40000]

delta_c_all[[1]] <- delta_c_all[[1]][delta_c_all[[1]] > -40000]
delta_c_all[[2]] <- delta_c_all[[2]][delta_c_all[[2]] > -40000]
delta_c_all[[3]] <- delta_c_all[[3]][delta_c_all[[3]] > -40000]
delta_c_all[[4]] <- delta_c_all[[4]][delta_c_all[[4]] > -40000]
delta_c_all[[5]] <- delta_c_all[[5]][delta_c_all[[5]] > -40000]
delta_c_all[[6]] <- delta_c_all[[6]][delta_c_all[[6]] > -40000]
delta_c_all[[7]] <- delta_c_all[[7]][delta_c_all[[7]] > -40000]
delta_c_all[[8]] <- delta_c_all[[8]][delta_c_all[[8]] > -40000]
delta_c_all[[9]] <- delta_c_all[[9]][delta_c_all[[9]] > -40000]

denstrip_dat <- do.call(cbind, delta_c_all[1:6])
denstrip_all <- do.call(cbind, delta_c_all[1:9])
# denstrip_dat <- denstrip_all

x11(width = 14)
# png(filename = here::here("plots", "density_strips.png"),
#     width = 1000, height = 700,
#     units = "px", res = 300)
# tiff(filename = here::here("plots", "density_strips.tiff"),

LABELS <- c(
  "Observed amikacin vs bedaquiline; equal length of stay",
  "Observed amikacin vs bedaquiline; 10% reduction in admission",
  "Observed amikacin vs bedaquiline; 33% reduction in admission",
  "6 months amikacin vs bedaquiline; equal length of stay",
  "6 months amikacin vs bedaquiline; 10% reduction in admission",
  "6 months amikacin vs bedaquiline; 33% reduction in admission"#,
  # "8 months amikacin and 6 months bedaquiline durations; equal length of stay",
  # "8 months amikacin and 6 months bedaquiline durations; 90% length of stay",
  # "8 months amikacin and 6 months bedaquiline durations; 66% length of stay"
)

tiff(filename = here::here("plots", "density_strips_bw.tiff"),
     width = 3000, height = 2000,
     units = "px", res = 300)

par(mai = c(1, 5, 1, 0.5))
# caterplot(denstrip_dat,
mycaterplot(denstrip_dat,
            bty = "n",
            quantiles = list(outer = c(0.025, 0.975), inner = c(0.16, 0.84)),
            denstrip = TRUE,
            val.lim = c(-20000, 30000),
            reorder = FALSE,
            # labels = names(dataJags),
            labels = LABELS,
            style = "plain",
            col = "grey",
            cex.labels = 1)
abline(v = 0, lty = 2, col = "grey")
mtext(bquote("Mean total cost difference," ~ c[bdq] - c[amik] ~ "(£)"), side = 1, line = 3)
arrows(5000, 6.5, 15000, 6.5, length = 0.08)
arrows(-5000, 6.5, -15000, 6.5, length = 0.08)
text(x = 10000, y = 7, "Amikacin cheaper", )
text(x = -10000, y = 7, "Bedaquiline cheaper")

dev.off()


save(denstrip_dat, file = "data/denstrip_dat.RData")
save(denstrip_all, file = "data/denstrip_all.RData")


#################################
## posterior predictive checks ##
#################################

library(gridExtra)
library(bayesplot)
# library(cowplot)
# library(lattice)
# library(grid)
# library(ggthemes)

ppc_dens_list <- list()

for (i in 1:9) {

  ppc_dens_list$c1[[i]] <-
    ppc_dens_overlay(y = exp(dataJags[[i]]$c1),
                     yrep = exp(out[[i]]$BUGSoutput$sims.list$c1.rep[1:50, ])) +
    xlim(0, 2e5)

  ppc_dens_list$c2[[i]] <-
    ppc_dens_overlay(y = exp(dataJags[[i]]$c2),
                     yrep = exp(out[[i]]$BUGSoutput$sims.list$c2.rep[1:50, ])) +
    xlim(0, 2e5)
}


theme_grid <- theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
                    plot.margin = margin(2,2,2,2, "cm"))
margin <- theme(plot.margin = unit(c(0,0,0,0), "cm"))

ppc_c1_grid <-
  grid.arrange(
    ppc_dens_list$c1[[1]] + labs(subtitle = "rho0_obs") + theme_grid + margin,
    ppc_dens_list$c1[[2]] + labs(subtitle = "rho0.1_obs") + theme_grid + margin,
    ppc_dens_list$c1[[3]] + labs(subtitle = "rho0.33_obs") + theme_grid + margin,
    ppc_dens_list$c1[[4]] + labs(subtitle = "rho0_6mo_days") + theme_grid + margin,
    ppc_dens_list$c1[[5]] + labs(subtitle = "rho0.1_6mo_days") + theme_grid + margin,
    ppc_dens_list$c1[[6]] + labs(subtitle = "rho0.33_6mo_days") + theme_grid + margin,
    ppc_dens_list$c1[[7]] + labs(subtitle = "rho0_8mo_days") + theme_grid + margin,
    ppc_dens_list$c1[[8]] + labs(subtitle = "rho0.1_8mo_days") + theme_grid + margin,
    ppc_dens_list$c1[[9]] + labs(subtitle = "rho0.33_8mo_days") + theme_grid + margin,
    ncol = 3)

ppc_c2_grid <-
  grid.arrange(
    ppc_dens_list$c2[[1]] + labs(subtitle = "rho0_obs") + theme_grid + margin,
    ppc_dens_list$c2[[2]] + labs(subtitle = "rho0.1_obs") + theme_grid + margin,
    ppc_dens_list$c2[[3]] + labs(subtitle = "rho0.33_obs") + theme_grid + margin,
    ppc_dens_list$c2[[4]] + labs(subtitle = "rho0_6mo_days") + theme_grid + margin,
    ppc_dens_list$c2[[5]] + labs(subtitle = "rho0.1_6mo_days") + theme_grid + margin,
    ppc_dens_list$c2[[6]] + labs(subtitle = "rho0.33_6mo_days") + theme_grid + margin,
    ppc_dens_list$c2[[7]] + labs(subtitle = "rho0_8mo_days") + theme_grid + margin,
    ppc_dens_list$c2[[8]] + labs(subtitle = "rho0.1_8mo_days") + theme_grid + margin,
    ppc_dens_list$c2[[9]] + labs(subtitle = "rho0.33_8mo_days") + theme_grid + margin,
    ncol = 3)

# export
ggsave(ppc_c1_grid,
       filename = here::here("plots", "ppc_c1_grid.png"),
       width = 20, height = 20, units = "cm")

ggsave(ppc_c2_grid,
       filename = here::here("plots", "ppc_c2_grid.png"),
       width = 20, height = 20, units = "cm")

## by centre
# i <- 1L
# ppc_stat_grouped(y = exp(dataJags[[i]]$c1),
#                  yrep = exp(out[[i]]$BUGSoutput$sims.list$c1.rep),
#                  group = dataJags[[i]]$child,
#                  stat = "median") + xlim(0, 1e5)
#
# ppc_stat_grouped(y = exp(dataJags[[i]]$c2),
#                  yrep = exp(out[[i]]$BUGSoutput$sims.list$c2.rep),
#                  group = dataJags[[i]]$child,
#                  stat = "median") + xlim(0, 1e5)


xx <-
  m.centre_long2 %>%
  group_by(j, variable) %>%
  summarise(mean_cost = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            L95 = quantile(value, probs = 0.025),
            U95 = quantile(value, probs = 0.975))
