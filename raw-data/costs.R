
#' ---
#' title: "bdq amik cost analyis: costs.R"
#'
#' author: "N Green"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     keep_md: TRUE
#' ---

# input created from raw-data/data-clean.R
# using costs per centre (x4) for each component of the pathway
# total costs per patient for each scenario for min, mean, max
# unit costs are calculated and saved as data_input in data/data.RData.
# plots of hospital length of stay by cost are generated.
# [may need to split by self-injector/not self-inj OPAT]


library(dplyr)
library(ggplot2)

load(file = here::here("raw-data", "data_cleaned.RData"))


################
# centre costs #
################
# CENTRE_COSTID: 1.Imperial, 2.St Georges, 3.Birmingham, 4.Royal Free

costs <- list()

costs$raw <-
  data.frame(
    centre_costid =    c(1,2,3,4),
    centre_patientid = c(1,4,2,3),
    # bed = 400, #min:400; max:1200 #per day
    bed_under20d = rep(2429, 4),
    bed_over20d = rep(208, 4),
    PICC =        c(69,   NA,   NA,  200),
    hick =        c(0,    NA,   NA,  200),
    OPAT =        c(102,  148,  100, NA), # 24 hrs
    blood_amak =  c(7.5,  7.5,  NA,  29.56),
    blood_liver = c(3.8,  2.65, NA,  11.82),
    blood_renal = c(3.19, 1.62, NA,  9.45),
    hear =        c(NA,   NA,   33,  34.5),
    ECG =         c(NA,   25,   NA,  15),
    # ECG = c(153, 25, 124, 15), #remove larger values (@Amber)
    amik =        c(2*14.40, 2*11.52, NA, NA)
  )

costs$raw$blood <- costs$raw$blood_liver + costs$raw$blood_renal
costs$raw$lines <- (costs$raw$PICC + costs$raw$hick)/2

centre_scaling <- data.frame(centre_costid = c(1, 2, 3, 4),
                             scaling = c(1.2417, 1.2125, 1.0449, 1.2196))

costs$raw <- merge(costs$raw, centre_scaling, by = "centre_costid")

costs$raw <- mutate(costs$raw,
                    bed_under20d = bed_under20d*scaling,
                    bed_over20d = bed_over20d*scaling)

# cost_bdq_6mo <- 18000 #£
cost_bdq_6mo <- 18700 #£
cost_bdq_8mo <- cost_bdq_6mo/6*8 #£

cost_bdq <- list(
  # 6 months
  rho0.obs =         cost_bdq_6mo,
  rho0.1.obs =       cost_bdq_6mo,
  rho0.33.obs =      cost_bdq_6mo,
  rho0.6mo_days =    cost_bdq_6mo,
  rho0.1.6mo_days =  cost_bdq_6mo,
  rho0.33.6mo_days = cost_bdq_6mo,
  # 8 months
  rho0.8mo_days =    cost_bdq_8mo,
  rho0.1.8mo_days =  cost_bdq_8mo,
  rho0.33.8mo_days = cost_bdq_8mo)

# NHS tariff admission inj
# £2423 for days <19 then add (£208 x each additional day)

## ideal OPAT:
# - monthly hearing tests
# - weekly blood tests
# - daily nurse visit
# - 8 months injectables including admission

max_weeks <- 200

# weekly then fortnightly
# nurse_freq_weeks <- rep(7, max_weeks)
# blood_freq_weeks <- c(2,1,1,1, rep(c(0,1), max_weeks)) # inj only
blood_freq_weeks <- rep(1, max_weeks) # weekly

MAX_WEEKS <- 100
ECG_freq_weeks <- c(1,0,1, rep(c(0,0,0,1), MAX_WEEKS)) #monthly

costs$mean <- replace_NA_cols(costs$raw, mean) %>% round(2)
costs$max <- replace_NA_cols(costs$raw, max) %>% round(2)
costs$min <- replace_NA_cols(costs$raw, min) %>% round(2)

save(costs, file = here::here("data", "costs.RData"))
save(cost_bdq, file = here::here("data", "cost_bdq.RData"))


#################
# patient costs #
#################

costs_sa <- costs
costs_sa$raw <- NULL #remove so dont use in calcs

data_input <- list()

for (i in names(data)) {

  data_input[[i]] <- list()

  for (j in names(costs_sa)) {

    data_joined <- merge(data[[i]], costs_sa[[j]],
                         by = "centre_patientid")

    # amikacin observed, bundled OPAT ---------------------------------------
    data_joined <-
      data_joined %>%
      mutate(
        cost_inj =
          bed_under20d +                                        # fixed hospital cost
          pmax(0, total_admission_days_inj - 19)*bed_over20d +  # variable hospital cost
          num_hickman_line*hick +
          num_picc_line*PICC +
          OPAT_days_inj*OPAT +
          IV_days*amik +
          num_hearing_tests*hear +
          num_ECG_inj*ECG)

    # bedaquiline, oral pill -------------------------------------------------
    data_joined <-
      data_joined %>%
      mutate(
        cost_pill =
          bed_under20d +                                        # fixed hospital cost
          pmax(0, total_admission_days_bdq - 19)*bed_over20d +  # variable hospital cost
          num_blood_bdq*(blood_renal + blood_liver) +
          num_ECG_bdq*ECG +
          cost_bdq[[i]]
        )

    data_input[[i]][[j]] <- select(data_joined,
                                   centre_patientid, cost_inj, cost_pill)
  }
}

save(data_input, file = here::here("data", "data.RData"))


###########
## plots ##
###########

# hospital length of stay vs cost per person
# grid of histograms
# by scenarios

library(ggplot2)
library(purrr)

Tx_cost <- function(Tx_days,
                    costs) {

  adm_days <- 1:Tx_days
  res <-
    costs$raw$bed_under20d[[1]] +
    pmax(0, adm_days - 19)*costs$raw$bed_over20d[[1]] +
    costs$raw$OPAT[[1]]*(Tx_days - adm_days)
  res
}

Tx_cost_ <- map(data$rho0.obs$IV_days, Tx_cost, costs = costs)

png(filename = here::here("plots", "cost_by_los.png"))

plot(NULL,
     xlab = c("Hospital length of stay (days)"),
     ylab = c("Cost per person (£)"),
     xlim = c(0, 230), ylim = c(0, 60000),
     main = "Each patient's total cost by hospital LoS")
map(Tx_cost_, lines, col = alpha(rgb(0,0,0), 0.2), lwd = 2)
# lines(Tx_cost_6mo, col = "red", lwd = 4)
# lines(Tx_cost_8mo, col = "blue", lwd = 3)
abline(v = 19, lty = 2)
axis(side = 1, at = 19, labels = "19")

dev.off()


OP_days_bdq <- map(data, "OP_days_bdq")
par(mfrow = c(3,3))
for (i in seq_along(OP_days_bdq)) {
  hist(OP_days_bdq[[i]]/7, breaks = 50,
       main = names(OP_days_bdq)[i], xlab = "weeks")
}


