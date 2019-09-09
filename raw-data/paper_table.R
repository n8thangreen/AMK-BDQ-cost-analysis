
#' ---
#' title: "bdq-amik cost analysis: paper_table.R"
#'
#' author: "N Green"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     keep_md: TRUE
#' ---

# Input data summary table for paper
#
#                Mean costs (per day/unit)   Observed AK       Ideal AK             0% reduction BDQ   10% reduction BDQ    33% reduction BDQ
# -------------  --------------------------  ----------------  -------------------  -----------------  -------------------  -----------------------
# Admission^     2865.43, 245.37             61 [15.5, 92]     61 [15.5, 92]        61 [15.5, 92]      54.9 [13.95, 82.8]   40.87 [10.385, 63.315]
# Lines          117.25                      1 [1, 1]          1.05 [0.75, 1.505]   0                  0                    0
# ECG            20                          3                 3                    7 [7, 7]           7 [7, 7]             7 [7, 7]
# Hearing test   33.75                       3 [1, 6]          8                    0                  0                    0
# Blood tests    10.84                       0*                0*                   8 [8, 8]           8 [8, 8]             8 [8, 8]
# OPAT           116.67                      91 [0, 149.5]     107 [76, 152.5]      0                  0                    0
# AK (days)      23.04                       160 [91.5, 187]   168                  -                  -                    -
# BDQ (days)     107.14                      -                 -                    168                168                  168


suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(knitr))

load(here::here("raw-data", "data_cleaned.RData"))

# select statistic function
avg_fn <- median

# interquartile range
iqr <- function(x) quantile(as.numeric(x, units = "days"), probs = c(0.25, 0.5, 0.75))
table_row_names <- c("total_admission_days_bdq",
                     "num_lines",
                     "num_ECG_bdq",
                     # "num_ECG_inj",
                     "num_hearing_tests",
                     "num_blood_bdq",
                     "OPAT_days_inj",
                     "IV_days")

# format as median (LQ, UQ)

rho0.obs_stats <- lapply(data$rho0.obs[, table_row_names], FUN = iqr)
rho0.6mo_days_stats <- lapply(data$rho0.6mo_days[, table_row_names], FUN = iqr)
rho0.1.6mo_days_stats <- lapply(data$rho0.1.6mo_days[, table_row_names], FUN = iqr)
rho0.33.6mo_days <- lapply(data$rho0.33.6mo_days[, table_row_names], FUN = iqr)

iqr_format <- function(dat) lapply(dat,
                                   function(x) paste(x[2], " [", x[1], ", ", x[3], "]", sep = ""))

usage_per_patient <- list()
usage_per_patient$rho0.obs_stats <- do.call(rbind, iqr_format(rho0.obs_stats))
usage_per_patient$rho0.6mo_days_stats1 <- do.call(rbind, iqr_format(rho0.6mo_days_stats))
usage_per_patient$rho0.6mo_days_stats2 <- do.call(rbind, iqr_format(rho0.6mo_days_stats))
usage_per_patient$rho0.1.6mo_days_stats <- do.call(rbind, iqr_format(rho0.1.6mo_days_stats))
usage_per_patient$rho0.33.6mo_days <- do.call(rbind, iqr_format(rho0.33.6mo_days))

usage_per_patient_tab <-
  do.call(cbind, usage_per_patient) %>%
  `rownames<-`(NULL)

# substitute assumed values
usage_per_patient_tab[2, c(3,4,5)] <- 0    # lines
usage_per_patient_tab[3, c(1,2)] <- 3      # ECG
usage_per_patient_tab[4, c(3,4,5)] <- 0    # hear
usage_per_patient_tab[5, c(1,2)] <- "0*"   # blood
usage_per_patient_tab[6, c(3,4,5)] <- 0    # OPAT
usage_per_patient_tab[7, c(3,4,5)] <- "-"  # amik

usage_per_patient_tab[7, 2] <- 24*7 # amikacin
usage_per_patient_tab[4, 2] <- 8    # amikacin
usage_per_patient_tab <- rbind(usage_per_patient_tab,
                               c("-", "-", 24*7, 24*7, 24*7)) # bdq

# mean cost across sites
data("costs")

mean_costs <-
  apply(costs$raw, 2, mean, na.rm = TRUE)[
    c("lines", "ECG", "hear", "blood", "OPAT")] %>%
  round(digits = 2)

mean_costs <-
  c("admission" = paste(round(mean(costs$raw$bed_under20d), 2), round(mean(costs$raw$bed_over20d), 2), sep = ", "),
    mean_costs,
    "amik" = 23.04,
    "bdq" = round(18000/(6*4*7), digits = 2))


# tidy labels
names(mean_costs)[names(mean_costs) == "lines"] <- "Lines"
names(mean_costs)[names(mean_costs) == "hear"] <- "Hearing test"
names(mean_costs)[names(mean_costs) == "blood"] <- "Blood tests"
names(mean_costs)[names(mean_costs) == "bdq"] <- "BDQ (days)"
names(mean_costs)[names(mean_costs) == "amik"] <- "AK (days)"
names(mean_costs)[names(mean_costs) == "admission"] <- "Admission^"


combined_table_dat <- cbind(mean_costs,
                            usage_per_patient_tab)

colnames(combined_table_dat) <- c("Mean costs (per day/unit)",
                                  "Observed AK",
                                  "Ideal AK",
                                  "0% reduction BDQ",
                                  "10% reduction BDQ",
                                  "33% reduction BDQ")


suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))

print(
  combined_table_dat %>%
  kable() %>%
  kable_styling("striped") %>%
  add_header_above(c(" " = 1, " " = 1, "Observed treatment" = 1, "24 weeks treatment" = 4)) %>%
  add_header_above(c(" " = 1, " " = 1, "Usage per patient (days/number of tests; median [IQR])" = 5))
)
