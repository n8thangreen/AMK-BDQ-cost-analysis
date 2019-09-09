
#' ---
#' title: "bdq amik cost analyis: data-clean.R"
#'
#' author: "N Green"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     keep_md: TRUE
#' ---

# Takes the IPD raw data from PTB.csv and data.csv
# and replaces NA, calculates times between events
# calculates exected values as proportions
# and then creates a list of costs for each scenario
# where amik and bdq have different durations
# this is saved in data_cleaned.RData


library(readr)
library(here)
library(dplyr)
library(purrr)

PTB_data <-
  read_csv(here::here("raw-data", "PTB.csv"),
           # read_csv(here::here("raw-data", "test_data.csv"),
           col_types = cols("Study number" = col_character(),
                            PTB = col_factor(levels = c(1, 0))))
raw_data <-
  read_csv(here::here("raw-data", "data.csv"),
           # read_csv(here::here("raw-data", "test_data.csv"),
           col_types = cols(firstivi = col_factor(levels = c("capreomycin", "streptomycin", "amikacin")),
                            hospital_discharge = col_datetime(format = "%d/%m/%Y"),
                            startIVs = col_datetime(format = "%d/%m/%Y"),
                            stopIVs = col_datetime(format = "%d/%m/%Y")))

# replace NA with 0
raw_data <-
  raw_data %>%
  mutate(
    num_Ak_trough_levels = ifelse(is.na(num_Ak_trough_levels), 0, num_Ak_trough_levels),
    num_line_complication = ifelse(is.na(num_line_complication), 0, num_line_complication),
    num_picc_line = ifelse(is.na(num_picc_line), 0, num_picc_line),
    num_hickman_line = ifelse(is.na(num_hickman_line), 0, num_hickman_line),
    readmissions_days = ifelse(is.na(readmissions_days), 0, readmissions_days),
    readmission_line_complication = ifelse(is.na(readmission_line_complication), 0, readmission_line_complication),
    IV_status = ifelse(is.na(IV_status), "none", IV_status),
    num_hearing_tests = ifelse(is.na(num_hearing_tests), 0, num_hearing_tests))

# time between events
raw_data <-
  raw_data %>%
  mutate(
    IV_days = difftime(stopIVs, startIVs, units = "days"),
    admission_days_inj = difftime(hospital_discharge, startIVs, units = "days"),
    OPAT_weeks_obs = pmax(0, as.numeric(difftime(stopIVs, hospital_discharge, units = "weeks"))))

raw_data <- na.omit(raw_data)

# raw_data <-
#   raw_data %>%
#   merge(PTB_data, by = "Study number") %>%
#   dplyr::filter(PTB == 0, # PTB only
#                 !"Study number" %in% c("RF10", "RF13")) # >0 days inj only

## expected rates

hickman_per_week <- with(raw_data,
                         sum(num_hickman_line)/sum(OPAT_weeks_obs))

picc_per_week <- with(raw_data,
                      sum(num_picc_line)/sum(OPAT_weeks_obs))

raw_data <- dplyr::rename(raw_data, centre_patientid = centre)

Tx_days <- list('obs' = raw_data$IV_days,
                '6mo_days' = rep(6*4*7, nrow(raw_data)),
                '8mo_days' = rep(8*4*7, nrow(raw_data)))

inj_scenarios <- names(Tx_days)

hear_freq_weeks <- c(1, rep(c(0,0,0,1), 100)) #monthly

data <- vector("list")
prop_less <- c('rho0' = 0,
               'rho0.1' = 0.1,
               'rho0.33' = 0.33)

MAX_WEEKS <- 100
blood_freq_weeks <- c(1,0,1, rep(c(0,0,0,1), MAX_WEEKS)) #monthly

k <- 1

for (j in inj_scenarios) {

  # amikacin
  inj_data <-
    raw_data %>%
    mutate(
      total_admission_days_inj = admission_days_inj + readmissions_days,
      OPAT_days_inj = Tx_days[[j]] - total_admission_days_inj,
      OPAT_days_inj = pmax(0, OPAT_days_inj),
      OPAT_weeks_inj = OPAT_days_inj/7,                                  # convert to weeks
      inj_scenario = j,
      num_hearing_tests = ifelse(inj_scenario == "obs",                  # fixed number of hearing tests
                                 yes = num_hearing_tests,
                                 no = 8), #guidelines (@Amber)
                                 # no = count_events(hear_freq_weeks,
                                 #                   OPAT_weeks_inj)),
      num_hickman_line = ifelse(inj_scenario == "obs",                   # expected number of hickman
                                yes = num_hickman_line,                  # using pooled rate
                                no = OPAT_weeks_inj*hickman_per_week),
      num_picc_line = ifelse(inj_scenario == "obs",                      # expected number of picc
                             yes = num_picc_line,                        # using pooled rate
                             no = OPAT_weeks_inj*picc_per_week),
      num_lines = num_picc_line + num_hickman_line,
      num_ECG_inj = 3 #guidelines (@Amber)                               # fixed number of ECG
    )

  for (i in prop_less) {

    data[[k]] <-
      inj_data %>%
      mutate(
        admission_excess_bdq = admission_days_inj*i,                       # how much time avoided in hospital
        admission_days_bdq = admission_days_inj*(1 - i),                   # scale down hospital duration

        total_admission_days_bdq = admission_days_bdq + readmissions_days,
        total_admission_weeks_bdq = total_admission_days_bdq/7,            # convert to weeks

        # OP_days_bdq = OPAT_days_inj + admission_excess_bdq,              # OPAT days same as amik inj
        OP_days_bdq = Tx_days$`6mo_days` - total_admission_days_bdq,       # OPAT for fixed total Tx of 6 months
        OP_days_bdq = pmax(0, OP_days_bdq),

        OPAT_weeks_obs = round(OPAT_weeks_obs, 2),
        num_hickman_line = round(num_hickman_line, 2),
        num_picc_line = round(num_picc_line, 2),
        num_lines = round(num_lines, 2),
        total_admission_weeks_bdq = round(total_admission_weeks_bdq, 2),

        num_blood_bdq = 8, #guidelines (@Amber)                            # fixed values
                           # count_events(freq_weeks = blood_freq_weeks,
                           #              num_weeks  = OP_days_bdq/7),
        num_ECG_bdq = 7) #guidelines (@Amber)

    k <- k + 1
  }
}

names(data) <-
  expand.grid(names(prop_less), inj_scenarios) %>%
  apply(1, paste, collapse = ".")


##########
## save ##
##########

save(data, file = here::here("raw-data", "data_cleaned.RData"))

# export to Excel

library(openxlsx)

OUT <- createWorkbook()
for (i in names(data))
{
  addWorksheet(wb = OUT, sheetName = i)
  writeData(OUT, sheet = i, x = data[[i]])
}
saveWorkbook(wb = OUT,
             file = here::here("raw-data", "data_cleaned.xlsx"),
             overwrite = TRUE)

