
#' ---
#' title: "cost analyis: paper_output_summary_table.R"
#'
#' author: "N Green"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     keep_md: TRUE
#' ---

# Paper cost-effectiveness output summary table
#
#           | Scenario      | Cost AK(sd; £)	| Cost BDQ(sd; £)	| Cost BDQ-cost AK(sd; £)	| P(BDQ cost-saving) | Max cost BDQ for P(cost-saving)≥0.95
# Observed  |               |
# AK vs BDQ	| 0%  reduction |
#           | 10% reduction |
#           | 33% reduction |
# 24 weeks  | 0%  reduction |
# AK vs BDQ | 10% reduction |
#           | 33% reduction |
#


library(dplyr)
library(knitr)
library(kableExtra)


## read-in density strip data from scripts.R
load("data/denstrip_dat.RData")


# AK_stats <- map_df(summary_table, function(x) x[1, ])
# BDQ_stats <- map_df(summary_table, function(x) x[2, ])

# delta_stats <- apply(denstrip_dat, 2, summary)
delta_stats <- psych::describe(denstrip_dat)

thresh <- apply(denstrip_dat, 2, function(x) sum(x < 0)/length(x))

## max(unit bdq cost : p(bdq cost-saving) >= 0.95)
thresh_unitbdq <- NULL
for (bdq in seq(4000, 18000, by = 100)) {

  thresh_unitbdq <- rbind(thresh_unitbdq,
                          c(bdq = bdq, apply(denstrip_dat, 2, function(x) sum(x - 18000 + bdq < 0)/length(x))))
}

cbdq_0.95 <-
  as.data.frame(thresh_unitbdq) %>%
  melt(id.vars = "bdq") %>%
  group_by(variable) %>%
  filter(value >= 0.95) %>%
  filter(bdq == max(bdq))


summary_long <-
  m.centre_long2 %>%
  group_by(j, variable) %>%
  summarize(mean = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            Uq = quantile(value, probs = 0.75, na.rm = TRUE),
            Lq = quantile(value, probs = 0.25, na.rm = TRUE)
  )

AK_stats <- summary_long %>% filter(variable == "amik")
BDQ_stats <- summary_long %>% filter(variable == "bdq")

out_tab <-
  data.frame(
    check.names = FALSE,
    "scenario" = names(summary_table),

    # "cost AK" = paste(round(AK_stats$median, 0), " (", round(AK_stats$sd, 0), ")", sep = ""),
    # "cost BDQ" = paste(round(BDQ_stats$median, 0), " (", round(BDQ_stats$sd, 0), ")", sep = ""),

    "cost AK (£)" = paste(round(AK_stats$mean, 0), " (", round(AK_stats$sd, 0), ")", sep = ""),
    "cost BDQ (£)" = paste(round(BDQ_stats$mean, 0), " (", round(BDQ_stats$sd, 0), ")", sep = ""),

    "c1 - c0 (£)" = paste(round(delta_stats$mean, 0), " (", round(delta_stats$sd, 0), ")", sep = ""),
    "P(BDQ cost saving)" = round(thresh, 2),
    "cost BDQ for P(cost-saving) >= 95% (£)" = cbdq_0.95$bdq
  )



print(
   out_tab %>%
    kable() %>%
    kable_styling("striped"))
