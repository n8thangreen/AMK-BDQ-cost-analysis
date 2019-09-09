
#' ---
#' title: "cost analyis: twoway-analysis.R"
#'
#' author: "N Green"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     keep_md: TRUE
#' ---


library(dplyr)

load(here::here("data", "costs.RData"))
load(here::here("data", "cost_bdq.RData"))

load(file = here::here("raw-data", "data_cleaned.RData"))
t_adm <- as.numeric(median(data$rho0.obs$total_admission_days_inj))

c_under20 <- costs$mean[1, "bed_under20d"]
c_over20 <- costs$mean[1, "bed_over20d"]
days_Tx <- 6*7*4
c_ECG <- costs$mean[1, "ECG"]
c_renal <- costs$mean[1, "blood_renal"]
c_liver <- costs$mean[1, "blood_liver"]
c_hick <- costs$mean[1, "hick"]
c_picc <- costs$mean[1, "PICC"]
c_OPAT <- 102
c_hear <- costs$mean[1, "hear"]
c_bdq <- 18000

n_hear <- 1
n_picc <- 1
n_hick <- 0

inputs <- expand.grid(rho = seq(0, 1, 0.01),
                      c_bdq = seq(9000, c_bdq, 100))

for (i in seq_len(nrow(inputs))) {

  c_total_bdq <- c_under20 + max(inputs$rho[i]*t_adm - 19, 0)*c_over20 + max(days_Tx - inputs$rho[i]*t_adm, 0)/7*(c_ECG + c_renal + c_liver) + inputs$c_bdq[i]

  c_total_amik <- c_under20 + max(t_adm - 19, 0)*c_over20 + c_hick*n_hick + c_picc*n_picc + c_OPAT*max(days_Tx - t_adm, 0) + c_hear*n_hear

  inputs$c_total[i] <- c_total_amik - c_total_bdq
}

library(ggplot2)
library(zoo)

ggplot(inputs, aes(rho, c_bdq)) + geom_tile(aes(fill = c_total)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_bw()

# https://www.r-bloggers.com/controlling-heatmap-colors-with-ggplot2/
## use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
color_palette <- colorRampPalette(c("#d53e4f", "#fc8d59", "#fee08b", "#ffffbf", "#e6f598", "#99d594", "#3288bd"))(20)
inputs$c_interv <- findInterval(inputs$c_total, seq(-30000, 20000, by = 1000))#, all.inside = TRUE)
label_text <- cut(inputs$c_total, seq(-30000, 20000, by = 1000)) %>% levels()

ggplot(inputs, aes(x = rho, y = c_bdq, fill = factor(c_interv))) +
  # geom_tile(colour = "black") +
  geom_tile() +
  # scale_fill_gradient(low = "white", high = "steelblue") +
  scale_fill_manual(values = color_palette, name = "", labels = label_text) +
  theme_bw()
