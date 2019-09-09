---
title: "cost analyis"

author: "N Green"
date: "2018-10-16"
output:
  html_document:
    keep_md: TRUE
---


```r
library(readr)
library(here)
```

```
## here() starts at C:/Users/ngreen1/Google Drive/bedaquilie amikacin cost analysis/costAnalysis
```

```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
data <-
  read_csv(here::here("raw-data/data.csv"),
           col_types = cols(firstivi = col_factor(levels = c("capreomycin", "streptomycin", "amikacin")),
                            startIVs = col_datetime(format = "%d/%m/%Y"),
                            stopIVs = col_datetime(format = "%d/%m/%Y")))

data <-
  data %>%
  mutate(IV_days = stopIVs - startIVs,
         # adm_to_startIV = startIVs - admission,
         # OPAT_days = IV_days - adm_to_startIV
         num_Ak_trough_levels = ifelse(is.na(num_Ak_trough_levels), 0, num_Ak_trough_levels),
         num_line_complication  = ifelse(is.na(num_line_complication), 0, num_line_complication),
         num_picc_line  = ifelse(is.na(num_picc_line), 0, num_picc_line),
         num_hickman_line  = ifelse(is.na(num_hickman_line), 0, num_hickman_line),
         days_in_hospital  = ifelse(is.na(days_in_hospital), 0, days_in_hospital),
         num_hearing_tests  = ifelse(is.na(num_hearing_tests), 0, num_hearing_tests)
  )

str(data)
```

```
## Classes 'tbl_df', 'tbl' and 'data.frame':	100 obs. of  17 variables:
##  $ Study number                    : chr  "IK 1" "NP 2" "ST 3" "SA 4" ...
##  $ studynum                        : int  1 2 3 4 5 6 7 8 9 10 ...
##  $ centre                          : int  4 4 4 4 4 4 4 4 4 4 ...
##  $ admission                       : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ firstivi                        : Factor w/ 3 levels "capreomycin",..: 3 1 1 1 3 1 1 1 1 1 ...
##  $ startIVs                        : POSIXct, format: "2008-01-16" "2009-06-13" ...
##  $ stopIVs                         : POSIXct, format: "2008-05-02" "2009-10-14" ...
##  $ days_in_hospital                : num  73 4 197 104 240 26 78 81 35 127 ...
##  $ num_hearing_tests               : num  0 0 0 0 12 5 1 4 1 0 ...
##  $ num_Ak_trough_levels            : num  2 0 0 0 25 0 20 0 0 0 ...
##  $ num_picc_line                   : num  2 1 3 1 2 2 1 1 2 0 ...
##  $ num_hickman_line                : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ readmission_line_complication   : int  NA NA NA NA NA 1 NA NA NA NA ...
##  $ total_duration_readmissions_days: int  NA NA NA NA NA 21 NA NA NA NA ...
##  $ num_line_complication           : num  1 0 2 0 2 1 0 0 1 0 ...
##  $ IV_status                       : chr  NA NA NA NA ...
##  $ IV_days                         :Class 'difftime'  atomic [1:100] 107 123 197 94 168 187 18 183 179 87 ...
##   .. ..- attr(*, "units")= chr "days"
```

```r
save(data, file = here::here("data", "data.RData"))
```


---
title: "data-clean.R"
author: "ngreen1"
date: "Tue Oct 16 16:01:29 2018"
---
