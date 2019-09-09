
#
replace_NA_cols <- function(dat, fn)
  apply(dat, 2,
        function(x) ifelse(is.na(x), fn(x, na.rm = TRUE), x))
