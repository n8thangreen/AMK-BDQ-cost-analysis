
#' count_events
#'
#' @param freq_weeks
#' @param num_weeks
#'
#' @return
#' @export
#'
#' @examples
count_events <- function(freq_weeks,
                         num_weeks) {
  num_weeks %>%
    pmax(0) %>%
    map(seq_len) %>%
    map_dbl(~sum(freq_weeks[.]))
}
