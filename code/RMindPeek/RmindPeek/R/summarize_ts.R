#' Summarize TS data
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
summarize_ts <- function(data) {
  data_sum <- data %>%
    group_by(region, state, stimulus, concentration, frame, FS, MS) %>%
    summarize(
      mean = mean(dff),
      median = median(dff),
      sd = sd(dff),
      n = n(),
      sem = sd / sqrt(n)
    ) %>%
    mutate(
      stim = paste(stimulus, concentration),
      id = paste(region, state, stimulus, concentration)
    )

  return(data_sum)
}
