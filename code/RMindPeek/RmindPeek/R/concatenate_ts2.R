#' concatenate TS
#' concatenate time series data by changing the frame
#'
#' @param stimuli stimnuli in required order
#' @param data
#'
#' @return
#' @export
#'
#' @examples
concatenate_ts2 <- function(data, stimuli = c('nostim 0','h2o 0','sucrose 200mM','yeast 10%')) {
  data <- data %>%
    mutate(stim = paste(stimulus, concentration)) %>%
    filter(stim %in% stimuli)
  range <- range(data$frame)
  message("Range is: ", paste(range, collapse = ":"))
  len <- length(range[1]:range[2])

  for (i in 1:length(stimuli)){
    stimulus_i <- stimuli[i]
    data$frame[which(data$stim == stimulus_i)] <-
      data$frame[which(data$stim == stimulus_i)] + len*(i-1)
  }
  return(data)
}
