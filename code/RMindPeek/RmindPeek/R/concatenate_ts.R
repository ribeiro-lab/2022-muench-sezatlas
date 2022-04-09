#' concatenate TS
#' concatenate time series data by changing the frame
#' @param data
#' @param type
#'
#' @return
#' @export
#'
#' @examples
concatenate_ts <- function(data, include_nostim = FALSE) {
  range <- range(data$frame)
  message("Range is: ", paste(range, collapse = ":"))
  len <- length(range[1]:range[2])

  if (include_nostim == FALSE){
    data <- data %>%
      mutate(stim = paste(stimulus, concentration))
    data$frame[which(data$stim == 'h2o 0')] <-
      data$frame[which(data$stim == 'h2o 0')] + len
    data$frame[which(data$stim == 'sucrose 200mM')] <-
      data$frame[which(data$stim == 'sucrose 200mM')] + len*2
    data$frame[which(data$stim == 'yeast 10%')] <-
      data$frame[which(data$stim == 'yeast 10%')] + len*3
  } else {
    data <- data %>%
      mutate(stim = paste(stimulus, concentration))
    data$frame[which(data$stim == 'nostim 0')] <-
      data$frame[which(data$stim == 'nostim 0')] + len
    data$frame[which(data$stim == 'h2o 0')] <-
      data$frame[which(data$stim == 'h2o 0')] + len*2
    data$frame[which(data$stim == 'sucrose 200mM')] <-
      data$frame[which(data$stim == 'sucrose 200mM')] + len*3
    data$frame[which(data$stim == 'yeast 10%')] <-
      data$frame[which(data$stim == 'yeast 10%')] + len*4
  }
    return(data)
}
