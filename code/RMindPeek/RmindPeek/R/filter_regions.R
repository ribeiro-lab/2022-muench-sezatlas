
#' filter_region
#' Filter regions based on: "response probablility"
#'
#' @param data DFF time series dataframe
#' @param method filtering method. P == (response probability) or R for responses above the SD rersponse threshold
#' @param SD multiplicator. How many times pre-stimulus SD to count as a response?
#' @param P probability threshold for filter.
#' @param bg Frames for calculating background SD.
#' @param res Frames of response time window.
#'
#' @return vector of regions that pass filter in at least one internal state
#' @export
#'
#' @examples
filter_regions <-
  function(data,
           method = "P",
           SD = 4,
           P = 0.5,
           bg = -5:0,
           res = 1:5) {
    data_peak <-
      peak_response(data, method = "mean", frames = res, bgf = bg, sd_mp = SD)

    if (method == "R") {
      data_peak <- filter(data_peak, response == TRUE)

      all_regions <- unique(data$region)
      regions <- unique(data_peak$region)

      message(
        "filter_ts: Returning ",
        length(regions),
        " of ",
        length(all_regions),
        " regions"
      )
    }


    if (method == "P") {
      response_p <-
        data_peak %>% group_by(region, state, stimulus, concentration, FS, MS) %>%
        summarize(
          n = n(),
          sum = sum(response),
          p = sum / n,
          response = mean(peak)
        ) %>%
        filter(n > 5) %>%
        mutate(
          stim = paste(stimulus, concentration),
          is = paste(MS, FS),
          id = paste(region, state, stimulus, concentration)
        ) %>%
        group_by(region) %>%
        mutate(response_norm = response / max(response, na.rm = T)) %>%
        ungroup() %>%
        filter(p >= P)

      all_regions <- unique(data$region)
      regions <- unique(response_p$region)

      message(
        "filter_ts: Returning ",
        length(regions),
        " of ",
        length(all_regions),
        " regions"
      )

    }

data_filtered <- subset(data, region %in% regions)

message(
  "\n\nfilter_ts: Returning ",
  nrow(data_filtered),
  " of ",
  nrow(data),
  " rows"
)

return(data_filtered)

  }
