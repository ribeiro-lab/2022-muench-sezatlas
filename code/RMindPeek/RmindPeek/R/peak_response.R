#' Calculate peak responses
#'
#' @param data dff data (TODO: generalize to other data)
#' @param method calculation method, defaults to "mean" of ´frames´
#' @param frames frames to perform the calculation on
#' @param bgf background frames
#' @param sd_mp multiplyer for SD threshold response detection (peak > sd_bg * sd_mp)
#' use same parameters as in `filter_regions`.
#'
#' @return data frame containing peak values
#' @export
#'
#' @examples
peak_response <-
  function(data,
           method = "mean",
           frames = 1:5,
           bgf = -5:0,
           sd_mp = 2) {
    if (method == "mean") {
      data_sum <- data %>%
        dplyr::group_by(region,
                 state,
                 stimulus,
                 concentration,
                 fly,
                 recording,
                 FS,
                 MS) %>%
        dplyr::summarize(
          peak = mean(dff[frame %in% frames]),
          sd_peak = sd(dff[frame %in% frames]),
          bg = mean(dff[frame %in% bgf]),
          sd_bg = sd(dff[frame %in% bgf])
        ) %>%
        dplyr::mutate(response = ifelse(abs(peak) > (abs(bg) + (sd_bg * sd_mp)), TRUE, FALSE)) %>%
        dplyr::ungroup()

      return(data_sum)
    }

  }
