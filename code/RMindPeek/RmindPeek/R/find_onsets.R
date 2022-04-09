#' Find stimulus onsets
#'
#' @param data dataframe (of dff data, for now (data_t))
#' @param frame_filter max frame to look for onset
#' @param sd_multi value to multiply bg sd by to define onset detection threshold
#' @param bg background frames for noise SD level calculation
#'
#' @return
#' @export
#'
#' @examples
find_onsets <-
  function(data,
           frame_filter = 30,
           bg = -5:-1,
           sd_multi = 4
           ) {
    # calculate stimulus onsets -----------------------------------------------
    data_onsets <-
      data %>% filter(frame %in% min(frame):frame_filter) %>%
      filter(
        stimulus %ni% c("nostim")
      ) %>%
      group_by(fly, recording, region, stimulus, concentration) %>%
      arrange(frame) %>%
      mutate(
        mean_bg = mean(dff[frame %in% bg]),
        sd_bg = sd(dff[frame %in% bg]),
        on = (dff > (mean_bg + (sd_bg * sd_multi))) &
          frame > max(bg),
        onset = min(frame[which(on == TRUE)]),
        onset = ifelse(is.infinite(onset), NA, onset)
      ) %>%
      ungroup() %>%
      filter(frame == 0) %>%
      select(region, state, fly, stimulus, concentration, stimulus_nr, recording, FS, MS, onset)

    return(data_onsets)
  }
