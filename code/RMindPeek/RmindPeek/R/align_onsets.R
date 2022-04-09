#' Align to stimulus onset
#' TODO: base alignment on mean of region 66 & 67
#'
#' @param data dataframe (of dff data, for now (data_t))
#' @param align_region region to align to (usually sensilla input region)
#' @param frame_filter max frame to look for onset
#' @param sd_multi value to multiply bg sd by to define onset detection threshold
#' @param onset_intended intended stimulus onset (to shift nostim and when no response is detected)
#' @param cutframes frames to reurn after shifting
#' @param bg background frames for noise SD level calculation
#' @param value_column name of the column that contains the values
#'
#' @return
#' @export
#'
#' @examples
align_onsets <-
  function(data,
           align_region = 66,
           frame_filter = 15,
           bg = 1:8,
           sd_multi = 4,
           onset_intended = 10,
           cutframes = -5:45,
           value_column = 'dff') {
    # calculate stimulus onsets -----------------------------------------------
    data_onsets <-
      data %>% filter(region == align_region, frame %in% 0:frame_filter) %>%
      filter(
        stimulus %ni% c("nostim")
      ) %>%
      group_by(fly, stimulus, concentration, recording) %>%
      arrange(frame) %>%
      mutate(
        mean_bg = mean(get(value_column)[frame %in% bg]),
        sd_bg = sd(get(value_column)[frame %in% bg]),
        on = (get(value_column) > (mean_bg + (sd_bg * sd_multi))) &
          frame > 6,
        onset = min(frame[which(on == TRUE)]),
        onset_raw = onset,
        onset = ifelse(is.infinite(onset), onset_intended, onset)
      ) %>%
      ungroup() %>%
      filter(frame == 0) %>%
      select(region, state, fly, stimulus, concentration, recording, onset, onset_raw)

    # correct onsets ----------------------------------------------------------
    if (any(duplicated(data_onsets$recording) == TRUE)) {
      message("Duplicates detected")
    }

    # correct onsets ----------------------------------------------------------
    if (any(is.infinite(data_onsets$onset_raw))) {
      warning("Onsets could not be detected in all cases, check 'detected_onset column' for 'Inf'. Respective responses were shifted by 'intended_onset' value.")
    }

    data_corrected <- data %>%
      mutate(
        group = paste(recording, stimulus, concentration, region, MS, FS),
        frame_old = frame,
        detected_onset = NaN
      )

    for (i in 1:dim(data_onsets)[1]) {
      indx <- which(data_corrected$recording == data_onsets$recording[i])
      data_corrected[indx, ]$frame <-
        data[indx, ]$frame - (data_onsets$onset[i] - 1) # remove one to have frame 1 as the first response frame, not frame 0

      indx2 <- which(data_corrected$recording == data_onsets$recording[i] & data_corrected$region == align_region)
      data_corrected[indx2, ]$detected_onset <-
        data_onsets$onset_raw[i]
    }


    # shift "nostim" to 0 -----------------------------------------------------

    data_corrected[data_corrected$stimulus == "nostim", "frame"] <-
      data_corrected[data_corrected$stimulus == "nostim", "frame"] - onset_intended

    data_corrected <- data_corrected %>% filter(frame %in% cutframes)

    return(data_corrected)
  }
