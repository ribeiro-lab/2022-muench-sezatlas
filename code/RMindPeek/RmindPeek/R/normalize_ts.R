#' normalize time-series data
#'
#' Normalizes TS data e.g. by calculating DF/F.
#'
#' @param data input data (format as in data_t)
#' @param bg_frames background frames for f0 calculation.
#' @param method for now only 'dff'
#' @param value_column
#' @param region region id for "one_region" option
#'
#' @return
#' @export
#'
#' @examples
normalize_ts <- function(data, bg_frames = 4:9, value_column = 'value', method = 'dff', regionx = NULL){
if (method == 'dff') {
  data <- data %>%
    dplyr::group_by(recording, region) %>%
    dplyr::mutate(f0 = mean(get(value_column)[frame %in% bg_frames]), dff = (get(value_column) - f0) / f0) %>%
    dplyr::ungroup()

  na_recs <- unique(data$recording[which(is.na(data$dff))])
  if (length(na_recs) > 0) {
    na_animals <- unique(data$fly[which(is.na(data$dff))])
    na_regions <- unique(data$region[which(is.na(data$dff))])

    data <- data[-which(data$recording %in% na_recs & data$region %in% na_regions),]
    warning(paste0('Created NAs in ', length(na_recs), ' recordings from ', length(na_animals), ' animals. (', paste(na_animals, collapse = ', '),  '), recordings were excluded. Regions involved: ', paste(na_regions, collapse = ', ')))

  }

}

if (method == 'z-score') {
  data <- data %>%
    dplyr::group_by(recording, region) %>%
    dplyr::mutate(z = (get(value_column) - mean(get(value_column))) / sd(get(value_column))) %>%
    dplyr::ungroup()
  #
  # na_recs <- unique(data$recording[which(is.na(data$dff))])
  # na_animals <- unique(data$fly[which(is.na(data$dff))])
  # na_regions <- unique(data$region[which(is.na(data$dff))])
  #
  # data <- data[-which(data$recording %in% na_recs & data$region %in% na_regions),]
}


if (method == 'shift-baseline') {
  data <- data %>%
    dplyr::group_by(recording, region) %>%
    dplyr::mutate(f0 = mean(get(value_column)[frame %in% bg_frames]), value_shifted = get(value_column) - f0) %>%
    dplyr::ungroup()
}

if (method == 'norm01') {
  data <- data %>%
    dplyr::group_by(recording, region) %>%
    dplyr::mutate(f0 = mean(get(value_column)[frame %in% bg_frames]), normalized = get(value_column) - f0, normalized = normalized / max(normalized)) %>%
    dplyr::ungroup()
}

  if (method == 'one_region') {
    data <- data %>%
      dplyr::group_by(recording) %>%
      dplyr::mutate(f0 = mean(get(value_column)[frame %in% bg_frames & region == regionx]), normalized = (get(value_column) - f0) / f0) %>%
      dplyr::ungroup()
  }


  return(data)
}
