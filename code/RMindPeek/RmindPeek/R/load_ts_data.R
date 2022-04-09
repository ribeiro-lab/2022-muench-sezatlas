#' Load time series data
#'
#' @param path path containing ts data
#' @param ds names of the data sets to load (defaults to "all")
#' @param stim_exclude stimuli to exclude (character vector of format "stimulus concentration")
#' @param stim_include stimuli to exclude (character vector of format "stimulus concentration")
#'
#' @return data.frame with all ts data cleaned up
#' @export
#'
#' @examples
load_ts_data <- function(path = "", ds = NA, stim_exclude = NA, stim_include = NA) {
  file_source <- path
  files <- list.files(file_source, include.dirs = FALSE, pattern='.csv')
  data  <- data.frame()
  for (i in files) {
    data <- rbind(data,read.csv(paste0(file_source,i)))
  }

  # fix some flies
  fx <- which(is.na(data$stimulus))
  data$stimulus[fx] <- data$concentration[fx]
  data$concentration[fx] <- NA

  if(any(is.na(data$concentration))) {
    data[is.na(data$concentration),]$concentration <- 0
  }


  # clean up double stimulations within 1 animal
  tmp <- data %>% dplyr::filter(region == 1) %>%
    dplyr::select(state, fly, stimulus, concentration, stimulus_nr, recording) %>%
    dplyr::group_by(state, fly, stimulus, concentration) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    dplyr::filter(n > 1)

  if(nrow(tmp) > 0) {
    for (i in 1:nrow(tmp)) {
      stim_nrs <- unique(dplyr::filter(data, state == tmp$state[i], fly == tmp$fly[i], stimulus == tmp$stimulus[i], concentration == tmp$concentration[i])$stimulus_nr)
      stim_rm  <- tail(sort(stim_nrs),-1)
      indx_rm  <- which(data$state == tmp$state[i] &
                          data$fly == tmp$fly[i] &
                          data$stimulus == tmp$stimulus[i] &
                          data$concentration == tmp$concentration[i] &
                          data$stimulus_nr %in% stim_rm)
      data <- data[-indx_rm, ]
    }
  }


  data_t <- data %>% tidyr::gather(key = frame, value = value, -c(region:recording)) %>%
    dplyr::mutate(frame = as.numeric(substr(frame,2,99)))

  # correcting stimulus names & concentrations typos
  levels(data_t$stimulus) <- tolower(levels(data_t$stimulus))
  if (any(levels(data_t$stimulus) == 'surose')) levels(data_t$stimulus)[which(levels(data_t$stimulus) == 'surose')] <- 'sucrose'
  if (any(levels(data_t$stimulus) == 'sucroe')) levels(data_t$stimulus)[which(levels(data_t$stimulus) == 'sucroe')] <- 'sucrose'
  if (any(levels(data_t$stimulus) == 'h20')) levels(data_t$stimulus)[which(levels(data_t$stimulus) == 'h20')] <- 'h2o'

  data_t$FS <- NA
  data_t$MS <- NA
  data_t$FS[grep('ff', data_t$state)]     <- 'fed'
  data_t$FS[grep('10d', data_t$state)]    <- 'deprived'
  data_t$MS[grep('mated', data_t$state)]  <- 'mated'
  data_t$MS[grep('virgin', data_t$state)] <- 'virgin'

  data_t$region_index <- data_t$region

  # apply filters
  if(!any(is.na(stim_exclude)) & !any(is.na(stim_include))) {
    stop("Define EITHER stim_include OR stim_exclude, not both.")
  }

  n_row <- dim(data_t)[1]

  if(!any(is.na(ds))) {
    d_sets <- data_t$state
    data_t <- data_t[which(d_sets %in% ds),]
    message("Excluded: ", paste(unique(d_sets)[which(unique(d_sets) %ni% ds)], collapse = ", "))
    message("Returning: ", paste(unique(data_t$state), collapse = ", "))
  }


if(!any(is.na(stim_exclude))) {
    stim_conc <- paste(data_t$stimulus,data_t$concentration)
    data_t <- data_t[which(stim_conc %ni% stim_exclude),]
    message("Excluded: ", paste(unique(stim_conc)[which(unique(stim_conc) %in% stim_exclude)], collapse = ", "))
    message("Returning: ", paste(unique(paste(data_t$stimulus,data_t$concentration)), collapse = ", "))
  }

  if(!any(is.na(stim_include))) {
    stim_conc <- paste(data_t$stimulus,data_t$concentration)
    data_t <- data_t[which(stim_conc %in% stim_include),]
    message("Excluded: ", paste(unique(stim_conc)[which(unique(stim_conc) %ni% stim_include)], collapse = ", "))
    message("Returning: ", paste(unique(paste(data_t$stimulus,data_t$concentration)), collapse = ", "))
  }

  n_row_2 <- dim(data_t)[1]
  message("Returning ",n_row_2, "rows.\n", "Removed ", n_row - n_row_2, " rows.")
  return(droplevels(data_t))
}
