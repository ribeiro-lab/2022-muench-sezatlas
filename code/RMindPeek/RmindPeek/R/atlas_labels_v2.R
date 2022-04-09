#' Load atlas labels
#'
#' loads atlas labels for atlas v2
#'
#' @return atlas_labels variable
#' @export
#'
#' @examples
load_atlas_labels_v2 <- function() {
  assign(
    "atlas_labels",
    c('25'  = 'AMS1 left',
      '26'  = 'AMS1 right',
      '58'  = 'PMS4 left',
      '59'  = 'PMS4 right',
      '68'  = 'bitter',
      '21'  = 'borboleta outer left',
       '9'  = 'borboleta outer right',
      '60'  = 'borboleta inner left',
      '57'  = 'borboleta inner right',
      '62'  = 'borboleta outer ventral left',
      '63'  = 'borboleta outer ventral right',
      '19'  = 'MN 1',
      '73'  = 'MN 2',
      '74'  = 'MN 3',
      '71'  = 'MN 4'
    ),
    envir = .GlobalEnv)
}

