#' Load atlas labels
#'
#' @return atlas_labels variable
#' @export
#'
#' @examples
load_atlas_labels <- function() {
  assign(
    "atlas_labels",
    c('110' = 'peg projections',
      '111' = 'sensilla projections',
      '120' = 'borboleta',
      '121' = 'inner borboleta',
      '122' = 'outer borboleta',
      '1'   = 'central 1',
      #'42'  = 'central 2',
      '87'  = 'central 3',
      '93'  = 'pegs left',
      '94'  = 'pegs right',
      '66'  = 'sensilla left',
      '67'  = 'sensilla right',
      '20'  = 'bitter',
      '36'  = 'borbo outer left',
      '83'  = 'borbo outer right',
      '8'   = 'borbo inner left',
      '46'  = 'borbo inner right',
      '42'  = 'MN 1',
      '82'  = 'MN 2',
      '97'  = 'MN 3'
    ),
    envir = .GlobalEnv)
}

