#' %ni%
#' Inverse of %in%
#'
#' @param x
#' @param table
#'
#' @return
#' @export
#'
#' @examples
`%ni%` <- function (x, table)
{
  match(x, table, nomatch = 0) == 0
}

#' Dahaniel custom ggplot theme
#'
#' @param base_size
#' @param base_family
#'
#' @return
#' @export
#'
#' @examples
theme_dahaniel2 <- function(base_size = 5, base_family = "Arial", base_line_size = base_size/22, base_rect_size = base_size/22) {
  require(ggplot2)
  #require(ggthemes)
  theme_bw(base_size = base_size, base_family = base_family, base_line_size = base_line_size, base_rect_size =base_rect_size) +
    theme(
      axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey80"),
      panel.grid = element_blank(),
      axis.line = element_blank(), #element_line(colour = "black", size = .2),
      strip.background = element_blank(),
      strip.text = element_text(face = 'bold', size = 6),
      legend.key.size = unit(.3,"line"))
}
