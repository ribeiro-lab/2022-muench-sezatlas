#' Plot data onto atlas representation
#'
#' @param data
#' @param atlas
#' @param value_column
#' @param region_column
#' @param value_name
#' @param plot_legend
#' @param fill_color
#' @param line_color line color, set to "MAP" in order to map to the value column
#' @param line_width
#' @param label
#' @param alpha alpha value of the polygons
#' @param label_color
#' @param label_size
#'
#' @return
#' @export
#'
#' @examples
plot_atlas <- function(data = NULL,
                       atlas = sez_atlas,
                       value_column = "value",
                       region_column = "region",
                       value_name = "value",
                       plot_legend = FALSE,
                       fill_color = "grey70",
                       line_color = "grey30",
                       alpha = 1,
                       line_width = .2,
                       # plot_limits = NULL,
                       label = FALSE,
                       label_color = "grey20",
                       label_size = 2
                       ){

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop(
      "ggplot2 is required for atlas plotting, please install via
       install.packages('ggplot2')",
      call. = FALSE
    )

  if (!requireNamespace("dplyr", quietly = TRUE))
    stop(
      "dplyr is required for atlas plotting, please install via
       install.packages('dplyr')",
      call. = FALSE
    )

  if (!is.null(data)){
    atlas$value <- data[[value_column]][match(atlas$region, data$region)]
  } else {
    atlas$value <- atlas$region
  }

  # if (!is.null(data) & is.null(plot_limits)) {
  #   plot_limits = range(atlas$value, na.rm = T)
  # }

  plot <- ggplot2::ggplot(data = atlas)
  if (!is.null(data)) {
    if (line_color == "MAP") {
      plot <- plot + ggplot2::geom_polygon(ggplot2::aes(x = x, y = y, fill=value, color=value, group=group), size = line_width, alpha = alpha)
    } else {
      plot <- plot + ggplot2::geom_polygon(ggplot2::aes(x = x, y = y, fill=value, group=group), color=line_color, size = line_width, alpha = alpha)
      # plot <- plot + ggplot2::scale_fill_viridis_c(option = "C")
    }

  } else {
    plot <- plot + ggplot2::geom_polygon(ggplot2::aes(x = x, y = y, group=group), fill = fill_color, color=line_color, size = line_width, alpha = alpha)
  }
  plot <- plot + ggplot2::coord_fixed() +
    theme_dahaniel2(base_size = 5) +
    ggplot2::theme(axis.text =  ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          axis.text.x  = ggplot2::element_blank(),
          axis.text.y  = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank())

  if (plot_legend == FALSE) {plot <- plot + ggplot2::theme(legend.position = "none")}


  if (label == TRUE){
    label_centers <- atlas %>%
      group_by(region, id) %>%
      summarize(x = max(x)-((max(x)-min(x))/2), y = max(y)-((max(y)-min(y))/2))
    if (label_color == "MAP") {
      plot <- plot + geom_text(data = label_centers, aes(x = x, y = y, label = region, color = region), size = label_size)
    } else {
      plot <- plot + geom_text(data = label_centers, aes(x = x, y = y, label = region), size = label_size, color = label_color)
    }

  }


  return(plot)

}
