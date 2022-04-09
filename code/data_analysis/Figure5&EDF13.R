

# MS changes top 10 from Extended Data Figure 12 ----------------------------------------
# order_delta are calculated in EDF12.R
order_delta

plot_data <- tail(order_delta, 8)
plot_data <-
  rbind(plot_data, order_delta[which(order_delta$region %in% c(3, 8)), ])
rng <- range(plot_data$delta)
plot_atlas(
  filter(plot_data),
  value_column = "delta",
  plot_legend = TRUE,
  line_color = "grey40",
  line_width = .1
) +
  scale_fill_gradientn(
    colours = c('white', ms_colors['mated']),
    na.value = 'grey100',
    limits = c(0, rng[2]),
    breaks = c(0, .1, .2),
    expression("" * Delta * "response")
  ) +
  theme(legend.position = 'right')

ggsave(
  paste0(out_path, 'delta_response_MS.pdf'),
  width = 160,
  height = 40,
  scale = plot_scale,
  units = plot_units,
  plot = last_plot()
)

# all peak maps -----------------------------------------------------------
# select all regions in specific slices to scale colors accordingly
slice9 <- c(1,4,7,8,11,13,14,15,16,17,19,23,24,27,30,31,32,34,36,39,40,41,42,44,46,50,52,54,66,67,69,70,71,72,73,74,75,76,77,79)
slice14 <- c(1,2,5,6,7,8,9,10,11,16,23,25,26,27,28,29,38,43,45,46,47,51,52,54,56,57,58,59,60,61,62,63,64,65,68,72,76,77,78,79,80)


data_peak_sum <- data_peak %>%
  filter(stimulus == 'yeast') %>%
  group_by(region, MS, FS, stimulus) %>%
  dplyr::summarize(
    mean = mean(peak),
    n = n(),
    sd = sd(peak),
    sem = sd / sqrt(n),
    median = median(peak),
    max = max(peak),
    .groups = 'drop'
  ) %>%
  ungroup()

for (s in c('slice9', 'slice14')) {
  plot_data <- data_peak_sum %>%
    filter(region %in% get(s))
  
  range = range(plot_data$mean, na.rm = TRUE)
  
  for (ms in c('virgin', 'mated')) {
    for (fs in c('fed', 'deprived')) {
      plot_atlas(
        filter(data_peak_sum, stimulus == 'yeast', MS == ms, FS == fs),
        value_column = 'mean',
        plot_legend = TRUE,
        line_color = "grey40",
        line_width = .1
      ) +
        scale_fill_gradientn(
          colors = rev(RColorBrewer::brewer.pal(11, "RdBu")),
          na.value = 'white',
          "response\nmagnitude",
          limits = c(-range[2], range[2]),
          breaks = c(-floor(range[2] * 10) / 10, 0, floor(range[2] * 10) / 10)
        ) +
        labs(title = paste(s, ms, fs, sep = ', ')) +
        theme(legend.position = 'bottom')
      ggsave(
        paste0(
          out_path,
          'map_dff_RdBu_slicescale_mean_',
          s,
          '_',
          ms,
          '_',
          fs,
          '.pdf'
        ),
        width = 160,
        height = 30,
        dpi = plot_dpi,
        scale = plot_scale,
        units = plot_units,
        plot = last_plot()
      )
    }
  }
}


