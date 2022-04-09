data_peak_sum <- data_peak %>%
  filter(MS == 'mated') %>%
  group_by(region, MS, FS, stimulus) %>%
  summarize(
    mean = mean(peak),
    n = n(),
    sd = sd(peak),
    sem = sd / sqrt(n),
    median = median(peak),
    max = max(peak)
  ) %>%
  ungroup()

for (fs in c('fed', 'deprived')) {
  mx = max(filter(data_peak_sum)[, measure], na.rm = TRUE)
  for (s in c('nostim', 'h2o', 'sucrose', 'yeast')) {
    plot_atlas(
      filter(data_peak_sum, stimulus == s, FS == fs),
      value_column = measure,
      plot_legend = TRUE,
      line_color = "grey20",
      line_width = .1
    ) +
      scale_fill_gradientn(
        colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")),
        limits = c(-mx, mx),
        breaks = c(round(-mx, 2), 0, round(mx, 2)),
        lab_df
      ) +
      theme(legend.position = 'right')
    ggsave(
      paste0(out_path, '/EDF4_', s, '_', fs, '.pdf'),
      width = 120,
      height = 25,
      dpi = 300,
      scale = 1,
      plot = last_plot(),
      unit = 'mm'
    )
  }
}