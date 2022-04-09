

# load data ---------------------------------------------------------------
ts_path <- "data/imaging/time_series_2%yeast/"

load_atlas_labels_v2()

data <- load_ts_data(
  path = ts_path
) %>%
  normalize_ts(method = "dff") %>%
  align_onsets(align_region = as.numeric(names(atlas_labels[which(atlas_labels == "PMS4 left")])), sd_multi = 4) %>%
  mutate(id = paste(region, state, stimulus, concentration))

data_peak <-
  peak_response(
    data,
    method = "mean",
    frames = 1:5,
    bgf = -5:0,
    sd_mp = SD_cutoff
  )


# boxplots yeast 2% -----------------------------------------------
data_peak <-
  peak_response(
    data,
    method = "mean",
    frames = 1:5,
    bgf = -5:0,
    sd_mp = SD_cutoff
  )

stat_wilcox <- data_peak %>%
  filter(MS == "mated") %>%
  group_by(region, stimulus, concentration) %>%
  wilcox_test(peak ~ FS) %>%
  add_significance()

regions <-  c(25, 58, 21)

data_peak$MS <- factor(data_peak$MS, levels = c("virgin", "mated"))
data_peak$FS <- factor(data_peak$FS, levels = c("fed", "deprived"))
data_peak$MSFS <- paste0(data_peak$MS, ":", data_peak$FS)
data_peak$MSFS <-
  factor(
    data_peak$MSFS,
    levels = c(
      'virgin:fed',
      'virgin:deprived',
      'mated:fed',
      'mated:deprived'
    )
  )
data_peak$stimulus <-
  factor(
    data_peak$stimulus,
    levels = c('nostim', 'h2o', 'sucrose', 'sucrose1', 'yeast', 'yeast1')
  )

plot_stats <-
  dplyr::filter(stat_wilcox, region %in% regions, stimulus == 'yeast')
if (nrow(plot_stats) != 0)
  plot_stats$y.position <- 1
plot_data  <-
  filter(data_peak, region %in% regions, stimulus %in% 'yeast')

range <- range(plot_data$peak, na.rm = T)
range <- range[2] - range[1]
step  <- range / 14

if (nrow(plot_stats) != 0) {
  for (r in unique(plot_stats$region)) {
    max_y <- max(filter(plot_data, region == r)$peak, na.rm = T)
    pos <- which(plot_stats$region == r)
    plot_stats$y.position[pos] <- seq(length(pos)) * step + max_y
  }
}

plot_stats$y.position <- c(.18, .34, .89)

atlas_labels_plot <- c('25' = 'AMS1\n25',
                       '58' = 'PMS4\n58',
                       '21' = 'outer borboleta\n21')

p <- ggplot(data = plot_data, aes(x = FS, y = peak)) +
  geom_boxplot(
    aes(color = FS),
    width = .5,
    outlier.alpha = 0,
    size = 1,
    outlier.shape = 16,
    outlier.size = 2,
    outlier.colour = "black"
  ) +
  ggbeeswarm::geom_beeswarm(
    aes(color = FS),
    size = 3.5,
    shape = 16,
    cex = 3,
    alpha = .6
  ) +
  scale_color_manual(values = fs_colors,
                     'metabolic\nstate',
                     labels = as_labeller(lab_state)) +
  scale_fill_manual(values = fs_colors, 'metabolic\nstate') +
  facet_wrap(. ~ region,
             scales = "free_y",
             labeller = as_labeller(atlas_labels_plot)) +
  labs(y = lab_df) +
  theme_dahaniel2(base_size = 24) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.border = element_blank(),
    strip.text = element_text(size = 20),
    legend.position = 'right'
  ) +
  expand_limits(y = 0) +
  coord_cartesian(clip = "off")

if (nrow(plot_stats) != 0) {
  p <-
    p + stat_pvalue_manual(
      plot_stats,
      label = "{formatC(p,2)}",
      size = 5,
      bracket.size = .6,
      tip.length = 0.00
    )
}
plot(p)
ggsave(
  paste0(out_path, "EDF10.pdf"),
  width = 80,
  height = 50,
  dpi = 300,
  scale = 3,
  units = 'mm',
  plot = last_plot()
)
write.csv(plot_data %>% select(region, FS, peak),
          paste0(out_path, "source-data_EDF10.csv"))
