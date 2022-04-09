# load data ---------------------------------------------------------------
ts_path <- "data/imaging/time_series_1Msucrose/"

load_atlas_labels_v2()

data <- load_ts_data(
  path = ts_path,
) %>%
  normalize_ts(method = "dff") %>%
  align_onsets(align_region = as.numeric(names(atlas_labels[which(atlas_labels == "PMS4 left")])), sd_multi = 4) %>%
  mutate(id = paste(region, state, stimulus, concentration)) #%>%

data_peak <-
  peak_response(
    data,
    method = "mean",
    frames = 1:5,
    bgf = -5:0,
    sd_mp = SD_cutoff
  )



# 1M ff vs dep ------------------------------------------------------------

data_peak <-
  peak_response(
    filter(data, stimulus %in% c('sucrose'), concentration %in% c("1M")),
    method = "mean",
    frames = 1:5,
    bgf = -5:0,
    sd_mp = SD_cutoff
  )
r = c(58, 59, 21, 9)
data_peak <- filter(data_peak, region %in% r)
data_peak$FS <- factor(data_peak$FS, levels = c('fed', 'deprived'))
data_peak$region <- factor(data_peak$region, levels = c(58, 59, 21, 9))

stat_wtest <- data_peak %>%
  group_by(region, stimulus, concentration) %>%
  wilcox_test(peak ~ FS) %>%
  add_significance()

plot_data <- data_peak

p <- ggplot(data = plot_data, aes(x = FS, y = peak)) +
  ggbeeswarm::geom_beeswarm(
    aes(color = FS),
    size = 3.5,
    shape = 16,
    cex = 3,
    alpha = .6
  ) +
  geom_boxplot(
    aes(color = FS),
    width = .5,
    outlier.alpha = 0,
    size = 1,
    outlier.shape = 16,
    outlier.size = 2,
    outlier.colour = "black",
    fill = NA
  ) +
  scale_color_manual(values = fs_colors,
                     'metabolic\nstate',
                     labels = as_labeller(lab_state)) +
  scale_fill_manual(values = fs_colors, 'metabolic\nstate') +
  facet_grid(. ~ region, scales = "free") +
  labs(y = lab_df) +
  theme_dahaniel2(base_size = 24) +
  theme(
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    strip.text = element_text(size = 20),
    legend.position = 'right'
  ) +
  coord_cartesian(clip = "off")
p <-
  p + stat_pvalue_manual(
    stat_wtest,
    label = "{formatC(p,2)}",
    size = 5,
    bracket.size = .6,
    tip.length = 0,
    y.position = 1.1
  )
p
ggsave(
  paste0(out_path, "/EDF9.pdf"),
  width = 80,
  height = 45,
  dpi = 300,
  scale = 3,
  units = 'mm',
  plot = p
)
write.csv(plot_data %>% select(region, FS, peak),
          paste0(out_path, "source-data_EDF9.csv"))
