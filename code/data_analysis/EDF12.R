# regions with sign. effects (wilcox) --------------------------------

data_peak <-
  peak_response(
    data,
    method = "mean",
    frames = 1:5,
    bgf = -5:0,
    sd_mp = SD_cutoff
  )

data_peak <-
  filter(data_peak, stimulus %ni% c('h2o', 'nostim', 'sucrose'))


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
  factor(data_peak$stimulus,
         levels = c('nostim', 'h2o', 'sucrose', 'yeast'))
data_peak$region <- factor(data_peak$region)


plot_data <- data_peak

#order
order_median <-
  plot_data %>% filter(stimulus == 'yeast', MS == 'virgin', FS == 'deprived') %>%
  group_by(region) %>%
  summarize(median = median(peak)) %>%
  arrange(median) %>%
  mutate(order_median = row_number())

order_delta_S6 <-
  plot_data %>% filter(stimulus == 'yeast', MS == 'mated') %>%
  group_by(region, FS) %>%
  summarize(median = median(peak)) %>%
  tidyr::spread(FS, median) %>%
  mutate(delta = deprived - fed) %>%
  arrange(delta) %>%
  ungroup %>%
  mutate(order_delta = dplyr::row_number()) %>%
  select(region, order_delta)

order_delta <- plot_data %>%
  filter(stimulus == 'yeast', FS == 'deprived') %>%
  group_by(region, MS) %>%
  summarize(median = median(peak)) %>%
  tidyr::spread(MS, median) %>%
  mutate(delta = mated - virgin) %>%
  arrange(delta) %>%
  ungroup %>%
  mutate(order_delta = dplyr::row_number()) %>%
  select(region, delta, order_delta)

plot_data <- merge(plot_data, order_delta)

plot_data <- filter(plot_data, FS == 'deprived')

stat_wilcox <- plot_data %>%
  filter(FS == "deprived") %>%
  group_by(region, stimulus) %>%
  wilcox_test(peak ~ MS) %>%
  add_significance() %>%
  add_xy_position(x = "region", dodge = 0.85) %>%
  filter(p.signif != 'ns')


p <-
  ggplot(data = filter(plot_data, FS == 'deprived'), aes(x = reorder(region, order_delta), y = peak)) +
  geom_point(
    aes(color = MS),
    position = position_dodge(width = .85),
    alpha = .4,
    shape = 16,
    size = 2
  ) +
  geom_boxplot(
    aes(color = MS),
    position = position_dodge(width = .85),
    size = .6,
    outlier.shape = 21,
    outlier.alpha = 1,
    outlier.fill = 'white',
    outlier.colour = 'grey60',
    alpha = 0
  ) +
  scale_color_manual(values = ms_colors, 'reproductive state') +
  scale_fill_manual(values = ms_colors, 'reproductive state') +
  facet_grid(stimulus ~ .) +
  labs(x = 'atlas region', y = lab_df) +
  theme_dahaniel2(base_size = 28) +
  theme(
    axis.text.x  = element_text(
      angle = 0,
      hjust = .5,
      vjust = .5
    ),
    strip.text = element_text(size = 24),
    panel.grid.major.x = element_line(size = .3, color = 'grey80'),
    panel.grid.major.y = element_line(size = .2, color = 'grey80'),
    panel.border = element_blank(),
    legend.position = 'bottom'
  ) +
  coord_flip()

p <-
  p + stat_pvalue_manual(
    stat_wilcox,
    label = "p.signif",
    x = 'region',
    size = 9,
    y.position = -.3,
    hjust = 0
  )


ggsave(
  paste0(out_path, "all_boxplots_MS.pdf"),
  width = 80,
  height = 280,
  dpi = 300,
  scale = 3,
  units = 'mm',
  plot = p
)

write.csv(
  filter(plot_data, FS == 'deprived', stimulus == 'yeast') %>%
    select(region, MS, peak) %>%
    arrange(MS, region),
  paste0(out_path, "/source-data_FigS11.csv")
)







# Table S3 pooled fdr --------------------------------
pooled_data <- data_peak

stat_wilcox_y <- pooled_data %>%
  filter(FS == 'deprived', stimulus == 'yeast') %>%
  group_by(region, stimulus) %>%
  wilcox_test(peak ~ MS) %>%
  add_significance()

stat_wilcox_fdr_y <- pooled_data %>%
  filter(stimulus == 'yeast') %>%
  group_by(region, stimulus) %>%
  wilcox_test(peak ~ MS) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()


write.csv(stat_wilcox_y,
          paste0(out_path, "MS_stat_nocor_deprived_yeast.csv"))
write.csv(stat_wilcox_fdr_y,
          paste0(out_path, "MS_stat_fdr_pooled_yeast.csv"))


delta_mated_y <-
  plot_data %>% filter(stimulus == 'yeast', FS == 'deprived') %>%
  group_by(region, MS) %>%
  summarize(median = median(peak)) %>%
  tidyr::spread(MS, median) %>%
  mutate(delta = mated - virgin) %>%
  ungroup


delta_pooled_y <- plot_data %>% filter(stimulus == 'yeast') %>%
  group_by(region, MS) %>%
  summarize(median = median(peak)) %>%
  tidyr::spread(MS, median) %>%
  mutate(delta = mated - virgin) %>%
  ungroup


write.csv(delta_mated_y,
          paste0(out_path, "MS_delta_deprived_yeast.csv"))
write.csv(delta_pooled_y,
          paste0(out_path, "MS_delta_pooled_yeast.csv"))
