# scatterplots ------------------------------------------------------------

plot_data <-
  data_peak %>% filter(MS == 'mated', stimulus %in% c('sucrose', 'yeast')) %>%
  select(region, FS, stimulus, fly, peak) %>%
  group_by(region, FS, stimulus) %>%
  summarize(peak = mean(peak)) %>%
  spread(FS, peak)

plot_data$highlight <- 'other'
plot_data$highlight[which(plot_data$region %in% c(9, 21, 57, 60, 62, 63))] <-
  'borboleta'

ggplot(plot_data, aes(y = deprived, x = fed)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    alpha = .2,
    linetype = 2,
    size = .3
  ) +
  geom_smooth(
    method = "lm",
    color = 'grey60',
    formula = y ~  x,
    size = .3
  ) +
  geom_point(
    alpha = .4,
    aes(color = highlight),
    size = .5,
    shape = 16
  ) +
  scale_color_manual(values = c('red', 'black')) +
  stat_regline_equation(
    label.y = .1,
    label.x = .16,
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~")),
    color = "grey50",
    size = 1.5
  ) +
  facet_wrap( ~ stimulus, nrow = 1, labeller = as_labeller(lab_stim)) +
  coord_fixed(xlim = c(-.05, .85),
              ylim = c(-.05, .85),
              clip = "off") +
  labs(y = 'response magnitudes\nin deprived flies', x = 'response magnitudes in fully fed flies') +
  scale_y_continuous(breaks = c(0, 0.4, 0.8)) +
  scale_x_continuous(breaks = c(0, 0.4, 0.8)) +
  theme_dahaniel2(base_size = 5) +
  theme(
    panel.border = element_blank(),
    legend.position = 'none',
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
ggsave(
  paste0(out_path, "scatterplot_FS.pdf"),
  plot = last_plot(),
  width = 44,
  height = 28,
  scale = plot_scale,
  units = plot_units
)

write.csv(
  filter(plot_data, stimulus == 'yeast') %>%
    select(region, deprived, fed) %>%
    arrange(region),
  paste0(out_path, "/source-data_Fig4a_yeast.csv")
)

write.csv(
  filter(plot_data, stimulus == 'sucrose') %>%
    select(region, deprived, fed) %>%
    arrange(region),
  paste0(out_path, "/source-data_Fig4a_sucrose.csv")
)

# borboleta responses --------------------------------

borboleta_regions <- c(9, 57, 63)
atlas_labels_plot <-
  gsub(" left", "", atlas_labels[as.character(borboleta_regions)])
plot_data <-
  filter(
    data_peak,
    region %in% borboleta_regions,
    stimulus %ni% c('h2o', 'nostim'),
    MS == 'mated'
  )



stat_wilcox <- plot_data %>%
  filter(MS == "mated") %>%
  group_by(region, stimulus) %>%
  wilcox_test(peak ~ FS) %>%
  add_significance()

# regions with effects
effect_regions <- borboleta_regions

plot_data$MS <- factor(plot_data$MS, levels = c("virgin", "mated"))
plot_data$FS <- factor(plot_data$FS, levels = c("fed", "deprived"))
plot_data$MSFS <- paste0(plot_data$MS, ":", plot_data$FS)
plot_data$MSFS <-
  factor(
    plot_data$MSFS,
    levels = c(
      'virgin:fed',
      'virgin:deprived',
      'mated:fed',
      'mated:deprived'
    )
  )
plot_data$stimulus <-
  factor(plot_data$stimulus,
         levels = c('nostim', 'h2o', 'sucrose', 'yeast'))

for (r in effect_regions) {
  plot_stats <- dplyr::filter(stat_wilcox, region == r)
  if (nrow(plot_stats) != 0)
    plot_stats$y.position <- 1
  plot_data_r  <- filter(plot_data, region == r)
  
  range <- range(plot_data_r$peak, na.rm = T)
  range <- range[2] - range[1]
  step  <- range / 14
  
  if (nrow(plot_stats) != 0) {
    for (s in unique(plot_stats$stimulus)) {
      max_y <- max(filter(plot_data_r, stimulus == s)$peak, na.rm = T)
      pos <- which(plot_stats$stimulus == s)
      plot_stats$y.position[pos] <- seq(length(pos)) * step + max_y
    }
  }
  
  atlas_lab <- atlas_labels_plot[as.character(r)]
  if (is.na(atlas_lab))
    atlas_lab <- ""
  
  p <- ggplot(data = plot_data_r, aes(x = FS, y = peak)) +
    
    geom_boxplot(
      aes(color = FS),
      width = .5,
      outlier.alpha = 0,
      size = .3,
      outlier.shape = 16,
      outlier.size = 3.5,
      outlier.colour = "black"
    ) +
    ggbeeswarm::geom_beeswarm(
      aes(color = FS),
      size = .5,
      shape = 16,
      alpha = .6,
      cex = 4
    ) +
    scale_color_manual(values = fs_colors,
                       'metabolic state',
                       labels = as_labeller(lab_state)) +
    scale_fill_manual(values = fs_colors, 'metabolic state') +
    facet_wrap(stimulus ~ ., nrow = 1, labeller = as_labeller(lab_stim)) +
    labs(subtitle = atlas_lab, y = lab_df) +
    theme_dahaniel2(base_size = 5) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      panel.border = element_blank(),
      legend.position = 'none'
    ) +
    coord_cartesian(clip = "off")
  
  if (nrow(plot_stats) != 0) {
    p <-
      p + stat_pvalue_manual(
        plot_stats,
        label = "{formatC(p, 2)}",
        size = 1.5,
        bracket.size = stat_bracket.size,
        tip.length = stat_tip.length
      )
  }
  plot(p)
  ggsave(
    paste0(out_path, "bp_selected_region", r, ".pdf"),
    width = 22,
    height = 30,
    dpi = plot_dpi,
    scale = plot_scale,
    units = plot_units,
    plot = last_plot()
  )
}
ggsave(
  paste0(out_path, "bp_selected_legend.pdf"),
  width = 40,
  height = 40,
  dpi = 300,
  scale = 3,
  units = 'mm',
  plot = last_plot() + theme(legend.position = 'bottom')
)

for (i in effect_regions) {
  write.csv(
    filter(plot_data, region == i, MS == 'mated') %>%
      select(region, FS, stimulus, peak) %>%
      arrange(region),
    paste0(out_path, paste0("/source-data_FigS4c_", i, ".csv"))
  )
}


# borboleta traces --------------------------------------------------------

plotdata <-
  filter(data,
         region %in% c(9, 57, 63),
         stimulus %in% c('yeast'),
         MS == 'mated') %>%
  group_by(region, frame, stimulus, FS) %>%
  summarize(
    n = n(),
    median = median(dff),
    mean = mean(dff),
    sd = sd(dff),
    sem = sd / sqrt(n)
  ) %>%
  ungroup()


ggplot(plotdata) +
  annotate(
    'segment',
    x = c(0, 10),
    xend = c(5, 20),
    y = -.06,
    yend = -.06,
    color = 'grey65',
    size = 1
  ) +
  geom_ribbon(aes(
    x = frame,
    ymin = mean - sem,
    ymax = mean + sem,
    fill = FS
  ), alpha = .2) +
  geom_line(aes(x = frame, y = mean, color = FS), size = .4) +
  facet_wrap(stimulus ~ region) +
  labs(y = lab_df, x = lab_ts) +
  scale_y_continuous(breaks = c(0, 0.3, .6), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 20, 40), expand = c(0, 0)) +
  scale_color_manual(values = fs_colors) +
  scale_fill_manual(values = fs_colors) +
  theme_dahaniel2(base_size = 5) +
  theme(
    legend.position = 'none',
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(
  paste0(out_path, "borboleta_traces_mean_right.pdf"),
  plot = last_plot(),
  width = 55,
  height = 25,
  scale = plot_scale,
  units = plot_units
)




# response onsets atlas ---------------------------------------------------
data_onsets <- data %>% find_onsets(sd_multi = 4, frame_filter = 10)

data_onsets_sum <- data_onsets %>%
  filter(is.na(onset) == F) %>%
  group_by(state, FS, MS, stimulus, concentration, region) %>%
  summarize(
    sd = sd(onset, na.rm = T),
    n = n(),
    median = median(onset, na.rm = T),
    mad = mad(onset, na.rm = T),
    onset = mean(onset, na.rm = T)
  ) %>%
  filter(n > 4)

range_onset <-
  filter(data_onsets_sum, stimulus == s, MS == ms)$onset %>% range()
for (fs in c('fed', 'deprived')) {
  plot_atlas(
    data = filter(data_onsets_sum, stimulus == "yeast", FS == fs, MS == 'mated'),
    value_column = 'onset',
    plot_legend = T,
    line_width = .1
  ) +
    scale_fill_gradientn(
      colors = RColorBrewer::brewer.pal(9, "Spectral"),
      limits = range_onset,
      values = c(0, .08, .1, .15, .2, .3, .5, .6, .8, 1),
      na.value = '#F7F7F7',
      "time to first response [s]"
    ) +
    theme(legend.position = "bottom")
  ggsave(
    plot = last_plot(),
    file = paste0(out_path, 'response-onset_atlas_', fs, '.pdf'),
    width = 80,
    height = 40,
    scale = plot_scale,
    units = plot_units
  )
}


# response onsets ---------------------------------------------------------
# summarize data for plotting
data_onsets_sum <- data_onsets %>%
  filter(is.na(onset) == F, stimulus == 'yeast') %>%
  group_by(state, FS, MS, stimulus, concentration, region) %>%
  summarize(
    sd = sd(onset, na.rm = T),
    n = n(),
    sem = sd / sqrt(n),
    onset = mean(onset, na.rm = T)
  ) %>%
  filter(n > 4, MS == 'mated') %>%
  group_by(stimulus, region) %>%
  mutate(order_median = max(median),
         order_mean = max(onset)) %>%
  mutate(nstates = n(),
         new = if_else(nstates == 1 & "deprived" %in% FS, TRUE, FALSE)) %>%
  ungroup() %>% arrange(desc(order_mean)) %>%
  mutate(order_index_mean = row_number()) %>%
  arrange(onset) %>%
  mutate(order_all_mean = row_number())

data_onsets_sum$cat <- 'other'
data_onsets_sum$cat[which(data_onsets_sum$region %in% c(58, 59, 25, 26, 68))] <-
  'sensory'
data_onsets_sum$cat[which(data_onsets_sum$region %in% c(19, 44, 71, 73, 74))] <-
  'motor'

data_onsets$cat <- 'other'
data_onsets$cat[which(data_onsets$region %in% c(58, 59, 25, 26, 68))] <-
  'sensory'
data_onsets$cat[which(data_onsets$region %in% c(19, 44, 71, 73, 74))] <-
  'motor'

# shadow plot
data_onsets_sum$shadow <- "no"
data_shadow <-
  data_onsets_sum %>% filter(FS == 'fed', nstates == 2) %>%
  mutate(shadow = "yes", FS = 'deprived')
data_onsets_sum <- rbind(data_onsets_sum, data_shadow)

ggplot(
  filter(data_onsets_sum, stimulus == 'yeast'),
  aes(
    x = onset,
    y = reorder(factor(region), order_index_mean),
    color = cat,
    group = shadow
  )
) +
  annotate(
    'polygon',
    x = c(0, 0, 5, 5) ,
    y = c(-Inf, Inf, Inf, -Inf),
    fill = 'grey55',
    alpha = .2
  ) +
  geom_linerange(aes(
    xmin = onset - sem,
    xmax = onset + sem,
    alpha = shadow
  ), size = .22) +
  geom_point(size = 1, aes(
    alpha = shadow,
    shape = new,
    stroke = .3
  )) +
  geom_text(
    data = filter(data_onsets_sum, shadow == "no"),
    aes(label = n),
    x = 8,
    size = 1.5,
    color = 'grey55'
  ) +
  scale_shape_manual(values = c(16, 1), guide = 'none') +
  scale_alpha_manual(values = c(1, .3), guide = 'none') +
  scale_color_manual(values = pal_cat, '') +
  coord_cartesian(xlim = c(-1, 8)) +
  labs(x = 'time to first response [s]', y = 'atlas region') +
  facet_wrap(~ FS, labeller = as_labeller(lab_state)) +
  theme_dahaniel2() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 0))
ggsave(
  plot = last_plot(),
  file = paste0(out_path, 'F4g_response-onset.pdf'),
  width = 50,
  height = 75,
  scale = plot_scale,
  units = plot_units
)


regions <-
  filter(data_onsets_sum, stimulus == 'yeast')$region %>% unique()

write.csv(
  filter(
    data_onsets,
    stimulus == 'yeast',
    region %in% regions,!is.na(onset),
    MS == 'mated'
  ) %>%
    select(region, FS, onset) %>%
    arrange(FS, region),
  paste0(out_path, "/source-data_Fig4g.csv")
)


# ### deprived vs ff just yeast ------------------------------------------------
library(corrr)
library(corrplot)

data_c <-
  concatenate_ts2(data, c("nostim 0", "h2o 0", "sucrose 200mM", "yeast 10%"))
data_c <- filter(data_c, stimulus == 'yeast')

cor_mated10ds <- data_c %>%
  ungroup() %>%
  filter(state == 'mated_10d_sucrose') %>%
  select(fly, region, frame, dff) %>%
  arrange(frame) %>%
  spread(frame, dff)

cor_matedff <- data_c %>%
  ungroup() %>%
  filter(state == 'mated_ff') %>%
  select(fly, region, frame, dff) %>%
  arrange(frame) %>%
  spread(frame, dff)

all_cors_dep <- data.frame()
for (a in unique(cor_mated10ds$fly)) {
  animal_a <- filter(cor_mated10ds, fly == a, region != 35)
  animal_a_mat <- as.matrix(animal_a[, -c(1, 2)])
  rownames(animal_a_mat) <- animal_a$region
  tidy_cors <- na.omit(t(animal_a_mat)) %>%
    correlate(diagonal = 1) %>%
    as_cordf(diagonal = 1) %>%
    stretch() %>%
    mutate(animal = a)
  all_cors_dep <- rbind(all_cors_dep, tidy_cors)
}

all_cors_dep <- all_cors_dep %>% group_by(x, y) %>%
  summarize(r = mean(r)) %>%
  ungroup() %>%
  retract() %>%
  column_to_rownames("term") %>%
  as.matrix()

pal <-
  colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))(200)

ddist <- as.dist(1 - all_cors_dep)
hc <- hclust(ddist)
regionorder <- hc$labels[hc$order]

all_cors_dep <- all_cors_dep[hc$order, hc$order]

write.csv(all_cors_dep,
          paste0(out_path, "/source-data_Fig4e_deprived.csv"))

svg(
  height = 7,
  width = 7,
  file = paste0(out_path, 'F4_dep_deporder_yeast.svg')
)
cplot <-
  corrplot(
    all_cors_dep,
    type = "full",
    order = "original",
    method = "color",
    addrect = F,
    tl.col = "black",
    tl.srt = 90,
    tl.cex = 1,
    diag = T,
    tl.pos = 'n',
    col = pal,
    title = NULL
  )
dev.off()

png(
  height = 1200,
  width = 1200,
  file = paste0(out_path, 'F4_dep_deporder_yeast.png'),
  type = 'cairo'
)
cplot <-
  corrplot(
    all_cors_dep,
    type = "full",
    order = "original",
    method = "color",
    addrect = F,
    tl.col = "black",
    tl.srt = 90,
    tl.cex = 1,
    diag = T,
    tl.pos = 'n',
    col = pal,
    title = NULL
  )
dev.off()

pdf(
  height = 4,
  width = 4,
  file = paste0(out_path, 'F4_legend.pdf')
)
cplot <-
  corrplot(
    all_cors_dep,
    type = "full",
    order = "original",
    method = "color",
    addrect = F,
    tl.col = "black",
    tl.srt = 90,
    tl.cex = 1,
    diag = T,
    tl.pos = 'n',
    col = pal,
    title = NULL
  )
dev.off()

png(
  height = 1200,
  width = 1200,
  file = paste0(out_path, 'F4_dep_deporder_yeast_Supplement.png'),
  type = 'cairo'
)
cplot <-
  corrplot(
    all_cors_dep,
    type = "full",
    order = "original",
    method = "color",
    addrect = F,
    tl.col = "black",
    tl.srt = 90,
    tl.cex = 1.2,
    diag = T,
    tl.pos = 'full',
    col = pal,
    title = NULL
  )
dev.off()



all_cors_ff <- data.frame()
for (a in unique(cor_matedff$fly)) {
  animal_a <- filter(cor_matedff, fly == a, region != 35)
  animal_a_mat <- as.matrix(animal_a[, -c(1, 2)])
  rownames(animal_a_mat) <- animal_a$region
  tidy_cors <- na.omit(t(animal_a_mat)) %>%
    correlate(diagonal = 1) %>%
    as_cordf(diagonal = 1) %>%
    stretch() %>%
    mutate(animal = a)
  all_cors_ff <- rbind(all_cors_ff, tidy_cors)
}

all_cors_ff <- all_cors_ff %>% group_by(x, y) %>%
  summarize(r = mean(r)) %>%
  ungroup() %>%
  retract() %>%
  column_to_rownames("term") %>%
  as.matrix()

all_cors_ff <- all_cors_ff[hc$order, hc$order]

write.csv(all_cors_ff, paste0(out_path, "/source-data_Fig4e_ff.csv"))

svg(
  height = 7,
  width = 7,
  file = paste0(out_path, 'F4_ff_deporder_yeast.svg')
)
cplot <-
  corrplot(
    all_cors_ff,
    type = "full",
    order = "original",
    method = "color",
    addrect = F,
    tl.col = "black",
    tl.srt = 90,
    tl.cex = 1,
    diag = T,
    tl.pos = 'n',
    col = pal,
    title = NULL
  )
dev.off()

png(
  height = 1200,
  width = 1200,
  file = paste0(out_path, 'F4_ff_deporder_yeast.png'),
  type = 'cairo'
)
cplot <-
  corrplot(
    all_cors_ff,
    type = "full",
    order = "original",
    method = "color",
    addrect = F,
    tl.col = "black",
    tl.srt = 90,
    tl.cex = 1,
    diag = T,
    tl.pos = 'n',
    col = pal,
    title = NULL
  )
dev.off()

png(
  height = 1200,
  width = 1200,
  file = paste0(out_path, 'F4_ff_deporder_yeast_Supplement.png'),
  type = 'cairo'
)
cplot <-
  corrplot(
    all_cors_ff,
    type = "full",
    order = "original",
    method = "color",
    addrect = F,
    tl.col = "black",
    tl.srt = 90,
    tl.cex = 1.2,
    diag = T,
    tl.pos = 'full',
    col = pal,
    title = NULL
  )
dev.off()



# calc difference
all_cors_diff <- all_cors_dep - all_cors_ff

png(
  height = 1200,
  width = 1200,
  file = paste0(out_path, 'F4_diff_deporder_yeast.png'),
  type = "cairo"
)
corrplot(
  all_cors_diff,
  type = "full",
  order = "original",
  method = "color",
  addrect = F,
  tl.col = "black",
  tl.srt = 90,
  tl.cex = 1.2,
  diag = T,
  tl.pos = 'full',
  col = pal,
  title = NULL
)
dev.off()
