# number of responsive regions --------------------------------------------
load_atlas_labels_v2()
data_peak <-
  peak_response(
    data,
    method = "mean",
    frames = 1:5,
    bgf = -5:0,
    sd_mp = SD_cutoff
  )

data_count <- data_peak %>%
  na.omit() %>%
  group_by(fly, MS, FS, stimulus) %>%
  summarize(sum = sum(response)) %>%
  mutate(MSFS = paste(MS, FS)) %>%
  ungroup()

stat_wilcox <- data_count %>%
  filter(MS == 'mated') %>%
  group_by(stimulus, MS) %>%
  wilcox_test(sum ~ FS) %>%
  add_significance()

data_count$MSFS <-
  factor(data_count$MSFS, levels = rev(unique(data_count$MSFS)))
data_count$stimulus  <-
  factor(data_count$stimulus,
         levels = c('nostim', 'h2o',  'sucrose', 'yeast'))
stat_wilcox$stimulus <-
  factor(stat_wilcox$stimulus,
         levels = c('nostim', 'h2o',  'sucrose', 'yeast'))
data_count$FS <-
  factor(data_count$FS, levels = c('fed', 'deprived'))

ggplot(data = filter(data_count, MS == 'mated'), aes(FS, y = sum)) +
  ggbeeswarm::geom_beeswarm(
    aes(color = FS),
    alpha = .4,
    size = .5,
    shape = 16,
    cex = 4
  ) +
  stat_summary(
    aes(color = FS),
    geom = "errorbar",
    width = .3,
    size = .3,
    fun.data = function(x) {
      return(c(
        y = median(x),
        ymin = as.numeric(quantile(x)[2]),
        ymax = as.numeric(quantile(x)[4])
      ))
    },
    show.legend = FALSE
  ) +
  stat_summary(
    aes(color = FS),
    geom = "crossbar",
    width = .5,
    size = .3,
    fun.data = function(x) {
      return(c(
        y = median(x),
        ymin = median(x),
        ymax = median(x)
      ))
    },
    show.legend = TRUE
  ) +
  facet_grid(. ~ stimulus) +
  scale_color_manual(values = fs_colors, 'metabolic state', labels = lab_state) +
  scale_fill_manual(values = fs_colors, 'metabolic state', labels = lab_state) +
  facet_wrap(stimulus ~ ., nrow = 1, labeller = as_labeller(lab_stim)) +
  labs(y = 'number of\nresponsive regions') + #title = paste("Region", r),
  theme_dahaniel2(base_size = 5, base_family = 'Arial') +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 6),
    legend.position = 'bottom',
    legend.key.size = unit(.3, "line")
  ) +
  stat_pvalue_manual(
    stat_wilcox,
    label = "{formatC(p, 2)}",
    size = stat_size,
    bracket.size = stat_bracket.size,
    tip.length = stat_tip.length,
    y.position = 60
  )
ggsave(
  paste0(out_path, "responsive_regions.pdf"),
  width = 40,
  height = 30,
  dpi = plot_dpi,
  scale = plot_scale,
  plot = last_plot(),
  units = plot_units
)
write.csv(
  filter(data_count, MS == 'mated') %>% select(stimulus, FS, sum) %>% arrange(stimulus, FS),
  paste0(out_path, "/source-data_Fig3a.csv")
)

# SVM ---------------------------------------------------------------------

# predict stimulus per state ----------------------------------------------
states <-  c('mated_ff', 'mated_10d_sucrose')
for (s in states) {
  svmdata <- data_peak %>% ungroup() %>%
    filter(state == s,
           region != 35,
           stimulus %in% c("yeast", "sucrose")) %>%
    mutate(label = paste(stimulus)) %>%
    arrange(fly, region, stimulus, concentration) %>%
    mutate(stim_state = paste(stimulus, state, fly),
           MSFS = paste(MS, FS)) %>%
    select(region, fly, stim_state, MSFS, label, peak) %>%
    spread(key = region, value = peak)
  
  
  flies <- unique(svmdata$fly)
  combs <- combn(flies, 2)
  res <- data.frame()
  for (n in 1:dim(combs)[2]) {
    flies_n <- combs[, n]
    svm_test_data  <- svmdata[svmdata$fly %in% flies_n, ]
    svm_training_data <- svmdata[svmdata$fly %ni% flies_n, ]
    svm_training_matrix <- as.matrix(svm_training_data[, -c(1:4)])
    labels <- factor(svm_training_data$label)
    model <-
      svm(x = svm_training_matrix, y = labels, kernel = "linear")
    
    for (i in 1:nrow(svm_test_data)) {
      true <- svm_test_data$label[i]
      pred <-
        predict(model, svm_test_data[i, -c(1:4)], na.action = na.fail)
      res_i <-
        data.frame(
          fly = svm_test_data$fly[i],
          true = true,
          pred = pred,
          hit = true == pred
        )
      res <- rbind(res, res_i)
    }
  }
  
  confusion_matrix <- as.data.frame(table(res[, 2:3])) %>%
    group_by(true) %>%
    mutate(n = sum(Freq), frct = Freq / n) %>%
    ungroup()
  confusion_matrix$true <-
    factor(confusion_matrix$true, levels = sort(levels(confusion_matrix$true)))
  confusion_matrix$pred <-
    factor(confusion_matrix$pred, levels = sort(levels(confusion_matrix$pred)))
  
  ggplot(data = confusion_matrix,
         mapping = aes(y = pred,
                       x = true)) +
    geom_tile(aes(fill = frct)) +
    geom_text(
      aes(label = scales::percent(frct, 2)),
      vjust = .5,
      hjust = .5,
      size = 1.4
    ) +
    coord_fixed() +
    scale_fill_viridis_c(
      limits = c(0, 1),
      breaks = c(0, 1),
      'correct',
      option = 'D',
      labels = scales::percent
    ) +
    labs(subtitle = lab_internalstates_FS[s], y = 'predicted\nstimulation', x = 'true stimulation') +
    theme_dahaniel2(base_size = 5, base_family = 'Arial') +
    theme(
      axis.text.x = element_text(angle = 0, hjust = .5),
      axis.text.y = element_text(angle = 90, hjust = .5),
      legend.position = 'none'
    )
  ggsave(
    paste0(out_path, '/SVM__predict.stim_', s, '.pdf'),
    width = 24,
    height = 24,
    dpi = plot_dpi,
    scale = plot_scale,
    plot = last_plot(),
    units = plot_units
  )
}

# predict FS per stimulus ----------------------------------------------
stimulus <- as.character(unique(data$stimulus))
for (s in c('sucrose', 'yeast')) {
  svmdata <- data_peak %>% ungroup() %>%
    filter(state %in% c('mated_ff',
                        'mated_10d_sucrose'),
           region != 35,
           stimulus %in% s) %>%
    mutate(label = FS) %>%
    arrange(fly, region, stimulus, concentration) %>%
    mutate(stim_state = paste(stimulus, state, fly)) %>%
    select(region, fly, stim_state, label, peak) %>%
    spread(key = region, value = peak)
  
  
  flies <- unique(svmdata$fly)
  combs <- combn(flies, 4)
  res <- data.frame()
  
  for (n in 1:dim(combs)[2]) {
    flies_n <- combs[, n]
    svm_test_data  <- svmdata[svmdata$fly %in% flies_n, ]
    svm_training_data  <- svmdata[svmdata$fly %ni% flies_n, ]
    svm_training_matrix <- as.matrix(svm_training_data[, -c(1:3)])
    labels <- factor(svm_training_data$label)
    model <-
      svm(x = svm_training_matrix, y = labels, kernel = "linear")
    
    for (i in 1:nrow(svm_test_data)) {
      true <- svm_test_data$label[i]
      pred <-
        predict(model, svm_test_data[i, -c(1:3)], na.action = na.fail)
      res_i <-
        data.frame(
          fly = svm_test_data$fly[i],
          true = true,
          pred = pred,
          hit = true == pred
        )
      res <- rbind(res, res_i)
    }
  }
  
  confusion_matrix <- as.data.frame(table(res[, 2:3])) %>%
    group_by(true) %>%
    mutate(n = sum(Freq), frct = Freq / n) %>%
    ungroup()
  confusion_matrix$true <-
    factor(confusion_matrix$true, levels = rev(sort(levels(
      confusion_matrix$true
    ))))
  confusion_matrix$pred <-
    factor(confusion_matrix$pred, levels = rev(sort(levels(
      confusion_matrix$pred
    ))))
  
  
  p <- ggplot(data = confusion_matrix,
              mapping = aes(y = pred,
                            x = true)) +
    geom_tile(aes(fill = frct)) +
    geom_text(
      aes(label = scales::percent(frct, 2)),
      vjust = .5,
      hjust = .5,
      size = 1.4
    ) +
    coord_fixed() +
    scale_x_discrete(labels = lab_state) +
    scale_y_discrete(labels = lab_state) +
    scale_fill_viridis_c(
      limits = c(0, 1),
      breaks = c(0, 1),
      'correct\nclassifications',
      option = 'D',
      labels = scales::percent
    ) +
    labs(title = s, y = 'predicted\nmetabolic state', x = 'true metabolic state') +
    theme_dahaniel2(base_size = 5, base_family = 'Arial') +
    theme(
      axis.text.x = element_text(angle = 0, hjust = .5),
      axis.text.y = element_text(angle = 90, hjust = .5),
      legend.key.size = unit(.2, "line")
    )
  ggsave(
    filename = paste0(out_path, '/SVM_predict.FS_', s, '.pdf'),
    width = 24,
    height = 30,
    dpi = plot_dpi,
    scale = plot_scale,
    units = plot_units,
    plot = p + theme(legend.position = "none")
  )
}
legend <- get_legend(p + theme(legend.position = 'bottom'))
ggsave(
  plot = legend,
  filename = paste0(out_path, '/SVM_legend.pdf'),
  width = 24,
  height = 24,
  dpi = plot_dpi,
  scale = plot_scale,
  units = plot_units
)


# response probability plots -----------------------------------------------

cutoff <- .3
mv = 1
ms = 'mated'

response_p <-
  data_peak %>%
  filter(MS == ms) %>%
  group_by(region, stimulus, concentration, FS) %>%
  summarize(n = n(),
            sum = sum(response),
            p = sum / n) %>%
  filter(n > 5) %>%
  ungroup()

p_diff <- response_p %>% select(region, FS, stimulus, p) %>%
  spread(key = FS, value = p) %>%
  mutate(
    d_DmF = deprived - fed,
    count = case_when(d_DmF > 0 ~ 1, d_DmF < 0 ~ -1, d_DmF == 0 ~ 0),
    d_DmF = ifelse(abs(d_DmF) < cutoff, 0.0, d_DmF),
    region = factor(region)
  ) %>%
  group_by(region) %>%
  mutate(mean = mean(d_DmF)) %>%
  ungroup() %>%
  arrange(mean) %>%
  mutate(index = row_number())



max_val <- max(abs(range(p_diff$d_DmF, na.rm = TRUE)))
if (mv != 0)
  max_val <- mv

pal <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
cols <- c(pal[1:5], rep(pal[6], 3), pal[7:11])

p_diff$stimulus <-
  factor(p_diff$stimulus, levels = c('nostim', 'h2o', 'sucrose', 'yeast'))

ggplot(filter(p_diff)) +
  geom_tile(aes(
    y = reorder(region, index),
    x = stimulus,
    fill = d_DmF
  ), color = "white") +
  scale_fill_gradientn(
    colours = pal,
    limits = c(-max_val, max_val),
    expression("" * Delta * P[response]),
    breaks = c(-1, -0.3, 0, .3, 1)
  ) +
  coord_equal(ratio = .5) +
  theme_dahaniel2(base_size = 5) +
  scale_y_discrete() +
  scale_x_discrete(labels = lab_stim, position = 'top') +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = .5,
      hjust = 0
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position = "bottom"
  ) +
  labs(x = "", y = 'region')
ggsave(
  paste0(out_path, "response_probability_all.pdf"),
  width = 35,
  height = 81,
  dpi = plot_dpi,
  scale = plot_scale,
  units = plot_units,
  plot = last_plot()
)


plot_atlas(
  filter(p_diff, stimulus == 'yeast'),
  value_column = "d_DmF",
  plot_legend = TRUE,
  line_color = "grey20"
) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
    limits = c(-max_val, max_val),
    expression(ΔP[response])
  ) +
  theme(legend.position = 'none')
ggsave(
  paste0(out_path, 'deltaP_yeast.svg'),
  width = 550 / 300,
  height = 200 / 300,
  dpi = 300,
  scale = 5,
  plot = last_plot()
)

plot_atlas(
  filter(p_diff, stimulus == 'sucrose'),
  value_column = "d_DmF",
  plot_legend = TRUE,
  line_color = "grey20"
) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
    limits = c(-max_val, max_val),
    expression(ΔP[response])
  ) +
  theme(legend.position = 'none')
ggsave(
  paste0(out_path, 'deltaP_sucrose.svg'),
  width = 400 / 300,
  height = 200 / 300,
  dpi = 300,
  scale = 5,
  plot = last_plot()
)

plot_atlas(
  filter(p_diff, stimulus == 'h2o'),
  value_column = "d_DmF",
  plot_legend = TRUE,
  line_color = "grey20"
) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
    limits = c(-max_val, max_val),
    expression(ΔP[response])
  ) +
  theme(legend.position = 'none')
ggsave(
  paste0(out_path, 'deltaP_h2o.svg'),
  width = 400 / 300,
  height = 200 / 300,
  dpi = 300,
  scale = 5,
  plot = last_plot()
)

plot_atlas(
  filter(p_diff, stimulus == 'nostim'),
  value_column = "d_DmF",
  plot_legend = TRUE,
  line_color = "grey20"
) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
    limits = c(-max_val, max_val),
    expression(ΔP[response])
  ) +
  theme(legend.position = 'none')
ggsave(
  paste0(out_path, 'deltaP_nostim.svg'),
  width = 400 / 300,
  height = 200 / 300,
  dpi = 300,
  scale = 5,
  plot = last_plot()
)
