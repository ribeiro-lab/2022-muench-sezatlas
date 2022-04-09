# region categories -------------------------------------------------------

atlas_categories <- data.frame(
  region = 1:81,
  category = factor("other", levels = c("other", "input", "output")),
  name = NA
)

atlas_categories$category[c(25, 26, 58, 59, 68, 75, 60)] <- "input"
atlas_categories$category[c(44, 19, 74, 73, 71)] <- "output"

atlas_categories$name[c(25, 26)] <- 'AMS1'
atlas_categories$name[c(58, 59)] <- 'PMS4'
atlas_categories$name[c(68)] <- 'PMS2&3'
atlas_categories$name[19] <- 'MN1'
atlas_categories$name[44] <- 'MN2'
atlas_categories$name[74] <- 'MN3'
atlas_categories$name[73] <- 'MN4'
atlas_categories$name[71] <- 'MN5'

atlas_categories$name <-
  factor(atlas_categories$name,
         levels = c('MN1', 'MN2', 'MN3', 'MN4', 'MN5', 'AMS1', 'PMS4'))

plot_atlas(data = atlas_categories,
           value_column = "name",
           plot_legend = T,
) +
  scale_fill_manual(values = c(mn_cols, sens_cols), na.value = 'grey95') +
  theme(legend.position = "none")
ggsave(
  paste0(out_path, 'map_categories.pdf'),
  width = 200,
  height = 30,
  dpi = plot_dpi,
  scale = plot_scale,
  plot = last_plot(),
  units = 'mm'
)



# cluster -----------------------------------------------------------------

library(lsa)
library(ggdendro)
library(vegan)
library(dendextend)
data_cc <-
  concatenate_ts(filter(data, frame %in% c(-5:45), region != 35))

data_clust <-
  data_cc %>% mutate(group = paste(region, stimulus, MS, FS)) %>%
  filter(state == "mated_10d_sucrose", stimulus != "nostim", region != 35) %>%
  group_by(region, MS, FS, frame) %>%
  summarize(mean = mean(dff, na.rm = TRUE)) %>%
  ungroup()

data_rowmeans <- data_clust %>%
  group_by(region, MS, FS) %>%
  summarize(mean = mean(mean, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(mean) %>%
  mutate(region_index = dplyr::row_number()) %>%
  select(region, region_index) %>%
  ungroup()

data_dist <- data_clust %>% select(region, frame, mean) %>%
  spread(frame, mean) %>%
  column_to_rownames('region')
ddist <- as.dist(1 - cor(t(data_dist)))

hc <- hclust(ddist, "ward.D")

weights <- rep(100, 80)

weights[19] <- 1000000
weights[57] <- 100000
weights[27] <- 10000
weights[15] <- 10000
weights[64] <- 10000
weights[25] <- 1000

hc <- reorder(hc, weights, agglo.FUN = 'mean')
colors <-
  RColorBrewer::brewer.pal(n = 6, name = 'Set2')[c(3, 4, 1, 2, 6, 5)]

dend <- hc %>% as.dendrogram() %>%
  set("branches_k_color", k = 6, value = rev(colors)) %>%
  set("branches_lwd", .5) %>%
  set("labels_cex", .75)

ggdend <- as.ggdend(dend)
ggplot(ggdend, horiz = TRUE, offset_labels = -.3) +
  theme_dahaniel2(base_size = 0) +
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank()
  )

ggsave(
  paste0(
    out_path,
    "ts_all_heatmap_dep_noemptyand35_cor-ward-5-45_dendro.pdf"
  ),
  plot = last_plot(),
  width = 15,
  height = 90,
  scale = 3,
  units = "mm"
)

cluster_order <-
  data.frame(cluster_index = 1:80, region = hc$labels[hc$order])


# plot wit cluster order --------------------------------------------------
data_plot_all <-
  data_cc %>% mutate(group = paste(region, stimulus, MS, FS)) %>%
  filter(state == "mated_10d_sucrose", region != 35) %>%
  group_by(region, MS, FS, fly) %>%
  mutate(mean_for_index = mean(dff, na.rm = TRUE))

data_plot_all <-
  filter(data_plot_all, fly %ni% c('dm_R57C10xGC6s_181127a')) # no sucrose in this animal, skip for plot

data_plot_all <- merge(data_plot_all, cluster_order)

data_plot_all <- data_plot_all %>%
  group_by(region, MS, FS, frame) %>%
  arrange(mean_for_index, by_group = TRUE) %>%
  mutate(
    animal_index = dplyr::row_number(),
    plotindex = paste(cluster_index, animal_index, sep = '_'),
    region_plot_index = paste(region, plotindex, sep = '_')
  )

data_plot_all$plotindex <-
  factor(data_plot_all$plotindex, levels = gtools::mixedsort(unique(data_plot_all$plotindex)))
mx <- max(abs(data_plot_all$dff))

data_plot_all$region <-
  factor(data_plot_all$region, levels = rev(as.character(cluster_order$region)))

p <-
  ggplot(filter(data_plot_all, state == "mated_10d_sucrose", region != 35)) +
  geom_tile(aes(x = frame, y = animal_index, fill = dff), size = 0) +
  facet_grid(region ~ . , switch = "y") +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")),
    limits = c(-mx, mx),
    breaks = c(-1, 0, 1),
    values = c(0, .3, .35, .45, .5, .55, .65, .7, 1),
    lab_df
  ) +
  theme_dahaniel2(base_size = 24) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    strip.text.y.left = element_text(
      size = 10,
      face = 'plain',
      angle = 0,
      color = "grey50"
    ),
    panel.spacing.y = unit(c(0), "lines"),
    legend.position = 'right'
  ) +
  coord_cartesian(expand = FALSE) +
  labs(y = "", x = "")
p
ggsave(
  paste0(
    out_path,
    "ts_all_heatmap_dep_noemptyand35_cor-ward-5-45_facet.png"
  ),
  plot = p + theme(legend.position = 'right'),
  width = 60,
  height = 90,
  scale = 3,
  units = "mm",
  dpi = 600
)
ggsave(
  paste0(
    out_path,
    "ts_all_heatmap_dep_noemptyand35_cor-ward-5-45_facet_stim.png"
  ),
  plot = p + geom_vline(
    xintercept = c(52, 56, 62, 71, 103, 107, 113, 122, 154, 158, 164, 173),
    color = "#8A2BE2",
    alpha = .5
  ),
  width = 60,
  height = 90,
  scale = 3,
  units = "mm",
  dpi = 600
) # indicate stimulus times for alignment


# example traces ----------------------------------------------------------
plotdata <-
  filter(data,
         region %in% c(51, 58, 3),
         stimulus != 'nostim',
         state == "mated_10d_sucrose") %>%
  group_by(region, frame, stimulus) %>%
  summarize(
    n = n(),
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
    y = -.1,
    yend = -.1,
    color = 'grey65',
    size = .4
  ) +
  geom_ribbon(aes(
    x = frame,
    ymin = mean - sem,
    ymax = mean + sem
  ), fill = 'grey85') +
  geom_line(aes(x = frame, y = mean), color = 'grey30', size = .3) +
  facet_grid(region ~ stimulus) +
  labs(y = lab_df, x = lab_ts) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 20, 40), expand = c(0, 0)) +
  theme_dahaniel2(base_size = 5, base_family = "Arial") +
  theme(
    legend.position = 'none',
    panel.border = element_blank(),
    strip.text = element_text(size = 6),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(
  paste0(out_path, "example_traces.pdf"),
  plot = last_plot(),
  width = 35,
  height = 35,
  scale = plot_scale,
  units = plot_units,
  dpi = plot_dpi
)
