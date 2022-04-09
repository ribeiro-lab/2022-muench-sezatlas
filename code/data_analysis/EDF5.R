library(lsa)
library(corrr)
library(corrplot)

pal <-
  colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdBu')))(200)
# inter animal correlations -----------------------------------------------
data_pf <-
  filter(data_peak, MS == 'mated', FS == 'fed') %>% droplevels()
data_pd <-
  filter(data_peak, MS == 'mated', FS == 'deprived') %>% droplevels()

data_pf$flystate <-
  factor(data_pf$fly, levels = sort(unique(data_pf$fly)))
levels(data_pf$flystate) <-
  paste0("fly ", 1:length(unique(data_pf$fly)))

data_pd$flystate <-
  factor(data_pd$fly, levels = sort(unique(data_pd$fly)))
levels(data_pd$flystate) <-
  paste0("fly ", 1:length(unique(data_pd$fly)))

data_p <- rbind(data_pf, data_pd)

# yeast
cordata <-
  data_p %>% filter(stimulus %in% c('yeast'),
                    FS == "deprived",
                    MS == "mated",
                    region != 35) %>%
  select(flystate, region, peak) %>%
  spread(region, peak) %>%
  column_to_rownames("flystate") %>%
  t() %>%
  correlate(diagonal = F) %>%
  column_to_rownames("term") %>%
  as.matrix()

svglite::svglite(
  height = 5,
  width = 5,
  pointsize = 7,
  file = paste0(out_path, 'FigureS5a.svg'),
  fix_text_size = F
)
corrplot(
  cordata,
  type = "upper",
  order = "original",
  hclust.method = "single" ,
  method = "circle",
  addrect = 0,
  tl.col = "black",
  tl.srt = 90,
  tl.cex = 2.8,
  mar = c(1, 1, 1, 1),
  cl.cex = 2.2,
  cl.align.text = 'c',
  cl.ratio = 0.24,
  diag = F,
  col = pal,
  title = NULL,
)
dev.off()
write.csv(cordata , paste0(out_path, "source-data_FigS5a.csv"))


# shuffle
data_p_shuffle <-
  data_p %>% group_by(recording, stimulus, concentration) %>%
  mutate(region = sample(region)) %>%
  ungroup()

cordata <-
  data_p_shuffle %>% filter(stimulus %in% c('yeast'),
                            FS == "deprived",
                            MS == "mated",
                            region != 35) %>%
  select(flystate, region, peak) %>%
  spread(region, peak) %>%
  column_to_rownames("flystate") %>%
  t() %>%
  correlate(diagonal = F) %>%
  column_to_rownames("term") %>%
  as.matrix()

svglite::svglite(
  height = 5,
  width = 5,
  pointsize = 7,
  file = paste0(out_path, 'FigureS5b.svg'),
  fix_text_size = F
)
corrplot(
  cordata,
  type = "upper",
  order = "original",
  hclust.method = "single" ,
  method = "circle",
  addrect = 0,
  tl.col = "black",
  tl.srt = 90,
  tl.cex = 2.5,
  mar = c(1, 1, 1, 1),
  cl.cex = 2,
  cl.align.text = 'c',
  cl.ratio = 0.24,
  diag = F,
  col = pal,
  title = NULL,
)
dev.off()
write.csv(cordata , paste0(out_path, "source-data_FigS5b.csv"))


# nostim
cordata <-
  data_p %>% filter(stimulus %in% c('nostim'),
                    FS == "deprived",
                    MS == "mated",
                    region != 35) %>%
  select(flystate, region, peak) %>%
  spread(region, peak) %>%
  column_to_rownames("flystate") %>%
  t() %>%
  correlate(diagonal = F) %>%
  column_to_rownames("term") %>%
  as.matrix()

svglite::svglite(
  height = 5,
  width = 5,
  pointsize = 7,
  file = paste0(out_path, 'FigureS5c.svg'),
  fix_text_size = F
)
corrplot(
  cordata,
  type = "upper",
  order = "original",
  hclust.method = "single" ,
  method = "circle",
  addrect = 0,
  tl.col = "black",
  tl.srt = 90,
  tl.cex = 2.5,
  mar = c(1, 1, 1, 1),
  cl.cex = 2,
  cl.align.text = 'c',
  cl.ratio = 0.24,
  diag = F,
  col = pal,
  title = NULL,
)
dev.off()
write.csv(cordata , paste0(out_path, "source-data_FigS5c.csv"))



# inter animal correlations across stimuli --------------------------------
data_pf <-
  filter(data_peak, MS == 'mated', FS == 'fed') %>% droplevels()
data_pd <-
  filter(data_peak, MS == 'mated', FS == 'deprived') %>% droplevels()

data_pf$flystate <-
  factor(data_pf$fly, levels = sort(unique(data_pf$fly)))
levels(data_pf$flystate) <-
  paste0("fly ", 1:length(unique(data_pf$fly)), " fed")

data_pd$flystate <-
  factor(data_pd$fly, levels = sort(unique(data_pd$fly)))
levels(data_pd$flystate) <-
  paste0("fly ", 1:length(unique(data_pd$fly)), " deprived")

data_p <- rbind(data_pf, data_pd)

cordata <-
  data_p %>% filter(stimulus %in% c('yeast'), MS == "mated", region %ni% c(35)) %>%
  select(flystate, region, peak) %>%
  spread(region, peak) %>%
  column_to_rownames("flystate") %>%
  t() %>%
  correlate(diagonal = T) %>%
  column_to_rownames("term") %>%
  as.matrix()

svglite::svglite(
  height = 8,
  width = 8,
  pointsize = 7,
  file = paste0(out_path, 'FigureS5e_upper-o.svg'),
  fix_text_size = F
)
corrplot(
  cordata,
  type = "upper",
  order = "hclust",
  hclust.method = "single" ,
  method = "circle",
  addrect = NULL,
  tl.col = "black",
  tl.srt = 90,
  tl.cex = 2.8,
  mar = c(1, 1, 1, 1),
  cl.cex = 2.2,
  cl.align.text = 'c',
  cl.ratio = 0.24,
  diag = F,
  col = pal,
  title = NULL,
)
dev.off()
write.csv(cordata , paste0(out_path, "source-data_FigS5e.csv"))



# p1 vs p2 vs p2 shuffled (peak) ----------------------------------------------------------------
data_peak <-
  peak_response(
    data,
    method = "mean",
    frames = 1:5,
    bg = -5:0,
    sd_mp = SD_cutoff
  )
data_peak2 <-
  peak_response(
    data,
    method = "mean",
    frames = 11:20,
    bg = -5:0,
    sd_mp = SD_cutoff
  )
data_peak_last <-
  peak_response(
    data,
    method = "mean",
    frames = 41:45,
    bg = -5:0,
    sd_mp = SD_cutoff
  )
data_peak2_s <-
  data_peak2 %>% group_by(state, stimulus, fly, FS, MS) %>%
  mutate(region = sample(region)) %>%
  ungroup()

all_correlations <- data.frame()
for (s in c("sucrose", "yeast")) {
  data_p1 <- data_peak %>% filter(stimulus == s)
  data_p2 <- data_peak2 %>% filter(stimulus == s)
  data_p2_s <- data_peak2_s %>% filter(stimulus == s)
  data_p_last <- data_peak_last %>% filter(stimulus == s)
  
  cor_p1 <- data_p1 %>%
    ungroup() %>%
    filter(state == 'mated_10d_sucrose', region %ni% c(35)) %>%
    mutate(stimulusfly = paste0(stimulus, "__", fly)) %>%
    select(stimulusfly, region, peak)
  
  cor_p2 <- data_p2 %>%
    ungroup() %>%
    filter(state == 'mated_10d_sucrose', region %ni% c(35)) %>%
    mutate(stimulusfly = paste0(stimulus, "__", fly),
           peak2 = peak) %>%
    select(stimulusfly, region, peak2)
  
  cor_last <- data_p_last %>%
    ungroup() %>%
    filter(state == 'mated_10d_sucrose', region %ni% c(35)) %>%
    mutate(stimulusfly = paste0(stimulus, "__", fly),
           last = peak) %>%
    select(stimulusfly, region, last)
  
  cor_p2_s <- data_p2_s %>%
    ungroup() %>%
    filter(state == 'mated_10d_sucrose', region %ni% c(35)) %>%
    mutate(stimulusfly = paste0(stimulus, "__", fly),
           peak2s = peak) %>%
    select(stimulusfly, region, peak2s)
  
  cor_p1p2 <- merge(cor_p1, cor_p2)
  cor_p1p2p2s <- merge(cor_p1p2, cor_p2_s)
  cor_p1p2p2s <- merge(cor_p1p2p2s, cor_last)
  
  # correlate peaks within animal
  all_cors <- data.frame()
  for (animal in unique(cor_p1p2p2s$stimulusfly)) {
    cor_data_animalx <- filter(cor_p1p2p2s, stimulusfly == animal)
    cor_animalx_p2 <-
      cor(cor_data_animalx$peak, cor_data_animalx$peak2)
    cor_animalx_p2s <-
      cor(cor_data_animalx$peak, cor_data_animalx$peak2s)
    cor_animalx_last <-
      cor(cor_data_animalx$peak, cor_data_animalx$last)
    
    res_animalx <- data.frame(
      animal = animal,
      correlation_p2 = cor_animalx_p2,
      correlation_p2s = cor_animalx_p2s,
      correlation_last = cor_animalx_last,
      stimulus = s
    )
    all_cors <- rbind(all_cors, res_animalx)
  }
  all_correlations <- rbind(all_correlations, all_cors)
}

all_correlations <-
  all_correlations %>% gather(cor, value, -animal,-stimulus)
label_cor <-
  c(
    'correlation_last' = 'post stimulation',
    'correlation_p2' = 'stimulation 2',
    'correlation_p2s' = 'stimulation 2 (shuffled)'
  )
all_correlations$cor <-
  factor(
    all_correlations$cor,
    levels = c("correlation_p2", "correlation_last", "correlation_p2s")
  )

ggplot(data = all_correlations, aes(x = cor, y = value, color = cor)) +
  geom_boxplot(width = .5,
               outlier.alpha = 0,
               size = .3) +
  ggbeeswarm::geom_beeswarm(
    size = .5,
    shape = 16,
    alpha = .6,
    cex = 5
  ) +
  labs(title = NULL, y = 'Pearson correlation') +
  facet_wrap( ~ stimulus, nrow = 1, labeller = as_labeller(lab_stim)) +
  scale_y_continuous(breaks = c(-.2, 0, .2, .4, .6, .8, 1)) +
  scale_x_discrete(labels = as_labeller(label_cor)) +
  scale_color_manual(
    'correlation between\n1st stimulus and:',
    labels = as_labeller(label_cor),
    values = wesanderson::wes_palette("Darjeeling1", 5)[c(4, 3, 5)]
  ) +
  theme_dahaniel2() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank()
  )
ggsave(
  paste0(out_path, "correlations_S5d.pdf"),
  width = 65,
  height = 40,
  dpi = 300,
  scale = 1,
  units = 'mm',
  plot = last_plot()
)
write.csv(
  all_correlations %>% select(stimulus, cor, value),
  paste0(out_path, "source-data_FigS5d.csv")
)
