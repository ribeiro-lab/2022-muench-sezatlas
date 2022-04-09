# load packages -----------------------------------------------------------
out_path = "output_final_codecheck/"

library(RmindPeek)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(e1071)

extrafont::loadfonts()

load_atlas_labels_v2()



# define general parameters -----------------------------------------------


SD_cutoff = 4

# colors
is_colors <-  ggthemes::few_pal(palette = "Medium")(4)
names(is_colors) <- c("fed", "deprived", "mated", "virgin")
fs_colors <-
  c(wesanderson::wes_palette("GrandBudapest1", 4, type = "discrete")[2],
    'grey50')
names(fs_colors) <- c("deprived", "fed")
fs_colors2 <-
  c(wesanderson::wes_palette("GrandBudapest1", 4, type = "discrete")[2],
    'grey70')
names(fs_colors2) <- c("protein\ndeprived", "fully fed")
ms_colors <-
  c(wesanderson::wes_palette("IsleofDogs1", 4, type = "discrete")[1],
    'grey30')
names(ms_colors) <- c("mated", "virgin")

# region category colors
# MNs
mn_cols <- c('#49A7CC', '#8ACFEB')
mn_cols <- colorRampPalette(mn_cols)(5)
# sensory
sens_cols <- c('#FCC765', '#FFE6B8')
sens_cols <- colorRampPalette(sens_cols)(3)

pal_cat <- c(mn_cols[1], sens_cols[1], 'grey55')
names(pal_cat) <- c('motor', 'sensory', 'other')

lab_df <- expression("" * Delta * "F/F")
lab_ts <- 'time [s]'
lab_stim <-
  c(
    'nostim' = 'no stimulus',
    'h2o' = 'water',
    'sucrose' = 'sucrose',
    'yeast' = 'yeast'
  )

lab_state <- c('fed' = 'fully fed', 'deprived' = 'protein deprived')
lab_internalstates <-
  c(
    'virgin_ff' = 'virgin, fully fed',
    'virgin_10d_sucrose' = 'virgin, protein deprived',
    'mated_ff' = 'mated, fully fed',
    'mated_10d_sucrose' = 'mated, protein deprived'
  )
lab_internalstates_FS <-
  c(
    'virgin_ff' = 'fully fed',
    'virgin_10d_sucrose' = 'protein deprived',
    'mated_ff' = 'fully fed',
    'mated_10d_sucrose' = 'protein deprived'
  )
lab_internalstates_MS <-
  c(
    'virgin_ff' = 'virgin',
    'virgin_10d_sucrose' = 'virgin',
    'mated_ff' = 'mated',
    'mated_10d_sucrose' = 'mated'
  )


names(fs_colors) <- c("deprived", "fed")

# plots
plot_dpi   = 300
plot_scale = 1
plot_units = 'mm'

stat_size = 1.5
stat_bracket.size = .15
stat_tip.length = 0

ts_path <- "data/imaging/time_series_main/"

# load data ---------------------------------------------------------------

load_atlas_labels_v2()

data <- load_ts_data(
  path = ts_path,
  ds = c(
    "mated_10d_sucrose",
    "mated_ff",
    "virgin_10d_sucrose",
    "virgin_ff"
  ),
  stim_include = c("sucrose 200mM", "nostim 0", "h2o 0", "yeast 10%")
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
