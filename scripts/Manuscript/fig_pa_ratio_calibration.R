## PA ratio vs. calibration/discrimination on the 2026 field ground truth - high-PA grid,
## field absences kept (highPA_grid.R). Rows: metric; columns: CV design. Points by
## resolution: mean across replicates (classic shuffles, spatial seeds, or LOPO folds),
## +/- 1 SD; lines: linear fit on log ratio through every individual model.
## Note the 1/3 arc-second spatial split holds out all of Cochetopa (leave-Cochetopa-out).
##
## Run from scripts/: Rscript Manuscript/fig_pa_ratio_calibration.R

setwd('scripts/')
library(tidyverse)
source('Manuscript/theme_manuscript.R')

cv_labels <- c(
  'classic, high-PA pool + field abs'       = 'Classic',
  'spatial knndm, high-PA pool + field abs' = 'Spatial (knndm)',
  'LOPO, high-PA pool + field abs'          = 'Leave-one-population-out'
)

d <- read.csv('../results/tables/2026-groundtruth-gridsearch-pointpredict-all.csv') |>
  filter(cv_structure %in% names(cv_labels)) |>
  mutate(
    ratio = as.numeric(sub('1:', '', pa_ratio)),
    resolution = factor(resolution, levels = c('3 arc-second', '1 arc-second', '1/3 arc-second')),
    cv = factor(cv_labels[cv_structure], levels = cv_labels)
  )

metrics <- c(
  brier_scaled  = 'Scaled Brier score',
  cox_intercept = 'Cox intercept',
  cox_slope     = 'Cox slope',
  roc_auc       = 'ROC-AUC'
)

long <- d |>
  pivot_longer(all_of(names(metrics)), names_to = 'metric', values_to = 'value') |>
  mutate(metric = factor(metrics[metric], levels = metrics))

s <- long |>
  group_by(resolution, cv, ratio, metric) |>
  summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), .groups = 'drop')

# reference line per panel: no skill (Brier), perfect calibration (Cox), chance (ROC)
ref <- tibble(metric = factor(metrics, levels = metrics), y = c(0, 0, 1, 0.5))

res_colors <- c('3 arc-second' = '#76B9F4', '1 arc-second' = '#2F7FC1', '1/3 arc-second' = '#0A477D')
res_shapes <- c('3 arc-second' = 16, '1 arc-second' = 17, '1/3 arc-second' = 15)
dodge <- position_dodge(width = 0.08)

p <- ggplot(s, aes(ratio, mean, colour = resolution, shape = resolution)) +
  geom_hline(data = ref, aes(yintercept = y), linetype = 'dashed', colour = 'grey60', linewidth = 0.4) +
  # fit through every individual model, linear in log ratio (x is log10-scaled, so
  # y ~ x here is y ~ log10(ratio)); LOPO has only 1:4 and 1:20, so a line is all it supports
  geom_smooth(data = long, aes(ratio, value, colour = resolution, group = resolution),
              inherit.aes = FALSE, method = 'lm', formula = y ~ x, se = FALSE, linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, linewidth = 0.4, position = dodge) +
  geom_point(size = 1.8, position = dodge) +
  facet_grid(metric ~ cv, scales = 'free_y') +
  scale_x_log10(breaks = c(1, 2, 4, 8, 15, 30, 50)) +
  scale_colour_manual(values = res_colors, name = NULL) +
  scale_shape_manual(values = res_shapes, name = NULL) +
  labs(x = 'Presence : Absence Ratio (1:x)', y = NULL) +
  theme_manuscript()

p

ggsave('Manuscript/figures/fig_pa_ratio_calibration.png', p, width = 7.5, height = 8, dpi = 300)
ggsave('Manuscript/figures/fig_pa_ratio_calibration.pdf', p, width = 7.5, height = 8)
