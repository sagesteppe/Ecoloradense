library(tidyverse)

metrics <- read.csv('all_fit_metrics.csv') |>
  filter(iteration == 1 & tree != 'results_populationgroupsplit') |>
  select(-cox_intercept_se, -cox_slope_se) |>
  mutate(
    ratio = stringr::str_extract(fname, 'PA1[:][:digit:]+'),
    ratio = stringr::str_extract(ratio, '[:digit:]+$'),
    ratio = as.numeric(ratio)
  ) |>
  as_tibble() |>
  pivot_longer(cols = pr_auc:cox_slope)

cln_names = data.frame(
  short = unique(metrics$name),
  full = c('PR-AUC', 'ROC_AUC', 'Brier Class',
   'Boyce Index', 'Cox Intercept', 'Cox Slope')
)

tree_cv <- read.csv('tree_cv_lookup.csv')


metrics = left_join(metrics, cln_names, by = c('name' = 'short')) |>
  select(-name, name = full)

## bring in the outer_split / inner_cv labels for each tree
metrics = left_join(metrics, tree_cv, by = 'tree') |>
  filter(outer_split != 'LOPO')

## brier score closer to 0 is better
## boyce index close to 1 is better
## cox slope closer to 1 is better
## cox intercept closer to 0 is better

## reference (ideal) value per metric, drawn as a green dashed line in each facet
ref_lines <- data.frame(
  name = c('PR-AUC', 'ROC_AUC', 'Brier Class',
           'Boyce Index', 'Cox Intercept', 'Cox Slope'),
  yintercept = c(1, 1, 0, 1, 0, 1)
)

ggplot(metrics, aes(x = ratio, y=value, group = resolution, colour = resolution)) +
  geom_hline(data = ref_lines, aes(yintercept = yintercept),
             colour = 'green', linetype = 'dashed', inherit.aes = FALSE) +
  geom_smooth() +
  facet_wrap(~name,scales = 'free_y') +
  theme_bw() +
  theme(strip.text = element_text(size = 14))


ggplot(metrics, aes(x = ratio, y=value, group = resolution, colour = resolution)) +
  geom_point(alpha = 0.05) +
  geom_hline(data = ref_lines, aes(yintercept = yintercept),
             colour = 'green', linetype = 'dashed', inherit.aes = FALSE) +
  geom_smooth() +
  facet_grid(name ~ outer_split + inner_cv, scales = 'free_y') +
  ## custom y-limits per metric row; oob_keep zooms the view (coord-style)
  ## without dropping points, so geom_smooth is still fit on all data
  ggh4x::facetted_pos_scales(y = list(
    name == 'PR-AUC'        ~ scale_y_continuous(limits = c(0.5, 1),  oob = scales::oob_keep),
    name == 'ROC_AUC'       ~ scale_y_continuous(limits = c(0.7, 1),  oob = scales::oob_keep),
    name == 'Brier Class'   ~ scale_y_continuous(limits = c(0, 0.2),  oob = scales::oob_keep),
    name == 'Boyce Index'   ~ scale_y_continuous(limits = c(-0.5, 1), oob = scales::oob_keep),
    name == 'Cox Intercept' ~ scale_y_continuous(limits = c(-1, 1), oob = scales::oob_keep),
    name == 'Cox Slope'     ~ scale_y_continuous(limits = c(0, 2),  oob = scales::oob_keep)
  )) +
  theme_bw() +
  theme(strip.text = element_text(size = 14))


## ---------------------------------------------------------------------------
## apples-to-apples CV comparison: 3arc only (1-3arc & 1arc lack adequate data
## for the two new partitioner mappings). metric (rows) x strategy (cols).
## knndm kept only as an out-of-competition diagnostic control.
## ---------------------------------------------------------------------------
metrics_3arc <- metrics |>
  filter(resolution == '3arc') |>
  mutate(
    strategy = case_when(
      outer_split == 'random_stratified'      ~ 'random',
      outer_split == 'knndm_class_stratified' ~ 'knndm (diagnostic)',
      str_detect(producer, 'spatial_mapping') ~ 'spatial_mapping',
      str_detect(producer, 'LMPO')            ~ 'LMPO',
      .default = outer_split
    ),
    strategy = factor(strategy,
      levels = c('random', 'spatial_mapping', 'LMPO', 'knndm (diagnostic)'))
  ) |>
  filter(!is.na(strategy))

ggplot(metrics_3arc, aes(x = ratio, y = value)) +
  geom_point(alpha = 0.05) +
  geom_hline(data = ref_lines, aes(yintercept = yintercept),
             colour = 'green', linetype = 'dashed', inherit.aes = FALSE) +
  geom_smooth() +
  facet_grid(name ~ strategy, scales = 'free_y') +
  ggh4x::facetted_pos_scales(y = list(
    name == 'PR-AUC'        ~ scale_y_continuous(limits = c(0.5, 1),  oob = scales::oob_keep),
    name == 'ROC_AUC'       ~ scale_y_continuous(limits = c(0.7, 1),  oob = scales::oob_keep),
    name == 'Brier Class'   ~ scale_y_continuous(limits = c(0, 0.2),  oob = scales::oob_keep),
    name == 'Boyce Index'   ~ scale_y_continuous(limits = c(-0.5, 1), oob = scales::oob_keep),
    name == 'Cox Intercept' ~ scale_y_continuous(limits = c(-1, 1),   oob = scales::oob_keep),
    name == 'Cox Slope'     ~ scale_y_continuous(limits = c(0, 2),    oob = scales::oob_keep)
  )) +
  labs(x = 'P:A ratio', y = NULL,
       title = 'CV strategy comparison (3arc, iteration 1) -- OWN-HOLDOUT (optimistic)') +
  theme_bw() +
  theme(strip.text = element_text(size = 14))


## ---------------------------------------------------------------------------
## THE TRUTH: every model scored on the common 2026 field groundtruth (n=283,
## prevalence 0.339, identical for all strategies). This is the apples-to-apples
## external benchmark; the own-holdout plots above are each strategy's own test
## set and are optimistic. x = P:A ratio, metric (rows) x strategy (cols).
## ---------------------------------------------------------------------------
gt_long <- read.csv('all_fit_metrics.csv') |>
  filter(iteration == 1, resolution == '3arc') |>
  left_join(tree_cv, by = 'tree') |>
  mutate(
    strategy = case_when(
      outer_split == 'random_stratified'      ~ 'random',
      outer_split == 'knndm_class_stratified' ~ 'knndm (diagnostic)',
      str_detect(producer, 'spatial_mapping') ~ 'spatial_mapping',
      str_detect(producer, 'LMPO')            ~ 'LMPO',
      .default = outer_split
    ),
    strategy = factor(strategy,
      levels = c('random', 'spatial_mapping', 'LMPO', 'knndm (diagnostic)'))
  ) |>
  filter(!is.na(strategy)) |>
  select(PAratio, strategy,
         `ROC-AUC`        = gt2026_roc_auc,
         `PR-AUC`         = gt2026_pr_auc,
         `Brier (scaled)` = gt2026_brier_scaled,
         `Cox Intercept`  = gt2026_cox_intercept,
         `Cox Slope`      = gt2026_cox_slope) |>
  pivot_longer(-c(PAratio, strategy), names_to = 'name', values_to = 'value') |>
  mutate(name = factor(name, levels = c('ROC-AUC', 'PR-AUC', 'Brier (scaled)',
                                         'Cox Intercept', 'Cox Slope')))

ref_gt <- data.frame(
  name = factor(c('ROC-AUC', 'PR-AUC', 'Brier (scaled)', 'Cox Intercept', 'Cox Slope'),
                levels = levels(gt_long$name)),
  yintercept = c(1, 1, 1, 0, 1)
)

ggplot(gt_long, aes(x = PAratio, y = value)) +
  geom_point(alpha = 0.05) +
  geom_hline(data = ref_gt, aes(yintercept = yintercept),
             colour = 'green', linetype = 'dashed', inherit.aes = FALSE) +
  geom_smooth() +
  facet_grid(name ~ strategy, scales = 'free_y') +
  ggh4x::facetted_pos_scales(y = list(
    name == 'ROC-AUC'        ~ scale_y_continuous(limits = c(0.4, 1),  oob = scales::oob_keep),
    name == 'PR-AUC'         ~ scale_y_continuous(limits = c(0.3, 1),  oob = scales::oob_keep),
    name == 'Brier (scaled)' ~ scale_y_continuous(limits = c(-1, 1),   oob = scales::oob_keep),
    name == 'Cox Intercept'  ~ scale_y_continuous(limits = c(-2, 3),   oob = scales::oob_keep),
    name == 'Cox Slope'      ~ scale_y_continuous(limits = c(-0.5, 1.5), oob = scales::oob_keep)
  )) +
  labs(x = 'P:A ratio', y = NULL,
       title = 'Field-validated: 2026 groundtruth common test set (3arc, iteration 1)') +
  theme_bw() +
  theme(strip.text = element_text(size = 14))
