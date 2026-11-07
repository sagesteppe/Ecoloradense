PROJ_ROOT <- '/home/sagesteppe/Documents/Ecoloradense'

library(tidyverse)
set.seed(27)

## Is INDEPENDENT (2026 field ground-truth) performance predictable from HELD-OUT (test
## split) performance, for each evaluation metric? One fixed-effects model per metric:
##   gt_<m> ~ oos_<m> + domain + resolution + tree          (tree = CV strategy)
## controlling for domain/resolution/strategy, we ask whether the holdout score carries
## unique signal about the GT score (partial R^2 of oos_<m>), and whether that holdout->GT
## slope differs by CV strategy (oos_<m>:tree interaction). Mirrors the original nAP scratch,
## extended to ROC-AUC, Cox slope and Cox intercept.

oos        <- read.csv(file.path(PROJ_ROOT, 'results_summary', 'all_fit_metrics.csv'))
gt_metrics <- read.csv(file.path(PROJ_ROOT, 'results_summary', 'all_fit_2026gt.csv'))

## nAP = (AP - baseline prevalence) / (1 - baseline). OOS (holdout) baseline = the fit's own
## prevalence 1/(1+ratio); GT baseline = the field prevalence of that domain.
nAP <- function(ap, prev) (ap - prev) / (1 - prev)
oos$nAP        <- nAP(oos$pr_auc,        1 / (1 + oos$PAratio))
gt_metrics$nAP <- nAP(gt_metrics$pr_auc, gt_metrics$prevalence)

## metrics to test (same model structure for each). Ideal values noted for reference only.
metrics <- c(nAP = 'nAP', roc_auc = 'ROC-AUC',
             cox_slope = 'Cox slope', cox_intercept = 'Cox intercept')

## three CV designs: classic (random), and the tau split-replicate spatial / LMPO trees.
keep_trees <- c('Results_random_iter1a_high_pa',
                'Results_spatialmapping_splittau_phase2',
                'Results_LMPO_splittau_phase2')

## one OOS row per model (tree + fname); joined onto each of that model's GT rows.
analyse_metric <- function(mcol){
  oos_sel <- setNames(oos[, c('tree', 'fname', mcol)], c('tree', 'fname', 'oos_val'))
  gt_sel  <- gt_metrics[, c('tree', 'fname', 'resolution', 'PAratio', 'seed',
                            'domain', 'calibration', mcol)]
  names(gt_sel)[names(gt_sel) == mcol] <- 'gt_val'

  joined <- merge(gt_sel, oos_sel, by = c('tree', 'fname'), all.x = TRUE) |>
    filter(domain != 'Extrapolation',
           tree %in% keep_trees,
           calibration == 'raw',
           is.finite(gt_val), is.finite(oos_val))

  ## one seed per (tree x resolution x ratio x domain) cell -> avoid pseudoreplication
  dat <- joined |>
    group_by(tree, resolution, PAratio, domain) |>
    slice_sample(n = 1) |>
    ungroup()

  model <- lm(gt_val ~ oos_val + domain + resolution + tree, data = dat)
  d1    <- drop1(model, test = 'F')
  pr2   <- function(t){ ss <- d1[t, 'Sum of Sq']; ss / (deviance(model) + ss) }
  ci    <- confint(model)['oos_val', ]
  int_p <- anova(model, update(model, . ~ . + oos_val:tree))[2, 'Pr(>F)']

  summary_row <- data.frame(
    metric       = unname(metrics[mcol]),
    n            = nobs(model),
    R2           = summary(model)$r.squared,
    slope        = coef(model)['oos_val'],
    ci_lo        = ci[1], ci_hi = ci[2],
    partial_R2   = pr2('oos_val'),
    p_oos        = d1['oos_val', 'Pr(>F)'],
    partial_R2_resolution = pr2('resolution'),
    partial_R2_domain     = pr2('domain'),
    partial_R2_tree       = pr2('tree'),
    interaction_p = int_p,
    row.names = NULL)

  ## carry the exact per-cell data the model used, for the facet plot
  dat$metric <- unname(metrics[mcol])
  list(summary = summary_row, data = dat)
}

fits <- lapply(names(metrics), analyse_metric)
res      <- do.call(rbind, lapply(fits, `[[`, 'summary'))
plot_dat <- do.call(rbind, lapply(fits, `[[`, 'data'))
plot_dat$metric <- factor(plot_dat$metric, levels = unname(metrics))

## --- facet plot: GT vs holdout, one panel per metric -----------------------------------
design_lab <- c(Results_random_iter1a_high_pa = 'Classic',
                Results_spatialmapping_splittau_phase2 = 'Spatial mapping',
                Results_LMPO_splittau_phase2 = 'LMPO')
plot_dat$Strategy <- design_lab[plot_dat$tree]

p <- ggplot(plot_dat, aes(oos_val, gt_val)) +
  geom_point(aes(colour = domain, shape = Strategy), alpha = 0.7, size = 1.8) +
  geom_smooth(method = 'lm', se = TRUE, colour = 'black', linewidth = 0.6) +
  facet_wrap(~ metric, scales = 'free') +
  labs(x = 'Held-out (test split) metric', y = 'Independent 2026 ground-truth metric',
       colour = 'Domain', shape = 'CV strategy',
       title = 'Held-out performance is a poor predictor of independent field performance',
       subtitle = 'One point per tree x resolution x PA-ratio x domain (raw; extrapolation excluded)') +
  theme_bw(base_size = 11) +
  theme(legend.position = 'bottom')

ggsave(file.path(PROJ_ROOT, 'results', 'oos_vs_groundtruth_byMetric.png'),
       p, width = 9, height = 7, dpi = 150)

options(width = 200)
fmtp <- function(p) ifelse(p < 1e-3, formatC(p, format = 'e', digits = 1), formatC(p, format = 'f', digits = 3))
print(transform(res,
                R2 = round(R2, 3), slope = round(slope, 3),
                ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3),
                partial_R2 = round(partial_R2, 3),
                p_oos = fmtp(p_oos),
                partial_R2_resolution = round(partial_R2_resolution, 3),
                partial_R2_domain     = round(partial_R2_domain, 3),
                partial_R2_tree       = round(partial_R2_tree, 3),
                interaction_p = round(interaction_p, 3)),
      row.names = FALSE)
