## Consolidate every per-fit 2026 FIELD-GROUNDTRUTH metric CSV under [Rr]esults_* into one
## wide table, SPLIT BY APPLICABILITY DOMAIN. The '*-2026gt.csv' files score each fit against
## the common external 2026 field points (vs the fit's own train/test holdout, which
## aggregate_all_metrics.R consolidates), and main's AOA scoring tags each metric row with
## SurfaceArea (Interpolation / Extrapolation, from the presence-referenced AOA) x Distance
## (<=270m / >270m / all). We split into three domains:
##   'Interp <=270m'  - inside AOA, within 270 m of a known occurrence
##   'Interp >270m'   - inside AOA, beyond 270 m
##   'Extrapolation'  - outside AOA (Distance = 'all')
##
## Each domain carries its own `prevalence` column (the field prevalence of that subset), so
## nAP is computed per domain against the real prevalence rather than 1/(1+ratio). Boyce is
## absent (needs presence/background, not the field P/A set).
##
## One row per (fit x domain), columns:
##   tree, resolution, iteration, PAratio, distOrder, seed, fname, domain, n,
##   pr_auc, roc_auc, brier_class, brier_scaled, prevalence,
##   cox_intercept, cox_intercept_se, cox_slope, cox_slope_se
##
## Run with wd = scripts/:  Rscript First_modelling/aggregate_2026gt_metrics.R
## Writes results_summary/all_fit_2026gt.csv.

suppressPackageStartupMessages(library(dplyr))
source('functions.R')

# Use functions.R's PROJ_ROOT if it exists on this machine, else the parent of scripts/ -
# so the aggregator runs on both the data machine and a laptop checkout.
PROJ <- if (dir.exists(PROJ_ROOT)) PROJ_ROOT else normalizePath(file.path(getwd(), '..'))

roots <- list.dirs(PROJ, recursive = FALSE)
roots <- roots[grepl('^[Rr]esults(_|$)', basename(roots))]

# (SurfaceArea, Distance) -> domain label; order sets the factor/plot order (near -> far -> out)
domain_defs <- list(
  list(sa = 'Interpolation', dist = '<=270m', label = 'Interp <=270m'),
  list(sa = 'Interpolation', dist = '>270m',  label = 'Interp >270m'),
  list(sa = 'Extrapolation', dist = 'all',    label = 'Extrapolation')
)

pull <- function(sub, name, col = 'estimate'){
  v <- sub[[col]][sub$metric == name]
  if(length(v) == 0) NA else v[1]
}
metric_row <- function(sub, pf, fname, root, label){
  data.frame(
    tree = basename(root), resolution = pf$resolution, iteration = pf$iteration,
    PAratio = pf$PAratio, distOrder = pf$distOrder, seed = pf$seed, fname = fname,
    domain = label, n = pull(sub, 'pr_auc', 'n'),
    pr_auc = pull(sub, 'pr_auc'), roc_auc = pull(sub, 'roc_auc'),
    brier_class = pull(sub, 'brier_class'), brier_scaled = pull(sub, 'brier_scaled'),
    prevalence = pull(sub, 'prevalence'),
    cox_intercept = pull(sub, 'cox_intercept'), cox_intercept_se = pull(sub, 'cox_intercept', 'std_error'),
    cox_slope = pull(sub, 'cox_slope'), cox_slope_se = pull(sub, 'cox_slope', 'std_error'),
    stringsAsFactors = FALSE)
}

rows <- list()
for(root in roots){
  tables_dir <- file.path(root, 'tables')
  if(!dir.exists(tables_dir)) next
  csvs <- list.files(tables_dir, pattern = '-Seed[0-9]+-2026gt\\.csv$', full.names = TRUE)
  for(csv in csvs){
    fname <- sub('-2026gt$', '', sub('\\.csv$', '', basename(csv)))
    pf <- fname_parse(fname); if(is.na(pf$resolution)) next
    met <- tryCatch(read.csv(csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if(is.null(met) || !'metric' %in% names(met)) next
    if(all(c('SurfaceArea', 'Distance') %in% names(met))){
      for(d in domain_defs){
        sub <- met[met$SurfaceArea == d$sa & met$Distance == d$dist, , drop = FALSE]
        if(nrow(sub)) rows[[length(rows) + 1]] <- metric_row(sub, pf, fname, root, d$label)
      }
    } else {
      # older file without the AOA tags: one undifferentiated row
      rows[[length(rows) + 1]] <- metric_row(met, pf, fname, root, 'all')
    }
  }
}
out <- do.call(rbind, rows)

out_dir <- file.path(PROJ, 'results_summary')
dir.create(out_dir, showWarnings = FALSE)
out_path <- file.path(out_dir, 'all_fit_2026gt.csv')
write.csv(out, out_path, row.names = FALSE)

message(sprintf('wrote %d rows (fit x domain) across %d trees to %s',
                nrow(out), dplyr::n_distinct(out$tree), out_path))
print(dplyr::count(out, tree, domain, name = 'fits') |> as.data.frame())
