## Choose the "far" (landscape-scale) presence:absence ratio per resolution, from holdout
## data only - the 2026 field ground truth is kept out of selection entirely, so it stays
## an independent test of the chosen models.
##
## "Far" is scored on holdout presences vs. BACKGROUND absences only (field absences are
## the near-population question), reweighted to a realistic landscape prevalence.
## Done separately per CV design:
##   classic - each model's own holdout (Results_random_iter1a_high_pa/test_data),
##             field absences dropped by matching cell coordinates; 21 ratios x 10 shuffles.
##             SELECTS the ratio.
##   spatial - the common test set, background only (score_common_testset.R); 21 ratios x
##             3 seeds. Confirms. (1/3 arc-second = leave-Cochetopa-out.)
##   lopo    - presences vs background per held-out population; superseded by
##             the LMPO folds elsewhere in the pipeline.
##
## Selection rule, per resolution and design (fixed before looking at any results):
##   criterion   - mean scaled Brier at the target prevalence (a proper scoring rule)
##   constraints - median Cox intercept within +/- max_abs_intercept of 0, and mean ROC-AUC
##                 no more than roc_tolerance below the best ratio's
##   one-SE rule - the SMALLEST ratio passing the constraints whose mean is within one
##                 standard error (across replicates) of the best mean
##   model       - the replicate at that ratio whose scaled Brier is closest to the median
## Repeated at target prevalences 1-5% to show how stable the choice is; 1% is primary.
##
## Usage (working directory scripts/): Rscript First_modelling/select_far_ratio.R
##
## OUTCOME (2026-10-09): chosen ratios at 1% - classic 1:30-1:35, spatial 1:16-1:45 - but
## the choice tracks the target prevalence almost one-for-one (5% -> ~1:5-1:8), i.e. it just
## re-expresses the assumed prevalence. Feeds near_far_evaluation.R, which found nothing of
## use. Kept for the record; not used in the manuscript.

suppressPackageStartupMessages({
  library(sf)
  library(tidyverse)
  library(terra)
  library(ranger)
})
source('functions.R')

primary_target <- 0.01
targets <- c(0.01, 0.02, 0.03, 0.04, 0.05)
max_abs_intercept <- 0.5
roc_tolerance <- 0.02
resolutions <- c('3arc' = 90, '1arc' = 30, '1-3arc' = 10)
out_dir <- file.path(PROJ_ROOT, 'results_highPA_fieldabs', 'tables')

## 1. Classic: score every model's own holdout, background absences only ----------------
classic_root <- file.path(PROJ_ROOT, 'Results_random_iter1a_high_pa')
model_re <- '^(1-3arc|1arc|3arc)-Iteration1-PA1:([0-9.]+)DO:0-Seed(50[0-9]{2})\\.csv$'

field_abs_cells <- function(res_str){
  pres_file <- c('3arc' = '90m', '1arc' = '30m', '1-3arc' = '3m')[[res_str]]
  p <- st_read(file.path('..', 'data', 'Data4modelling', paste0(pres_file, '-presence-iter1.gpkg')), quiet = TRUE)
  fa <- p[p$Presenc == 0 & st_is(p, 'POINT'), ]
  # the stack's Longitude/Latitude layers are what modeller() stored in test_data, so
  # matching on them identifies field absences at the cell level
  e <- terra::extract(rastReader(paste0('dem_', res_str), file.path('..', 'data', 'spatial', 'processed'))[[c('Longitude', 'Latitude')]],
                      vect(fa))
  paste(round(e$Longitude, 3), round(e$Latitude, 3))
}

files <- list.files(file.path(classic_root, 'test_data'), pattern = model_re)
classic <- map_dfr(names(resolutions), function(res_str){
  fa_keys <- field_abs_cells(res_str)
  f_res <- files[startsWith(files, paste0(res_str, '-Iteration'))]
  message(sprintf('classic %s: scoring %d models', res_str, length(f_res)))
  map_dfr(f_res, function(f){
    m <- regmatches(f, regexec(model_re, f))[[1]]
    test <- read.csv(file.path(classic_root, 'test_data', f))
    is_field <- test$Occurrence == 0 & paste(round(test$Longitude, 3), round(test$Latitude, 3)) %in% fa_keys
    test <- test[!is_field, ]
    # at the floor ratio (field absences only, e.g. 1:1.08 at 3 arc-second) a holdout has
    # no background absences at all - nothing to score "far" on, so the model drops out
    if(length(unique(test$Occurrence)) < 2) return(NULL)
    prob <- predict(readRDS(file.path(classic_root, 'models', sub('\\.csv$', '.rds', f))), data = test)$predictions[, '1']
    map_dfr(targets, function(tp)
      bind_cols(tibble(design = 'classic', resolution = res_str, ratio = as.numeric(m[3]), replicate = as.integer(m[4]),
                       target_prevalence = tp, n_field_dropped = sum(is_field)),
                presenceScores(test$Occurrence, prob, target_prevalence = tp)$wide |>
                  select(brier_scaled, cox_intercept, roc_auc)))
  })
})

## 2. Spatial: the common test set, background only ----------------------------------------
spatial <- map_dfr(names(resolutions), function(res_str){
  read.csv(file.path(out_dir, paste0(res_str, '-common-testset-per-model.csv'))) |>
    filter(design == 'field_abs', test_subset == 'background_only', target_prevalence %in% targets) |>
    transmute(design = 'spatial', resolution = res_str, ratio = round(PAratio, 3), replicate = seed,
              target_prevalence, brier_scaled, cox_intercept, roc_auc)
})

## 3. Selection --------------------------------------------------------------------------
per_ratio <- bind_rows(classic, spatial) |>
  group_by(design, resolution, target_prevalence, ratio) |>
  summarise(n = n(), brier_mean = mean(brier_scaled), brier_se = sd(brier_scaled) / sqrt(n()),
            cox_int_median = median(cox_intercept), roc_mean = mean(roc_auc), .groups = 'drop')

# Constraints first define the admissible ratios; the one-SE rule then applies among them
# (best = best ADMISSIBLE ratio). Taking the best over all ratios instead lets an
# inadmissible best - scaled Brier often peaks where the Cox intercept has overshot -
# leave an empty one-SE window.
select_ratio <- function(d){
  admissible <- d |>
    filter(abs(cox_int_median) <= max_abs_intercept,
           roc_mean >= max(d$roc_mean) - roc_tolerance)
  if(nrow(admissible) == 0)
    return(tibble(best_ratio = NA_real_, best_brier = NA_real_, best_se = NA_real_,
                  chosen_ratio = NA_real_, n_admissible = 0L))
  best <- admissible[which.max(admissible$brier_mean), ]
  within_1se <- filter(admissible, brier_mean >= best$brier_mean - best$brier_se)
  tibble(best_ratio = best$ratio, best_brier = best$brier_mean, best_se = best$brier_se,
         chosen_ratio = min(within_1se$ratio), n_admissible = nrow(admissible))
}

selection <- per_ratio |>
  group_by(design, resolution, target_prevalence) |>
  group_modify(~ select_ratio(.x)) |>
  ungroup()

# median-performing replicate at the chosen ratio, for every target (1% is the primary;
# 5% was added later for surface prediction - predict_far_surfaces.R)
median_model <- selection |>
  filter(!is.na(chosen_ratio)) |>
  pmap_dfr(function(design, resolution, target_prevalence, chosen_ratio, ...){
    r <- bind_rows(classic, spatial) |>
      filter(design == !!design, resolution == !!resolution, target_prevalence == !!target_prevalence,
             ratio == chosen_ratio)
    r[which.min(abs(r$brier_scaled - median(r$brier_scaled))), ] |>
      select(design, resolution, target_prevalence, chosen_ratio = ratio, median_replicate = replicate,
             brier_scaled, cox_intercept, roc_auc)
  })

## 4. LOPO: direction check, 1:4 vs 1:20 ---------------------------------------------------
lopo <- map_dfr(names(resolutions), function(res_str){
  read.csv(file.path(PROJ_ROOT, 'results_highPA_populationsplit_fieldabs', 'tables',
                     paste0(res_str, '-LOPO-by-absence-type.csv'))) |>
    filter(subset == 'background vs presences @1%', !is.na(brier_scaled)) |>
    group_by(resolution, ratio) |>
    summarise(n_pops = n(), brier_mean = weighted.mean(brier_scaled, n_p),
              cox_int_median = median(cox_intercept), roc_mean = weighted.mean(roc_auc, n_p), .groups = 'drop')
})

write.csv(per_ratio, file.path(out_dir, 'far-ratio-per-ratio-scores.csv'), row.names = FALSE)
write.csv(selection, file.path(out_dir, 'far-ratio-selection.csv'), row.names = FALSE)
write.csv(median_model, file.path(out_dir, 'far-ratio-median-models.csv'), row.names = FALSE)
write.csv(lopo, file.path(out_dir, 'far-ratio-lopo-check.csv'), row.names = FALSE)

options(width = 200)
message('\n=== chosen far ratio by target prevalence (one-SE rule + constraints) ===')
selection |>
  select(design, resolution, target_prevalence, chosen_ratio) |>
  pivot_wider(names_from = target_prevalence, values_from = chosen_ratio, names_prefix = 'pi=') |>
  as.data.frame() |> print(row.names = FALSE)
message('\n=== median models at the chosen ratio ===')
print(as.data.frame(median_model), digits = 3, row.names = FALSE)
message('\n=== LOPO direction check (presences vs background @1%) ===')
print(as.data.frame(lopo), digits = 3, row.names = FALSE)
