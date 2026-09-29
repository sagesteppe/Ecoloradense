## Score every high-PA spatial-split model on ONE common test set: the entire test side of
## the cached spatial split (results_highPA/tables/{res}-split.rds) - all test presences
## and every test absence, field and background alike.
##
## Why: each model's own holdout (results_*/test_data/) contains only the absences drawn at
## its ratio, so its composition shifts with the ratio - mostly hard, near-population field
## absences at low ratios, mostly easy distant background at high ones - and high ratios
## partly look better because their test sets are easier. Here every ratio faces identical
## data. Scores are also split by absence type, to show how much of any ratio effect lives
## in the field absences vs. the background.
##
## Leak-free: the split partitions every point (all 18k absences included) once, before any
## PA subsampling; modeller() trains only on split$train_id rows, and a ratio only changes
## how many of those it draws. No spatial-split model - diluted (results_highPA/) or
## field-absence (results_highPA_fieldabs/), which share the split - has seen any test-side
## point. Only split 1's models (seeds 5001-5099) are scored; split 2's used another split.
## Classic-split models can't be scored this way: their split is redrawn per fit, so no
## point set is held out from all of them.
##
## Usage (working directory scripts/): Rscript First_modelling/score_common_testset.R 90

suppressPackageStartupMessages({
  library(sf)
  library(tidyverse)
  library(terra)
  library(ranger)
})
source('functions.R')

args <- commandArgs(trailingOnly = TRUE)
resolution <- as.numeric(args[1])
stopifnot('first argument must be a resolution: 90, 30 or 10' = resolution %in% c(90, 30, 10))
res_str <- res_string(resolution)
target_prevalences <- round(seq(0.01, 0.40, by = 0.01), 2)
roots <- c(diluted = 'results_highPA', field_abs = 'results_highPA_fieldabs')
p2proc <- file.path('..', 'data', 'spatial', 'processed')
p2dat  <- file.path('..', 'data', 'Data4modelling')

## 1. The modelling points, built exactly as highPA_grid.R builds them -------------------
abs <- bind_rows(
  st_read(file.path(p2dat, 'iter1-pa.gpkg'), quiet = TRUE) |> select(Occurrence, ID),
  st_read(file.path(p2dat, 'iter1-pa-highPA-reviewed.gpkg'), quiet = TRUE) |> select(Occurrence, ID)
)
pres_file <- c('90' = '90m', '30' = '30m', '10' = '3m')[[as.character(resolution)]]
x <- st_read(file.path(p2dat, paste0(pres_file, '-presence-iter1.gpkg')), quiet = TRUE) |>
  rename(Occurrence = Presenc) |>
  st_as_sf() |>
  mutate(field_abs = Occurrence == 0)
x <- bind_rows(x, mutate(abs, field_abs = FALSE))
x <- filter(x, st_is(x, 'POINT')) |>
  select(Occurrence, field_abs) |>
  mutate(.id = seq_len(n()))

split <- readRDS(file.path(PROJ_ROOT, 'results_highPA', 'tables', paste0(res_str, '-split.rds')))
stopifnot('split was built on different points than rebuilt here' =
            identical(sort(c(split$train_id, split$test_id)), x$.id))

## 2. Common test set + its predictors, extracted once ------------------------------------
test <- x[x$.id %in% split$test_id, ]
covs <- terra::extract(rastReader(paste0('dem_', res_str), p2proc), vect(test))
keep <- complete.cases(covs)
test <- test[keep, ]; covs <- covs[keep, ]
subsets <- list(
  all             = rep(TRUE, nrow(test)),
  field_abs_only  = test$Occurrence == 1 | test$field_abs,
  background_only = test$Occurrence == 1 | !test$field_abs
)
message(sprintf('%s common test set: %d presences, %d field absences, %d background absences.',
                res_str, sum(test$Occurrence == 1), sum(test$field_abs),
                sum(test$Occurrence == 0 & !test$field_abs)))

## 3. Score every split-1 spatial model ---------------------------------------------------
model_re <- paste0('^', res_str, '-Iteration1-PA1:([0-9.]+)DO:0-Seed(50[0-9]{2})\\.rds$')
models <- map_dfr(names(roots), function(design){
  f <- list.files(file.path(PROJ_ROOT, roots[[design]], 'models'), pattern = model_re, full.names = TRUE)
  m <- regmatches(basename(f), regexec(model_re, basename(f)))
  tibble(design = design, path = f,
         PAratio = as.numeric(map_chr(m, 2)), seed = as.integer(map_chr(m, 3)))
})
message(sprintf('Scoring %d models (%s).', nrow(models),
                paste(names(table(models$design)), table(models$design), sep = ': ', collapse = ', ')))

score_model <- function(design, path, PAratio, seed){
  prob <- predict(readRDS(path), data = covs, type = 'response')$predictions[, '1']
  map_dfr(names(subsets), function(s){
    idx <- subsets[[s]]
    truth <- test$Occurrence[idx]; p <- prob[idx]
    map_dfr(c(NA, target_prevalences), function(tp){
      sc <- presenceScores(truth, p, target_prevalence = if(is.na(tp)) NULL else tp)$wide
      bind_cols(tibble(design = design, PAratio = PAratio, seed = seed, test_subset = s,
                       target_prevalence = tp, n_presence = sum(truth == 1), n_absence = sum(truth == 0)),
                select(sc, -any_of(c('prevalence', 'observed_prevalence'))))
    })
  })
}
per_model <- pmap_dfr(select(models, design, path, PAratio, seed), score_model)

summary_tbl <- per_model |>
  group_by(design, test_subset, target_prevalence, PAratio) |>
  summarise(n_models = n(), across(c(pr_auc, roc_auc, brier_class, brier_scaled, cox_intercept, cox_slope,
                                     starts_with('dca_')),
                                   list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)),
                                   .names = '{.fn}_{.col}'),
            .groups = 'drop')

out_dir <- file.path(PROJ_ROOT, 'results_highPA_fieldabs', 'tables')
write.csv(per_model, file.path(out_dir, paste0(res_str, '-common-testset-per-model.csv')), row.names = FALSE)
write.csv(summary_tbl, file.path(out_dir, paste0(res_str, '-common-testset-summary.csv')), row.names = FALSE)

message('\n=== mean scaled Brier on the common test set (all absences), by target prevalence ===')
summary_tbl |>
  filter(test_subset == 'all', target_prevalence %in% c(0.01, 0.03, 0.05, 0.10, 0.34)) |>
  select(design, PAratio, target_prevalence, mean_brier_scaled) |>
  pivot_wider(names_from = target_prevalence, values_from = mean_brier_scaled, names_prefix = 'pi=') |>
  arrange(design, PAratio) |>
  as.data.frame() |>
  print(digits = 2, row.names = FALSE)
message('Written to ', out_dir)
