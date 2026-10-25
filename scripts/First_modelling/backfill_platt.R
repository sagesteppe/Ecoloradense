## Backfill Platt calibrators onto the W0 high-PA stash.
##
## The Results_*_iter1a_high_pa trees were fit at weight_scheme = 'W0', so modeller()
## wrote no -platt.rds and the GT path uses the identity. Cox slopes on the 2026 GT are
## 0.24-0.45 (overconfident). This script bolts a post-hoc Platt map onto each stashed
## caret model, WITHOUT refitting the published ranger model.
##
## Why not just `plattCalibrate(m$pred$obs, m$pred$\`1\`)`? Because modeller()'s caret call
## does not set classProbs = TRUE (metric = 'Accuracy'), so savePredictions = 'final' only
## stashed hard-class predictions. The caret model DOES carry the inner-CV fold indices
## (m$control$index/indexOut) and the training frame (m$trainingData with .outcome and
## .weights). So we replay the inner CV at m$bestTune only - fit ranger(probability = TRUE)
## on each fold's train rows, predict presence probabilities on its holdout rows, assemble
## OOF probs row-by-row, Platt-calibrate.
##
## Usage (wd = scripts/):
##   Rscript First_modelling/backfill_platt.R              # all high-PA trees
##   Rscript First_modelling/backfill_platt.R <tree> ...
##
## Writes <tree>/modelsTune/<fname>-platt.rds for every model that doesn't already have one.
## Prints a=intercept, b=slope per model.

suppressPackageStartupMessages({
  library(ranger)
})

p2proj <- Find(dir.exists, c('/media/steppe/hdd/EriogonumColoradenseTaxonomy',
                             '~/Documents/Ecoloradense'))
source(file.path(p2proj, 'scripts', 'suitability_functions.R'))  # plattCalibrate()

default_trees <- c('Results_LMPO_iter1a_high_pa',
                   'Results_random_iter1a_high_pa',
                   'Results_spatialknn_iter1a_high_pa',
                   'Results_spatialmapping_iter1a_high_pa')
args  <- commandArgs(trailingOnly = TRUE)
trees <- if(length(args) >= 1) args else default_trees

## For one caret model object, assemble OOF presence-class probs by replaying inner CV at
## bestTune. Returns a data.frame(rowIndex, obs, prob) with one row per holdout prediction
## (a row may appear in multiple folds under spatial_mapping/LMPO; we keep every holdout
## prediction and Platt-calibrate the pooled set, matching caret's own savePredictions
## behaviour).
oof_probs <- function(m){
  td <- m$trainingData
  stopifnot('.outcome' %in% names(td))
  td$.outcome <- factor(as.character(td$.outcome), levels = c('0', '1'))
  # Case weights were only attached for weight_scheme != 'W0'; many W0 trees have no
  # .weights column, which is fine - fall back to unit weights.
  has_cw <- '.weights' %in% names(td)
  if(!has_cw) td$.weights <- 1
  bt <- m$bestTune
  idx_in  <- m$control$index
  idx_out <- m$control$indexOut
  stopifnot(length(idx_in) == length(idx_out), length(idx_in) > 0)
  preds <- vector('list', length(idx_in))
  for(k in seq_along(idx_in)){
    tr <- td[idx_in[[k]], , drop = FALSE]
    ho <- td[idx_out[[k]], , drop = FALSE]
    if(nrow(ho) == 0 || length(unique(tr$.outcome)) < 2) next
    # class weights computed on this fold's train slice (modeller()'s formula)
    wts <- c(1 - sum(tr$.outcome == '0') / nrow(tr),
             1 - sum(tr$.outcome == '1') / nrow(tr))
    fit <- ranger::ranger(
      dependent.variable.name = '.outcome',
      data = tr[, setdiff(names(tr), '.weights'), drop = FALSE],
      mtry = bt$mtry, min.node.size = bt$min.node.size,
      splitrule = as.character(bt$splitrule),
      num.trees = m$finalModel$num.trees %||% 750,
      probability = TRUE, class.weights = wts, case.weights = tr$.weights,
      num.threads = 1, verbose = FALSE)
    p <- predict(fit, data = ho[, setdiff(names(ho), c('.outcome', '.weights')),
                                 drop = FALSE])$predictions[, '1']
    preds[[k]] <- data.frame(rowIndex = idx_out[[k]],
                             obs      = as.character(ho$.outcome),
                             prob     = p,
                             stringsAsFactors = FALSE)
  }
  do.call(rbind, preds)
}

`%||%` <- function(a, b) if(is.null(a)) b else a

for(tree in trees){
  mt_dir <- file.path(p2proj, tree, 'modelsTune')
  if(!dir.exists(mt_dir)){
    message('skip (no modelsTune/): ', tree); next
  }
  fs <- list.files(mt_dir, pattern = '\\.rds$', full.names = TRUE)
  fs <- grep('-(platt|trainDI)\\.rds$', fs, invert = TRUE, value = TRUE)
  message('\n=== ', tree, ' (', length(fs), ' models) ===')
  for(i in seq_along(fs)){
    f <- fs[i]
    out <- sub('\\.rds$', '-platt.rds', f)
    if(file.exists(out)) next
    m <- readRDS(f)
    oof <- tryCatch(oof_probs(m), error = function(e) e)
    if(inherits(oof, 'error')){
      message('  ERR ', basename(f), ': ', conditionMessage(oof)); next
    }
    if(nrow(oof) == 0 || length(unique(oof$obs)) < 2){
      message('  SKIP ', basename(f), ': OOF has <2 classes'); next
    }
    cal <- plattCalibrate(oof$obs, oof$prob)
    saveRDS(cal, out)
    cat(sprintf('%s  a=%.3f  b=%.3f  (n=%d)\n',
                basename(f), cal$coef['intercept'], cal$coef['slope'], nrow(oof)))
  }
}
