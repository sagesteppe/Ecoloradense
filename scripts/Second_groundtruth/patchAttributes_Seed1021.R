# Standalone runner for the patchAttributes() step of IdentifyPatches.Rmd,
# which is normally disabled there (eval = F). Extracted as-is so it can be
# run for the newly scoped threshold_masks/ (currently just Seed1021)
# without knitting the whole notebook.

library(terra)
library(tidyverse)
library(landscapemetrics)

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
source('functions.R')

p_suit <- file.path('..', 'results', 'suitability_maps')
f_suit <- file.path(p_suit, list.files(p_suit, pattern = 'Pr'))
f_aoa <- file.path(p_suit, list.files(p_suit, pattern = 'AOA'))

p_thresh <- file.path('..', 'results', 'threshold_masks')
f_thresh <- file.path(p_thresh, list.files(p_thresh))

input <- lapply(
  list(f_suit, f_aoa, f_thresh),
  \(x) data.frame(
    filepath = x,
    version = gsub('-thresh.*$|-Pr.*$|-AOA.*$', '', basename(x)))
  ) %>%
  purrr::reduce(inner_join, by = "version") %>%
  setNames( c('Pr', 'version', 'AOA', 'thresholds')) |>
  relocate(version, .before = 1)

rm(p_suit, f_suit, p_thresh, f_aoa, f_thresh)

stopifnot('Expected exactly one scoped version (threshold_masks/ should only hold Seed1021)' =
            nrow(input) == 1)

lapply(split(input, f = 1:nrow(input)), patchAttributes)
