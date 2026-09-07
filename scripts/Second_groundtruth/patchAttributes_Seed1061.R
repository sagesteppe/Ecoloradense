# Standalone runner for the patchAttributes() step of IdentifyPatches.Rmd,
# which is normally disabled there (eval = F). Extracted as-is so it can be
# run for the classic-split Seed1061 model without knitting the whole
# notebook. threshold_masks/ now holds both Seed1021 and Seed1061, so input
# is filtered explicitly rather than assumed to be a single row.

library(terra)
library(tidyverse)
library(landscapemetrics)

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
source('functions.R')

target_ver <- '1-3arc-Iteration1-PA1:1DO:0-Seed1061'

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
  relocate(version, .before = 1) |>
  filter(version == target_ver)

rm(p_suit, f_suit, p_thresh, f_aoa, f_thresh)

stopifnot('Expected exactly one row for target_ver' = nrow(input) == 1)

lapply(split(input, f = 1:nrow(input)), patchAttributes)
