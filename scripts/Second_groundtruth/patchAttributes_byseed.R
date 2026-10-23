# Standalone runner for the patchAttributes() step of IdentifyPatches.Rmd,
# which is normally disabled there (eval = F). Extracted so it can be run for
# any production seed without knitting the whole notebook. threshold_masks/
# holds multiple seeds; input is filtered to the one requested here.
#
# Pick which fitted model this run targets. Positional CLI arg or SEED env
# var; defaults to the original production seed.

library(terra)
library(tidyverse)
library(landscapemetrics)

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
source('functions.R')

.cli_args <- commandArgs(trailingOnly = TRUE)
SEED <- as.integer(if(length(.cli_args) >= 1) .cli_args[1] else Sys.getenv('SEED', unset = '1021'))
target_ver <- fname_build('1-3arc', 1, 1, 0, SEED)

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
