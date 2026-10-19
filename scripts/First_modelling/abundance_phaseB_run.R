## Phase B of abundance modelling: refit the full Phase A model stack
## separately for high-elevation and low-elevation training populations.
##
## Context (see `finding_density_model_extrapolation_limits.md`):
##   Phase A was fit on all training populations pooled. Scored against the
##   2026 field data, the two elevational phases turned out to be on ~10x
##   different count scales (occupied median 6 high vs 60.5 low), and the
##   pooled models earn their complexity only on the low-elev phase. Phase B
##   tests whether fitting separately per phase closes the gap on new
##   populations within each phase.
##
## Phase C (later, not here) will refit on 2024+2026 pooled and report only
## the test-split numbers for the final surface.
##
## Low-elevation set = CD only (Cochetopa Dome, 98 plots, 2793-2890 m,
## mean count 31.5). Every other training population is high-elevation
## (14 pops, 440 plots). Confirmed with user 2026-10-10.
##
## OUTPUT NAMING: phase strings land in `bn` before any -I substring so
## splitData()'s cache strip (functions.R L1878, `gsub('-I.*$', '', bn)`)
## produces distinct cache files per phase - the Phase A cache carries
## indices into the full 538-row df and would misalign on the filtered sets.
##
##   Phase A (pooled, already fit): bn = '1-3arc-Iteration1-PA1:1'
##   Phase B high-elev:             bn = '1-3arc-abundB-highelev'
##   Phase B low-elev:              bn = '1-3arc-abundB-lowelev'
##
## Model files land at results/count_models/models/<bn>-<family>.rds and
## follow the same conventions densityModeller() uses for Phase A.
##
## Run (from scripts/):
##   Rscript First_modelling/abundance_phaseB_run.R

## Package set mirrors the other densityModeller() drivers
## (test_phase1_density_scoped.R, compare_split_strategies.R). functions.R's
## poiss()/tweed()/gbs() call bare `multisession` (not future::multisession),
## so `library(future)` must be attached before densityModeller() runs.
suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(terra); library(brms); library(Metrics)
  library(parsnip); library(workflows); library(xgboost); library(bonsai); library(lightgbm)
  library(caret); library(future); library(CAST); library(gstat); library(twinning)
})

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
## LGBM non-spatial fit ships a ~2.5 GiB workflow (data embedded) to each
## future worker; the default 953.67 MiB cap trips tune_race_anova and the
## tune_grid fallback. XGB non-spat and LGBM spat stay under; only the new
## LGBM non-spat rung needs this.
options(future.globals.maxSize = 8 * 1024^3)
source('functions.R')
source('Census_uncertainty/density_bayes.R')
source('Census_uncertainty/density_conformal.R')

## Shared config: this list is the source of truth for Phase B's phase
## assignment and is read by `validate_phaseB_against_2026.R` too. Add
## populations here (not in two places) if the Low set grows.
LOW_ELEV_POPS <- c('CD')

## Phase A suitability raster: Phase B uses the exact same Pr.SuitHab column
## wrapper() extracts for Phase A, filtered by population. Not a new SDM.
phaseA_pr <- '1-3arc-Iteration1-PA1:1DO:0-Seed1021-Pr.tif'

fp <- file.path('..', 'results', 'count_models')

## 1. Assemble the plot-level df (same shape wrapper() hands densityModeller) --

df <- wrapper(phaseA_pr, return_early = TRUE)
coords <- sf::st_coordinates(df)
df$Longitude <- coords[, 1]
df$Latitude  <- coords[, 2]

undeclared_low <- setdiff(LOW_ELEV_POPS, unique(df$Lctn_bb))
if(length(undeclared_low) > 0){
  stop('LOW_ELEV_POPS references population(s) not in training data: ',
       paste(undeclared_low, collapse = ', '))
}

low_df  <- df[df$Lctn_bb %in% LOW_ELEV_POPS, ]
high_df <- df[!df$Lctn_bb %in% LOW_ELEV_POPS, ]

message(sprintf('Phase B split: LOW %d plots (%s); HIGH %d plots (%d populations: %s).',
                nrow(low_df),  paste(LOW_ELEV_POPS, collapse = ', '),
                nrow(high_df), dplyr::n_distinct(high_df$Lctn_bb),
                paste(sort(unique(high_df$Lctn_bb)), collapse = ', ')))

for(nm in c('Low', 'High')){
  d <- if(nm == 'Low') low_df else high_df
  message(sprintf('  %s set: n=%d, mean=%.2f, median=%.1f, max=%d, zeros=%d (%.0f%%).',
                  nm, nrow(d), mean(d$Prsnc_All), stats::median(d$Prsnc_All),
                  max(d$Prsnc_All), sum(d$Prsnc_All == 0),
                  100 * mean(d$Prsnc_All == 0)))
}

## 2. Fit both phase stacks ------------------------------------------------

## The phase string is inserted into the Phase A bn ahead of any -I so
## splitData's cache strip yields a phase-specific key. See header.
bn_low  <- '1-3arc-abundB-lowelev'
bn_high <- '1-3arc-abundB-highelev'

message('\n== Fitting Phase B LOW-elev stack (bn=', bn_low, ') ==')
low_result  <- densityModeller(low_df,  bn = bn_low,  fp = fp)

message('\n== Fitting Phase B HIGH-elev stack (bn=', bn_high, ') ==')
high_result <- densityModeller(high_df, bn = bn_high, fp = fp)

message('\nPhase B fits complete.')
message('  LOW  outputs: results/count_models/models/', bn_low,  '-*.rds')
message('  HIGH outputs: results/count_models/models/', bn_high, '-*.rds')
message('Next: Rscript Census_uncertainty/validate_phaseB_against_2026.R')
