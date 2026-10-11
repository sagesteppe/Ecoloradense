## Bootstrap-uncertainty prediction surface for the per-phase arithmetic-mean
## baseline. This baseline won (or tied) on both phase cuts in the 2026 field
## validation (High-phase MAE_occ 3.83 for the mean, Low-phase 45.68) ahead of
## every fitted ML candidate and ahead of brms Phase B's catastrophic GP.
## Including it as a predictive surface - with honest uncertainty bands -
## makes the manuscript's conditional-on-presence inversion operational, not
## just rhetorical.
##
## ESTIMATOR: the per-phase mean here is mean-of-plots (over train_df$Prsnc_All),
## NOT a mean of population means - the same quantity densityModeller()'s
## training-mean null reports and the Phase B per-phase constant rows surface.
## The matching bootstrap unit is therefore the PLOT, not the population.
## Resampling plots with replacement within each phase reproduces that
## estimator's sampling distribution.
##
## PHASE ASSIGNMENT ON THE RASTER: a per-cell DEM lookup. CD (Cochetopa Dome,
## the only Low-phase training population) tops out at 2890 m; the next-lowest
## training population (AP) starts at 2928 m. The 2900 m threshold cleanly
## separates them. Change THRESH_M at the top if you want a different cut.
##
## OUTPUTS (results/count_models/rasters/):
##   - abundB-arithmean-mean.tif       - per-cell point estimate (piecewise
##                                       constant: Low-phase mean on Low cells,
##                                       High-phase mean on High cells)
##   - abundB-arithmean-lower.tif      - 2.5% quantile across bootstrap draws
##   - abundB-arithmean-upper.tif      - 97.5% quantile
##   - abundB-arithmean-phase-mask.tif - 1 = High, 2 = Low (for reference)
## Plus results/tables/abundB-arithmean-bootstrap.csv: per-phase bootstrap
## quantiles, standard error, and n_plots - the numbers the manuscript cites.
##
## Run (from scripts/):
##   Rscript Census_uncertainty/arithmetic_mean_surface.R

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(terra)
})

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
source('functions.R')

## ---- Config ---------------------------------------------------------------

LOW_ELEV_POPS <- c('CD')          # mirrors abundance_phaseB_run.R
THRESH_M      <- 2900             # elevation threshold for raster phase assignment
B             <- 1000             # bootstrap draws
CI_LEVEL      <- 0.95             # for the lower/upper rasters
SEED          <- 1021
bn_phaseA     <- '1-3arc-Iteration1-PA1:1'

set.seed(SEED)

p2proc  <- file.path(PROJ_ROOT, 'data', 'spatial', 'processed')
out_dir <- file.path('..', 'results', 'count_models', 'rasters')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Training data -----------------------------------------------------

train_sf <- wrapper(paste0(bn_phaseA, 'DO:0-Seed1021-Pr.tif'), return_early = TRUE)
train_df <- sf::st_drop_geometry(train_sf)
train_df$Phase <- ifelse(train_df$Lctn_bb %in% LOW_ELEV_POPS, 'Low', 'High')

phase_summary <- train_df |>
  dplyr::group_by(Phase) |>
  dplyr::summarize(n_plots = dplyr::n(),
                   n_pops  = dplyr::n_distinct(Lctn_bb),
                   mean    = mean(Prsnc_All),
                   .groups = 'drop')
message('Training plots by phase:')
message(paste(utils::capture.output(print(as.data.frame(phase_summary), row.names = FALSE)),
              collapse = '\n'))

## ---- 2. Plot-level bootstrap --------------------------------------------

## For each phase, resample its plots with replacement B times and compute
## the mean per draw. Produces a B x 2 matrix of (High_mean, Low_mean) draws.
bootstrap_phase_mean <- function(x, B){
  n <- length(x)
  replicate(B, mean(x[sample.int(n, n, replace = TRUE)]))
}

draws <- data.frame(
  High = bootstrap_phase_mean(train_df$Prsnc_All[train_df$Phase == 'High'], B),
  Low  = bootstrap_phase_mean(train_df$Prsnc_All[train_df$Phase == 'Low'],  B)
)

alpha <- 1 - CI_LEVEL
boot_summary <- tibble::tibble(
  Phase    = c('High', 'Low'),
  n_plots  = c(sum(train_df$Phase == 'High'), sum(train_df$Phase == 'Low')),
  n_pops   = c(dplyr::n_distinct(train_df$Lctn_bb[train_df$Phase == 'High']),
               dplyr::n_distinct(train_df$Lctn_bb[train_df$Phase == 'Low'])),
  estimate = c(mean(train_df$Prsnc_All[train_df$Phase == 'High']),
               mean(train_df$Prsnc_All[train_df$Phase == 'Low'])),
  boot_mean  = c(mean(draws$High), mean(draws$Low)),
  boot_se    = c(stats::sd(draws$High), stats::sd(draws$Low)),
  lower      = c(stats::quantile(draws$High, probs = alpha / 2,     names = FALSE),
                 stats::quantile(draws$Low,  probs = alpha / 2,     names = FALSE)),
  upper      = c(stats::quantile(draws$High, probs = 1 - alpha / 2, names = FALSE),
                 stats::quantile(draws$Low,  probs = 1 - alpha / 2, names = FALSE)),
  B          = B,
  CI_level   = CI_LEVEL
)
cat('\nBootstrap summary (plot-level, B=', B, '):\n', sep = '')
print(as.data.frame(boot_summary), row.names = FALSE)

## ---- 3. Build surfaces ----------------------------------------------------

## DEM over the project's working extent. rastReader returns a SpatRaster
## with all covariate layers; elevation is the one we need for the phase mask.
rast_dat <- rastReader('dem_1-3arc', p2proc)
elev <- rast_dat[['elevation']]

phase_mask <- terra::classify(elev,
  rbind(c(-Inf, THRESH_M, 2),   # Low  -> 2
        c(THRESH_M,  Inf, 1)))  # High -> 1
names(phase_mask) <- 'Phase'    # 1 = High, 2 = Low

make_const_raster <- function(high_val, low_val){
  r <- terra::rast(phase_mask)
  r[] <- ifelse(terra::values(phase_mask) == 1, high_val,
         ifelse(terra::values(phase_mask) == 2, low_val, NA_real_))
  r
}

r_mean  <- make_const_raster(boot_summary$estimate[boot_summary$Phase == 'High'],
                             boot_summary$estimate[boot_summary$Phase == 'Low'])
r_lower <- make_const_raster(boot_summary$lower[boot_summary$Phase == 'High'],
                             boot_summary$lower[boot_summary$Phase == 'Low'])
r_upper <- make_const_raster(boot_summary$upper[boot_summary$Phase == 'High'],
                             boot_summary$upper[boot_summary$Phase == 'Low'])

names(r_mean)  <- 'arithmean'
names(r_lower) <- paste0('arithmean_q', formatC(100 * (alpha / 2),   format = 'g'))
names(r_upper) <- paste0('arithmean_q', formatC(100 * (1 - alpha/2), format = 'g'))

## ---- 4. Write ------------------------------------------------------------

tif <- function(nm) file.path(out_dir, paste0('abundB-arithmean-', nm, '.tif'))
terra::writeRaster(r_mean,     tif('mean'),       overwrite = TRUE)
terra::writeRaster(r_lower,    tif('lower'),      overwrite = TRUE)
terra::writeRaster(r_upper,    tif('upper'),      overwrite = TRUE)
terra::writeRaster(phase_mask, tif('phase-mask'), overwrite = TRUE,
                   datatype = 'INT1U')

out_csv <- file.path('..', 'results', 'tables', 'abundB-arithmean-bootstrap.csv')
write.csv(boot_summary, out_csv, row.names = FALSE)

message('Written:')
message('  ', tif('mean'))
message('  ', tif('lower'))
message('  ', tif('upper'))
message('  ', tif('phase-mask'))
message('  ', out_csv)
