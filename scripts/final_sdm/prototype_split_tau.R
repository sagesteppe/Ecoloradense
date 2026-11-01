## PROTOTYPE: does loosening the outer-split medoid (split_tau) turn reseeding into real
## train/test replicates for spatial_mapping / LMPO?
##
## Background. The spatial_mapping/LMPO outer split is quasi-deterministic: the extrapolation
## components are always held out, and the only seeded choice is kmed_pick() picking the exact
## argmin medoid per k-means cluster (Partitioner.R, Step 4). So reseeding barely moves the test
## set (~85% overlap). split_tau (added to partition_params) makes kmed_pick a seeded SOFT draw
## toward the cluster centre: tau = 0 is the old exact medoid, tau > 0 lets each seed pick a
## different-but-still-central point.
##
## This script rebuilds the partition for one resolution across tau x seeds and measures:
##   * pairwise Jaccard overlap of the presence TEST sets between seeds (lower = more replicate
##     variation unlocked), and
##   * frac_beyond_p95 (the split's env-stringency, already computed in partition_split) - how
##     much the "test points are env-extreme" guarantee softens as tau rises.
## Pick the smallest tau that drops Jaccard meaningfully below ~0.85 without collapsing
## frac_beyond_p95. NO models are fit here - this only probes the split.
##
## Needs the processed DEM rasters (../data/spatial/processed) + the Data4modelling gpkgs, so
## run it on the machine that has them.
##   Working dir scripts/:  Rscript final_sdm/prototype_split_tau.R [resolution] [flavor] [out.csv]
##     resolution: 90 | 30 | 10        (default 30)
##     flavor:     spatial_mapping | LMPO   (default spatial_mapping)

suppressPackageStartupMessages({
  library(sf); library(tidyverse); library(terra)
})
source('functions.R')                      # PROJ_ROOT + rastReader/res_string (via core_functions.R)
source('final_sdm/Partitioner.R')
source('final_sdm/partition_entrypoints.R')

args       <- commandArgs(trailingOnly = TRUE)
resolution <- if (length(args) >= 1) as.numeric(args[1]) else 30
flavor     <- if (length(args) >= 2) args[2] else 'spatial_mapping'
out_csv    <- if (length(args) >= 3) args[3] else
  file.path('final_sdm', sprintf('prototype_split_tau_%s_%s.csv', resolution, flavor))
stopifnot(resolution %in% c(90, 30, 10),
          flavor %in% c('spatial_mapping', 'LMPO'))

taus  <- c(0, 0.15, 0.25, 0.5, 1.0)
seeds <- 5001:5008                         # 8 seeds so pairwise Jaccard is well estimated
p2dat <- file.path('..', 'data', 'Data4modelling')
p2proc <- file.path('..', 'data', 'spatial', 'processed')

## --- data, exactly as highPA_grid.R builds x --------------------------------------------
abs <- bind_rows(
  st_read(file.path(p2dat, 'iter1-pa.gpkg'), quiet = TRUE) |> select(Occurrence, ID),
  st_read(file.path(p2dat, 'iter1-pa-highPA-reviewed.gpkg'), quiet = TRUE) |> select(Occurrence, ID))
pres_file <- c('90' = '90m', '30' = '30m', '10' = '3m')[[as.character(resolution)]]
x <- st_read(file.path(p2dat, paste0(pres_file, '-presence-iter1.gpkg')), quiet = TRUE) |>
  rename(Occurrence = Presenc) |> st_as_sf() |> mutate(field_abs = Occurrence == 0)
x <- bind_rows(x, mutate(abs, field_abs = FALSE))
x <- filter(x, st_is(x, 'POINT')) |> select(Occurrence, field_abs) |> mutate(.id = seq_len(n()))

rast_dat <- rastReader(paste0('dem_', res_string(resolution)), p2proc)
pres <- x[x$Occurrence == 1, ]
bg   <- x[x$Occurrence == 0, ]
build_fn <- if (flavor == 'LMPO') build_partition_lmpo else build_partition_mapping
message(sprintf('%s @ %s: %d presences, %d background.', flavor, pres_file, nrow(pres), nrow(bg)))

## --- build partitions across tau x seed -------------------------------------------------
runs <- list()
for (tau in taus) for (sd in seeds) {
  pp   <- partition_params(seed = sd, split_tau = tau)   # defaults otherwise (= production)
  part <- build_fn(pres, bg, rast_dat, pp)
  runs[[paste(tau, sd)]] <- list(
    tau = tau, seed = sd, test = sort(part$test_id),
    n_test = length(part$test_id),
    fbp95 = part$split_coverage$frac_beyond_p95)
  message(sprintf('  tau=%.2f seed=%d: n_test=%d  frac_beyond_p95=%.3f',
                  tau, sd, length(part$test_id), part$split_coverage$frac_beyond_p95))
}

## --- summarise: mean pairwise Jaccard + stringency per tau ------------------------------
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
summary_tbl <- do.call(rbind, lapply(taus, function(tau) {
  rs <- Filter(function(r) r$tau == tau, runs)
  tests <- lapply(rs, `[[`, 'test')
  pj <- combn(length(tests), 2, function(ij) jac(tests[[ij[1]]], tests[[ij[2]]]))
  data.frame(tau = tau,
             mean_jaccard = round(mean(pj), 3),      # 1 = identical test sets; lower = more variation
             min_jaccard  = round(min(pj), 3),
             mean_n_test  = round(mean(vapply(rs, `[[`, numeric(1), 'n_test')), 1),
             mean_frac_beyond_p95 = round(mean(vapply(rs, `[[`, numeric(1), 'fbp95')), 3))
}))

options(width = 200)
message('\n=== split_tau sweep (', flavor, ' @ ', resolution, '): test-set overlap vs stringency ===')
print(summary_tbl, row.names = FALSE)
write.csv(summary_tbl, out_csv, row.names = FALSE)
message('\nwrote ', out_csv,
        '\nRead: compare mean_jaccard at tau=0 (how much the split already moves on reseeding,',
        ' via kmeans alone) against higher tau. The smallest tau that drops mean_jaccard to a',
        ' useful spread (say <=0.7) while frac_beyond_p95 stays near its tau=0 value is the',
        ' usable looseness. If tau=0 is already low, reseeding alone may suffice.')
