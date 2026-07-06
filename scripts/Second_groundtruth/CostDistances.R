library(terra)
library(sf)
library(tidyverse)
library(gdistance)   
library(raster)     
library(stplanr)

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
source('functions.R')
# Builds isotropic and anisotropic least-cost distance surfaces for each
# model version at the three target resolutions (3arc, 1arc, 1-3arc).
#
# Outputs (per version × threshold combination):
#   results/cost_distances/{version}-{threshold}-costDist.tif      isotropic
#   results/cost_distances/{version}-{threshold}-costDistAniso.tif anisotropic
#
# The anisotropic surface at 1-3arc (~10m) is computationally expensive;
# gdistance builds a directed sparse matrix over the full AOA-masked domain.
# Consider cropping to a tighter bbox if memory is limiting.

p2proc <- '../data/spatial/processed'
p2res  <- '../results'

dir.create(file.path(p2res, 'cost_distances'), showWarnings = FALSE)

# -------------------------------------------------------------------------
# Enumerate versions for the three target resolutions (exclude 3m)
# Follows the same naming convention as IdentifyPatches.Rmd:
#   gsub('patches.tif', '', list.files(...)) gives version strings with a
#   trailing '-', e.g. '3arc-Iteration1-PA1:3DO:0-'
# -------------------------------------------------------------------------
versions <- gsub('patches.tif', '', list.files(file.path(p2res, 'patches')))
versions <- versions[grepl('^(1-3arc|1arc|3arc)-', versions)]

# Process coarsest to finest so errors surface cheaply
res_order <- c('3arc', '1arc', '1-3arc')
versions  <- versions[order(match(sub('-Iteration.*$', '', versions), res_order))]

# -------------------------------------------------------------------------
# Main loop — one pass per model version
# -------------------------------------------------------------------------
for (v in versions) {

  ver <- sub('-$', '', v)          # e.g. '3arc-Iteration1-PA1:3DO:0'
  res <- sub('-Iteration.*$', '', ver)      # e.g. '3arc'

  message('\n--- ', ver, ' ---')

  thresh_names <- c('spec_sens', 'equal_sens_spec', 'sensitivity')

  # Check whether any outputs are missing before loading anything heavy
  any_missing <- any(vapply(thresh_names, function(th) {
    !file.exists(file.path(p2res, 'cost_distances', paste0(ver, '-', th, '-costDist.tif'))) ||
    !file.exists(file.path(p2res, 'cost_distances', paste0(ver, '-', th, '-costDistAniso.tif')))
  }, logical(1)))

  if (!any_missing) {
    message('  all outputs exist, skipping.')
    next
  }

  # --- load layers only when work remains --------------------------------
  p_dem  <- file.path(p2proc, paste0('dem_', res))
  dem    <- file.path(p_dem, 'dem.tif')
  slope  <- file.path(p_dem, 'geomorphology', 'slope.tif')
  tri    <- file.path(p_dem, 'geomorphology', 'ruggedness.tif')
  veg    <- terra::rast(file.path(p_dem, 'Vegetation.tif'))
  forest <- veg[['Tree']]
  shrub  <- veg[['Shrub']]

  aoa_r     <- terra::rast(file.path(p2res, 'suitability_maps', paste0(ver, '-AOA.tif')))
  suit_path <- file.path(p2res, 'suitability_maps', paste0(ver, '-Pr.tif'))
  patches_r <- terra::rast(file.path(p2res, 'patches', paste0(v, 'patches.tif')))
  occ_tbl   <- read.csv(file.path(p2res, 'patch_summaries', paste0(v, 'occupied.csv')))

  # Transitions are built lazily — only on first threshold that needs them
  iso_trans   <- NULL
  aniso_trans <- NULL

  # --- loop over threshold types ----------------------------------------
  for (thresh_name in thresh_names) {

    out_iso   <- file.path(p2res, 'cost_distances',
                           paste0(ver, '-', thresh_name, '-costDist.tif'))
    out_aniso <- file.path(p2res, 'cost_distances',
                           paste0(ver, '-', thresh_name, '-costDistAniso.tif'))

    if (file.exists(out_iso) && file.exists(out_aniso)) {
      message('  ', thresh_name, ': outputs exist, skipping.')
      next
    }

    # Build transitions on first threshold that actually needs them
    if (is.null(iso_trans)) {
      message('  building resistance surfaces...')
      iso_resist <- buildCostSurface(
        slope = slope, tri = tri, suitability = suit_path,
        forest = forest, shrub = shrub,
        weights = c(slope = 1, tri = 1, suitability = 1, forest = 1, shrub = 0.5)
      ) |> terra::mask(aoa_r, maskvalues = 0)

      base_resist <- buildCostSurface(
        slope = slope, tri = tri, suitability = suit_path,
        forest = forest, shrub = shrub,
        weights = c(slope = 0, tri = 1, suitability = 1, forest = 1, shrub = 0.5)
      ) |> terra::mask(aoa_r, maskvalues = 0)

      message('  building isotropic transition...')
      iso_trans <- buildIsotropicTransition(resistance = iso_resist)
      rm(iso_resist); gc()

      message('  building anisotropic transition (slow at 1-3arc)...')
      aniso_trans <- buildAnisotropicTransition(dem = dem, resistance = base_resist)
      rm(base_resist); gc()
    }

    occ_ids <- occ_tbl$patch[occ_tbl$threshold == thresh_name]
    if (length(occ_ids) == 0) {
      message('  ', thresh_name, ': no occupied patches recorded, skipping.')
      next
    }

    message('  ', thresh_name, ': extracting source points...')
    sources <- getOccupiedBorderPts(patches_r[[thresh_name]], occ_ids)

    if (!file.exists(out_iso)) {
      message('  ', thresh_name, ': computing isotropic cost distance...')
      leastCostDist(iso_trans, sources) |>
        terra::writeRaster(out_iso, datatype = 'INT2U', overwrite = TRUE)
    }

    if (!file.exists(out_aniso)) {
      message('  ', thresh_name, ': computing anisotropic cost distance...')
      leastCostDistAniso(aniso_trans, sources) |>
        terra::writeRaster(out_aniso, datatype = 'INT2U', overwrite = TRUE)
    }

    message('  ', thresh_name, ': done.')
  }

  rm(iso_trans, aniso_trans, veg, forest, shrub, aoa_r, patches_r, occ_tbl)
  gc()
}