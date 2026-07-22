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
# At 1-3arc (~10m) the full study-area domain (~330M cells, spanning both
# the CB and Cochetopa Dome populations ~50km apart) exceeds the Matrix
# package's ~2.1B nonzero-entry sparse-matrix cap when gdistance::transition()
# builds the adjacency matrix — this is a hard 32-bit index limit, not just
# a memory/time problem. Since DesignGroundTruth.R already discards any
# candidate more than euclid_cap_m (2400m) from a known population, distances
# far from occupied patches are never used downstream anyway. So instead of
# building one transition over the full domain, `clusterRegionExtents()`
# splits occupied-patch border points into k geographically disjoint clusters
# (kmeans) and each region is cropped to its cluster's bbox + a 3000m buffer
# (safely past euclid_cap_m). Transitions are built and cost distances
# computed per region, then mosaicked back together with terra::merge().
# The mosaic only spans the clustered regions, not the full suitability grid
# (this bites hardest at 1-3arc, where clustering is doing the most work), so
# each written raster is extend()ed onto the matching -Pr.tif's exact grid
# before writing — a pad, not a resample, since they share origin/res.

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
  suit_full <- terra::rast(suit_path)   # full-domain template cost tiles are re-aligned to below
  patches_r <- terra::rast(file.path(p2res, 'patches', paste0(v, 'patches.tif')))
  occ_tbl   <- read.csv(file.path(p2res, 'patch_summaries', paste0(v, 'occupied.csv')))

  # Extract occupied-patch border points once per threshold — reused both
  # for region clustering below and as leastCostDist/leastCostDistAniso
  # sources in the main loop.
  border_pts_by_thresh <- setNames(lapply(thresh_names, function(th) {
    occ_ids <- occ_tbl$patch[occ_tbl$threshold == th]
    if (length(occ_ids) == 0) return(NULL)
    getOccupiedBorderPts(patches_r[[th]], occ_ids)
  }), thresh_names)

  all_border_pts <- dplyr::bind_rows(Filter(Negate(is.null), border_pts_by_thresh))

  if (nrow(all_border_pts) == 0) {
    message('  no occupied patches at any threshold, skipping.')
    rm(veg, forest, shrub, aoa_r, patches_r, occ_tbl, border_pts_by_thresh, all_border_pts)
    gc()
    next
  }

  # Split into geographically disjoint regions (CB vs. Cochetopa Dome), and
  # recursively bisect any region still too big for gdistance::transition()
  # to build safely — see clusterRegionExtents() header for why (a ~94M-cell
  # region crashed R with a 115GB+ memory spike). 8M cells confirmed safe on
  # this 125GB-RAM machine; the cap here (12M) is untested but still >7x
  # under the crash size — watch RSS on the first run at this setting, since
  # the 94M-cell run held flat until it suddenly blew up rather than
  # degrading gradually.
  region_exts  <- clusterRegionExtents(all_border_pts, k = 2, buffer_m = 5000,
                                        max_cells = 12e6, res_m = mean(terra::res(patches_r)))
  region_iso   <- vector('list', length(region_exts))
  region_aniso <- vector('list', length(region_exts))

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

    sources <- border_pts_by_thresh[[thresh_name]]
    if (is.null(sources) || nrow(sources) == 0) {
      message('  ', thresh_name, ': no occupied patches recorded, skipping.')
      next
    }

    iso_tiles   <- list()
    aniso_tiles <- list()

    for (i in seq_along(region_exts)) {

      reg_ext <- region_exts[[i]]
      src_reg <- sourcesInExt(sources, reg_ext)
      if (nrow(src_reg) == 0) next

      # Build this region's transitions lazily, on first threshold that
      # actually needs them, and reuse across the remaining thresholds.
      if (is.null(region_iso[[i]])) {
        message('  region ', i, ': building resistance surfaces...')
        iso_resist <- buildCostSurface(
          slope = terra::crop(terra::rast(slope), reg_ext),
          tri   = terra::crop(terra::rast(tri), reg_ext),
          suitability = terra::crop(terra::rast(suit_path), reg_ext),
          forest = terra::crop(forest, reg_ext),
          shrub  = terra::crop(shrub, reg_ext),
          weights = c(slope = 1, tri = 1, suitability = 1, forest = 1, shrub = 0.5)
        ) |> terra::mask(terra::crop(aoa_r, reg_ext), maskvalues = 0)

        base_resist <- buildCostSurface(
          slope = terra::crop(terra::rast(slope), reg_ext),
          tri   = terra::crop(terra::rast(tri), reg_ext),
          suitability = terra::crop(terra::rast(suit_path), reg_ext),
          forest = terra::crop(forest, reg_ext),
          shrub  = terra::crop(shrub, reg_ext),
          weights = c(slope = 0, tri = 1, suitability = 1, forest = 1, shrub = 0.5)
        ) |> terra::mask(terra::crop(aoa_r, reg_ext), maskvalues = 0)

        message('  region ', i, ': building isotropic transition...')
        region_iso[[i]] <- buildIsotropicTransition(resistance = iso_resist)
        rm(iso_resist); gc()

        message('  region ', i, ': building anisotropic transition...')
        region_aniso[[i]] <- buildAnisotropicTransition(
          dem = terra::crop(terra::rast(dem), reg_ext), resistance = base_resist)
        rm(base_resist); gc()
      }

      message('  ', thresh_name, ': region ', i, ' (', nrow(src_reg), ' source pts)...')
      # gdistance::transition() drops the CRS itself (see leastCostDist()),
      # so it's passed in from suit_full (shared CRS across the pipeline)
      # and reapplied to the output rather than recovered from the transition.
      iso_tiles[[length(iso_tiles) + 1]]     <- leastCostDist(region_iso[[i]], src_reg, crs = terra::crs(suit_full))
      aniso_tiles[[length(aniso_tiles) + 1]] <- leastCostDistAniso(region_aniso[[i]], src_reg, crs = terra::crs(suit_full))
    }

    if (!file.exists(out_iso) && length(iso_tiles) > 0) {
      message('  ', thresh_name, ': mosaicking isotropic cost distance...')
      # mosaicCostTiles(): overlapping regions resolved by min cost (not
      # tile order), rescaled to 0-65535 once globally (not per-tile), then
      # extended onto suit_full's grid (region clustering only shrinks the
      # extent to the occupied-patch neighbourhoods, not the full domain).
      mosaicCostTiles(iso_tiles, suit_full) |>
        terra::writeRaster(out_iso, datatype = 'INT2U', overwrite = TRUE)
    }

    if (!file.exists(out_aniso) && length(aniso_tiles) > 0) {
      message('  ', thresh_name, ': mosaicking anisotropic cost distance...')
      mosaicCostTiles(aniso_tiles, suit_full) |>
        terra::writeRaster(out_aniso, datatype = 'INT2U', overwrite = TRUE)
    }

    message('  ', thresh_name, ': done.')
  }

  rm(region_iso, region_aniso, border_pts_by_thresh, all_border_pts,
     veg, forest, shrub, aoa_r, suit_full, patches_r, occ_tbl)
  gc()
}
