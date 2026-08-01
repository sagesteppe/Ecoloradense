# =============================================================================
# Stratified ground-truth design for SDM verification — REAL DATA
#
# Real-data counterpart to sdm_groundtruth_design.R. Same downstream pipeline
# (region-orthogonalized cost anomaly R, GRTS w/ oversample, cLHS comparison),
# wired to actual rasters instead of the synthetic landscape. Resolution is
# whatever cfg$ver points at (currently 1arc, ~30m) - all messages/filenames
# below are tagged from cfg$ver, not hardcoded to one resolution.
#
# Differences from the synthetic version, and why:
#   - No D (suitability-disagreement) axis. 
#   - "Region" = named Site polygons from GroundTruch_Areas.gpkg, not a
#     synthetic tiling. Sites are unequal in size, so region allocation is
#     proportional-to-size with a floor, not flat n_per_region.
#   - E (Euclidean distance) is pixel-wise, but Cost (least-cost distance) was
#     built in CostDistances.R from occupied PATCH BORDER points, not from
#     individual occurrence points like E. That mismatch is intentional and is
#     exactly what R (the within-region residual of Cost on E) is designed to
#     surface — it is not something to reconcile away.
# =============================================================================

suppressPackageStartupMessages({
  library(terra)      # rasters
  library(sf)          # vector / sample frame
  library(dplyr)       # wrangling + group-wise residualization
  library(mgcv)        # gam smooth for R residualization
  library(spsurvey)    # grts()  (>= 5.0 interface)
  library(clhs)        # cost-constrained conditioned LHS (comparison)
  library(ggplot2)     # coverage diagnostics
})

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts/Second_groundtruth')
set.seed(42)

# -----------------------------------------------------------------------------
# 0.  CONFIG
# -----------------------------------------------------------------------------
cfg <- list(
  ver           = '1-3arc-Iteration1-PA1:2.7DO:0',  # model version
  cost_thresh   = 'spec_sens',                  # threshold used for cost-distance sources

  n_total_base   = 360,    # total base sites across ALL regions (budget, not per-region)
  min_per_region = 5,      # floor: every region gets at least this many base sites
                            # (capped at that region's available pixel count). NOTE:
                            # floor * n_regions must leave room in n_total_base for the
                            # proportional term, or every region just gets the floor.
  max_per_region = 60,     # ceiling: no region gets more than this many base sites,
                            # regardless of how much pixel-count share it would command
                            # (e.g. Cochetopa would otherwise swallow most of the budget).
                            # Excess budget from capped regions is NOT redistributed.
  size_exponent  = 0.5,    # allocate on avail^size_exponent, not raw avail. The 4 biggest
                            # sites hold 88% of all candidate pixels, so linear (=1) leaves
                            # mid-size regions (1-5k px) at just the floor. <1 compresses
                            # that dominance and shifts share toward the mid tier; 1 = linear.
  oversample_mult = 0.5,   # n_over = mult * n_base, capped at remaining available pixels

  w_coverage   = 1.0,      # flat base: guarantees span of S for identifiability
  w_discovery  = 1.0,      # bump in the contested mid-S band (find new occ.)
  w_far        = 0.8,      # reward high-E tail (break the "near" truncation)
  w_costanom   = 1.0,      # reward |R| (where the cost surface earns its keep)

  discovery_center = 0.5,  # mid suitability, on the within-region quantile scale
  discovery_sigma  = 0.18,

  r_cap_q = 0.90,          # stop upweighting cells above this within-region R quantile

  mindis = 30             # min spacing between sites, in CRS units (metres; UTM 13N)
)

# =============================================================================
# 1.  REAL INPUTS
# =============================================================================
p2res  <- '../../results'
p2data <- '../../data/Data4modelling'

suit_path  <- file.path(p2res, 'suitability_maps', paste0(cfg$ver, '-Pr.tif'))
aoa_path   <- file.path(p2res, 'suitability_maps', paste0(cfg$ver, '-AOA.tif'))
aniso_path <- file.path(p2res, 'cost_distances',
                        paste0(cfg$ver, '-', cfg$cost_thresh, '-costDistAniso.tif'))
iso_path   <- file.path(p2res, 'cost_distances',
                        paste0(cfg$ver, '-', cfg$cost_thresh, '-costDist.tif'))

stopifnot(
  'Suitability raster not found' = file.exists(suit_path),
  'AOA raster not found'         = file.exists(aoa_path),
  'No cost-distance raster found for this version/threshold. Run CostDistances.R first.' =
    file.exists(aniso_path) || file.exists(iso_path)
)

res_tag <- sub('-Iteration.*$', '', cfg$ver)   # e.g. "1arc" - used to tag outputs

S_mean <- terra::rast(suit_path)
aoa_r  <- terra::rast(aoa_path)

Cost <- if (file.exists(aniso_path)) {
  message('Using anisotropic cost distance.')
  terra::rast(aniso_path)
} else {
  message('Anisotropic cost distance not available; using isotropic.')
  terra::rast(iso_path)
}
# CostDistances.R now reapplies CRS itself (leastCostDist()/leastCostDistAniso()
# take a crs= arg, since gdistance::transition() drops it internally) before
# writing, so freshly-generated tifs already carry it. This line is now just a
# defensive fallback for tifs written before that fix (extent/res/dims already
# match S_mean exactly, so it's a metadata patch, not a resample).
terra::crs(Cost) <- terra::crs(S_mean)
stopifnot('Cost raster does not align with suitability raster' =
            terra::compareGeom(S_mean, Cost, stopOnError = FALSE))

# --- region: named Site polygons -> one (multi)polygon per site -> rasterize ---
region_v <- sf::st_read(
  file.path('..', '..', 'data', 'hikingTrails', 'GroundTruch_Areas.gpkg'),
  quiet = TRUE)

region_v$region_id <- as.integer(factor(region_v$Site))
region <- terra::rasterize(terra::vect(region_v), S_mean, field = 'region_id')

# --- known occurrences: source of Euclidean distance E ---
### NOTE: uses the full-precision 3m presence record regardless of model
### resolution, same convention as DesignGroundTruth.R.
occ <- sf::st_read(file.path(p2data, '3m-presence-iter1.gpkg'), quiet = TRUE) |>
  dplyr::filter(Presenc == 1) |>
  sf::st_transform(terra::crs(S_mean))

E <- terra::distance(S_mean, terra::vect(occ))

# --- mask everything to AOA (model's inference domain) + region (feasible-to-groundtruth) ---
S_mean <- terra::mask(S_mean, aoa_r, maskvalues = 0) |> terra::mask(region)
E      <- terra::mask(E,      aoa_r, maskvalues = 0) |> terra::mask(region)
Cost   <- terra::mask(Cost,   aoa_r, maskvalues = 0) |> terra::mask(region)

# --- patch: CostDistances.R only routes cost distance through cells within
# reach of known occurrences, so regions with no nearby occurrences (e.g.
# Upper Crossing) come back all-NA on Cost. That silently drops the whole
# region once as.data.frame(na.rm=TRUE) runs. Backfill with the max observed
# cost - these are, by construction, at least as hard to reach as anywhere
# already in the surface.
cost_max      <- terra::global(Cost, 'max', na.rm = TRUE)[1, 1]
cost_gap_mask <- is.na(Cost) & !is.na(region)
n_patched     <- terra::global(cost_gap_mask, 'sum', na.rm = TRUE)[1, 1]
if (n_patched > 0) {
  message(sprintf('Backfilling %s Cost-NA pixels within valid regions to max cost (%.1f).',
                  format(n_patched, big.mark = ','), cost_max))
  Cost <- terra::ifel(cost_gap_mask, cost_max, Cost)
}

# =============================================================================
# 2.  Assemble the per-pixel table over the mask
# =============================================================================
stk <- c(S_mean, E, Cost, region)
names(stk) <- c('S', 'E', 'Cost', 'region')

df <- as.data.frame(stk, cells = TRUE, xy = TRUE, na.rm = TRUE) |>
  dplyr::mutate(region = factor(region, labels = levels(factor(region_v$region_id)))) |>
  dplyr::filter(!is.na(region))
# reattach human-readable Site names for messaging/QA
site_lkp <- setNames(region_v$Site, region_v$region_id)
df$Site  <- site_lkp[as.character(df$region)]

message(sprintf('Pixels in frame: %s across %d regions (Sites)',
                format(nrow(df), big.mark = ','), nlevels(df$region)))

# =============================================================================
# 4.  Region-aware cost orthogonalization:  R = within-region residual of C on E
#     (same logic/guardrails as the synthetic script; tiny regions fall back
#     to a simple scaled Cost-E difference rather than a gam/lm smooth)
# =============================================================================
df <- df |>
  dplyr::group_by(region) |>
  dplyr::group_modify(~{
    if (nrow(.x) >= 30 && length(unique(.x$E)) > 10) {
      mm <- tryCatch(gam(Cost ~ s(E, k = 6), data = .x, method = 'REML'),
                     error = function(e) lm(Cost ~ poly(E, 2), data = .x))
      .x$R <- resid(mm)
    } else .x$R <- as.numeric(scale(.x$Cost - .x$E))
    .x
  }) |>
  dplyr::ungroup()

# =============================================================================
# 5.  Within-region standardization (rank/quantile)
# =============================================================================
qn <- function(v) (rank(v, ties.method = 'average') - 0.5) / length(v)
df <- df |>
  dplyr::group_by(region) |>
  dplyr::mutate(Sq = qn(S), Eq = qn(E), Rq = qn(R)) |>
  dplyr::ungroup()


# =============================================================================
# 6.  PRIORITY SURFACE — S/E/R only (no D term; nothing to disagree with)
# =============================================================================
build_priority <- function(d, k) {
  disc <- exp(-((d$Sq - k$discovery_center)^2) / (2 * k$discovery_sigma^2))
  far  <- d$Eq
  costanom <- abs(d$Rq - 0.5) * 2
  costanom[d$Rq > k$r_cap_q] <- abs(k$r_cap_q - 0.5) * 2

  p <- k$w_coverage +
       k$w_discovery * disc +
       k$w_far       * far +
       k$w_costanom  * costanom

  d$priority <- ave(p, d$region, FUN = function(z) z / sum(z)) * length(p)
  d$priority <- pmax(d$priority, 1e-6)
  d
}
df <- build_priority(df, cfg)

# =============================================================================
# 7.  GRTS DRAW — proportional-to-size allocation with a floor, since regions
#     range from 5 to ~14,600 candidate pixels. Flat n_per_region (as in the
#     synthetic demo) is not viable here — most small sites don't have enough
#     pixels to support it.
# =============================================================================
sframe <- sf::st_as_sf(df, coords = c('x', 'y'), crs = terra::crs(S_mean))

avail <- table(df$region)
lvls  <- names(avail)

n_regions <- length(lvls)
budget_for_prop <- max(cfg$n_total_base - cfg$min_per_region * n_regions, 0)
weight <- as.numeric(avail)^cfg$size_exponent
prop <- weight / sum(weight)

n_base <- cfg$min_per_region + round(prop * budget_for_prop)
n_base <- pmin(n_base, as.numeric(avail), cfg$max_per_region)  # availability + per-region cap
n_base <- setNames(n_base, lvls)

n_over <- pmin(round(n_base * cfg$oversample_mult), as.numeric(avail) - n_base)
n_over <- setNames(n_over, lvls)

alloc_tbl <- data.frame(
  region = lvls, Site = site_lkp[lvls], available = as.numeric(avail),
  n_base = n_base[lvls], n_over = n_over[lvls]
)
cat('\n--- Region allocation (proportional-to-size, floor = ', cfg$min_per_region, ') ---\n', sep = '')
print(alloc_tbl[order(-alloc_tbl$available), ], row.names = FALSE)
message(sprintf('\nRequested %d base + %d oversample across %d regions (target budget %d).',
                sum(n_base), sum(n_over), n_regions, cfg$n_total_base))

legacy_sites <- NULL   # e.g. st_read("historical_visited.gpkg") |> st_transform(crs(sframe))

dsgn <- spsurvey::grts(
  sframe       = sframe,
  n_base       = n_base,
  n_over       = n_over,
  stratum_var  = 'region',
  seltype      = 'proportional',
  aux_var      = 'priority',
  legacy_sites = legacy_sites,
  mindis       = cfg$mindis
)

base_sites <- dsgn$sites_base |> dplyr::mutate(set = 'base')
over_sites <- dsgn$sites_over |> dplyr::mutate(set = 'over')
message(sprintf('Drew %d base + %d oversample sites.', nrow(base_sites), nrow(over_sites)))

# =============================================================================
# 8.  cLHS COMPARISON — cost-constrained conditioned LHS on the same axes.
# =============================================================================
### NOTE: access_cost here is still the cost-anomaly proxy (Rq), same
### placeholder caveat as the synthetic script — swap in a real
### trailhead/road effort surface when available.
df$access_cost <- df$Rq
clhs_cols <- c('S', 'E', 'R')
total_n   <- sum(n_base)

clhs_res <- clhs(
  x        = as.data.frame(df[, c(clhs_cols, 'access_cost')]),
  size     = total_n,
  cost     = 'access_cost',
  iter     = 5000,
  simple   = FALSE,
  progress = FALSE
)
clhs_pts <- df[clhs_res$index_samples, ]

# =============================================================================
# 9.  Joint-coverage comparison: GRTS vs cLHS on the E-R plane + S marginal
# =============================================================================
# grts() carries S/E/Cost/Sq/Eq/Rq/priority straight through from sframe's
# attribute columns, so the drawn sites already have them - no re-join needed.
base_keyed <- base_sites |> sf::st_drop_geometry()

cover_plot <- dplyr::bind_rows(
  base_keyed |> dplyr::transmute(Eq, Rq, Sq, method = 'GRTS'),
  clhs_pts   |> dplyr::transmute(Eq, Rq, Sq, method = 'cLHS')
)

p_plane <- ggplot(cover_plot, aes(Eq, Rq)) +
  geom_point(data = df, color = 'grey88', size = 0.3) +
  geom_point(aes(color = method), size = 1.6) +
  geom_hline(yintercept = cfg$r_cap_q, linetype = 3) +
  facet_wrap(~method) +
  labs(title = sprintf('Realized joint coverage on the E-R plane (%s, real data)', res_tag),
       subtitle = 'grey = available pixels; dotted = R-pocket cap',
       x = 'Euclidean (within-region quantile)',
       y = 'Cost residual R (within-region quantile)') +
  theme_minimal() + theme(legend.position = 'none')

p_marg <- ggplot(cover_plot, aes(Sq, fill = method)) +
  geom_histogram(bins = 12, alpha = .6, position = 'identity') +
  labs(title = 'Suitability marginal coverage (identifiability axis)',
       x = 'Suitability (within-region quantile)', y = 'sites') +
  theme_minimal()

dir.create(file.path(p2res, 'GroundTruthSampling'), showWarnings = FALSE)
ggsave(file.path(p2res, 'GroundTruthSampling', paste0(cfg$ver, '-Test-coverage_E_R_plane.png')),
       p_plane, width = 9, height = 4.2, dpi = 130)
ggsave(file.path(p2res, 'GroundTruthSampling', paste0(cfg$ver, '-Test-coverage_S_marginal.png')),
       p_marg, width = 7, height = 3.6, dpi = 130)

# =============================================================================
# OUTPUTS for the field
# =============================================================================
out_base <- file.path(p2res, 'GroundTruthSampling', paste0(cfg$ver, '-Test-base.gpkg'))
out_over <- file.path(p2res, 'GroundTruthSampling', paste0(cfg$ver, '-Test-oversample.gpkg'))
### The ':' characters in cfg$ver (e.g. "PA1:3DO:0") make GDAL's driver-guessing
### misread the path as a "driver:connection" string - name the driver explicitly.
sf::st_write(base_sites, out_base, delete_dsn = TRUE, quiet = TRUE, driver = 'GPKG')
sf::st_write(over_sites, out_over, delete_dsn = TRUE, quiet = TRUE, driver = 'GPKG')

message('\nDone. Inspect: alloc_tbl, ', out_base, ', ', out_over,
        '\nand the coverage_E_R_plane / coverage_S_marginal PNGs in results/GroundTruthSampling/')
