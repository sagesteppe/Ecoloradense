# =============================================================================
# Stratified ground-truth design for SDM verification — REAL DATA
# Re-scoped to the narrowed manuscript lens.
#
# WHAT CHANGED vs the previous version, and why:
#
#   * Focal model = 1-3arc PA1:2.7 (cfg$ver). This script runs against ONE
#     product and draws ONE set of field plots; the coarser products are
#     evaluated on those SAME plots (goal 2), so nothing here is per-resolution.
#
#   * DISTANCE BANDS drive PURPOSE. Euclidean distance E to the nearest known
#     occurrence is cut into three bands:
#         band 1   E <  band_near              -> MODEL EVALUATION   (goal 1)
#         band 2   band_near <= E < band_far   -> DISCOVERY, near    (goal 3)
#         band 3   E >= band_far               -> DISCOVERY, far     (goal 3)
#     Breaks are absolute metres, set to multiples of the 3-arc pixel (~90 m)
#     so they also fall on cell edges of the coarser products:
#         band_near = 270  (3 px)   below this we assume cells are not
#                                    dispersal-limited, so a predicted-suitable
#                                    cell with no plants is model error, not
#                                    just unreachable ground. Tightened from
#                                    450m (a bumblebee-scale dispersal
#                                    distance): the taxon's pollinators are
#                                    mostly flies and solitary bees, whose
#                                    foraging/dispersal range runs shorter.
#         band_far  = 540  (6 px)   beyond it a new find is a new (sub)
#                                    population rather than an extension.
#
#   * THREE SEPARATE DRAWS, not one blended sample — one shared engine
#     (run_band_draw) called once per band. Evaluation and discovery want
#     opposite things from S, so they get different priority surfaces:
#         evaluation (band 1): SPAN the suitability gradient, so the confusion
#                              matrix / ROC-PR has support at both ends.
#         discovery  (bands 2-3): chase HIGH-S cells, where undetected plants
#                              are most likely.
#     Kept as separate output layers (never silently mixed), mirroring the
#     core/exploratory split in transect_placement.R. Splitting the draw also
#     guarantees real sample size in each band — the old within-region Eq
#     quantile did not.
#
#   * Cost-distance / cost-anomaly R axis and the cLHS comparison are REMOVED.
#     Neither serves goals 1/3 under this lens (both were in-house diagnostics),
#     and |R| upweighting actually points the wrong way for discovery, where
#     you want REACHABLE (low cost-distance) cells. A hook is left below where
#     least-cost reachability could re-enter the discovery priority if wanted.
#
# NOTE (goal 5): plant density / census is read on these same plots, so this
# draw is also the density sampling frame — no separate density design.
# =============================================================================

suppressPackageStartupMessages({
  library(terra)      # rasters
  library(sf)         # vector / sample frame
  library(dplyr)      # wrangling
  library(spsurvey)   # grts()  (>= 5.0 interface)
  library(ggplot2)    # coverage diagnostics
})

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
set.seed(42)

# -----------------------------------------------------------------------------
# 0.  CONFIG
# -----------------------------------------------------------------------------
cfg <- list(
  ver = '1-3arc-Iteration1-PA1:2.7DO:0',   # focal model version

  # --- distance bands (metres, Euclidean to nearest known occurrence) ---------
  band_near = 270,   # < this: evaluation frame (assumed not dispersal-limited)
  band_far  = 540,   # >= this: discovery-far frame (new populations)

  # --- per-draw budgets (base sites, before oversample) -----------------------
  #     These are the field-effort dials. Defaults preserve the old ~360 total,
  #     weighted toward evaluation since that is the headline goal. RETUNE to
  #     the real crew budget.
  n_eval      = 180,   # band 1  (goal 1)
  n_disc_near = 120,   # band 2  (goal 3)
  n_disc_far  = 60,    # band 3  (goal 3)

  # --- region allocation (shared across all three draws) ----------------------
  min_per_region = 5,     # floor: every region present in a band gets >= this
                          # many base sites (capped at that region's available
                          # pixel count in that band).
  max_per_region = 60,    # ceiling: no region gets more than this. Excess from
                          # capped regions is NOT redistributed.
  size_exponent  = 0.5,   # allocate on avail^exponent, not raw avail. <1
                          # compresses the dominance of the few huge Sites and
                          # shifts share toward the mid tier; 1 = linear.
  oversample_mult = 0.5,  # n_over = mult * n_base, capped at remaining pixels

  # --- priority shape ---------------------------------------------------------
  eval_bins    = 10,      # evaluation: flatten S over this many within-band
                          # quantile bins (inverse-frequency weighting)
  disc_exponent = 1.5,    # discovery: priority = Sq_band ^ exponent (>1 sharpens
                          # the pull toward the top of the suitability gradient)

  mindis = 30             # min spacing between sites, CRS units (m; UTM 13N)
)

# =============================================================================
# 1.  REAL INPUTS  (suitability, AOA, region, occurrences -> E)
#     Cost-distance rasters are no longer loaded.
# =============================================================================
p2res  <- '../../results'
p2data <- '../../data/Data4modelling'

suit_path <- file.path(p2res, 'suitability_maps', paste0(cfg$ver, '-Pr.tif'))
aoa_path  <- file.path(p2res, 'suitability_maps', paste0(cfg$ver, '-AOA.tif'))

stopifnot(
  'Suitability raster not found' = file.exists(suit_path),
  'AOA raster not found'         = file.exists(aoa_path)
)

res_tag <- sub('-Iteration.*$', '', cfg$ver)   # e.g. "1-3arc" - tags outputs

S_mean <- terra::rast(suit_path)
aoa_r  <- terra::rast(aoa_path)

# --- region: named Site polygons -> one (multi)polygon per site -> rasterize --
region_v <- sf::st_read(
  file.path('..', '..', 'data', 'hikingTrails', 'GroundTruch_Areas.gpkg'),
  quiet = TRUE)

region_v$region_id <- as.integer(factor(region_v$Site))
region <- terra::rasterize(terra::vect(region_v), S_mean, field = 'region_id')

# --- known occurrences: source of Euclidean distance E ------------------------
### NOTE: uses the full-precision 3m presence record regardless of model
### resolution, same convention as elsewhere in the pipeline.
occ <- sf::st_read(file.path(p2data, '3m-presence-iter1.gpkg'), quiet = TRUE) |>
  dplyr::filter(Presenc == 1) |>
  sf::st_transform(terra::crs(S_mean))

E <- terra::distance(S_mean, terra::vect(occ))

# --- mask to AOA (inference domain) + region (feasible-to-groundtruth) --------
S_mean <- terra::mask(S_mean, aoa_r, maskvalues = 0) |> terra::mask(region)
E      <- terra::mask(E,      aoa_r, maskvalues = 0) |> terra::mask(region)

# =============================================================================
# 2.  Per-pixel table + BAND assignment
# =============================================================================
stk <- c(S_mean, E, region)
names(stk) <- c('S', 'E', 'region')

df <- as.data.frame(stk, cells = TRUE, xy = TRUE, na.rm = TRUE) |>
  dplyr::mutate(region = factor(region, labels = levels(factor(region_v$region_id)))) |>
  dplyr::filter(!is.na(region))

site_lkp <- setNames(region_v$Site, region_v$region_id)
df$Site  <- site_lkp[as.character(df$region)]

## band from E: purpose follows directly from distance-to-nearest-occurrence.
df$band <- cut(df$E,
               breaks = c(-Inf, cfg$band_near, cfg$band_far, Inf),
               labels = c('evaluation', 'discovery_near', 'discovery_far'),
               right  = FALSE)          # [0,near) [near,far) [far,Inf)

message(sprintf('Pixels in frame: %s across %d Sites',
                format(nrow(df), big.mark = ','), nlevels(df$region)))
print(table(band = df$band))

# =============================================================================
# 3.  PRIORITY SURFACES — one per purpose, both computed WITHIN band x region
#     so "span S" / "high S" mean relative to the cells actually available in
#     that band (near-occurrence cells skew suitable; spanning the full-frame S
#     range would just miss most of them).
# =============================================================================
qn <- function(v) (rank(v, ties.method = 'average') - 0.5) / length(v)

## evaluation: flatten the S marginal. Bin within-band S into eval_bins quantile
## bins per region and weight each cell by 1 / (cells in its bin) -> the draw is
## pulled toward whatever S values are locally rare, giving support across the
## gradient (incl. predicted-unsuitable cells) instead of piling up on the
## suitable cells that dominate near occurrences.
coverage_priority <- function(d, k) {
  d |>
    dplyr::group_by(region) |>
    dplyr::mutate(
      bin = cut(qn(S), breaks = seq(0, 1, length.out = k$eval_bins + 1),
                include.lowest = TRUE, labels = FALSE),
      priority = 1 / ave(S, bin, FUN = length)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-bin)
}

## discovery: chase high suitability. priority = within-band S quantile raised
## to disc_exponent, so the top of the gradient gets most of the effort.
## HOOK: to fold in dispersal reachability, multiply priority by a reachability
## term here (e.g. a rescaled inverse least-cost distance) before normalizing.
discovery_priority <- function(d, k) {
  d |>
    dplyr::group_by(region) |>
    dplyr::mutate(priority = qn(S) ^ k$disc_exponent) |>
    dplyr::ungroup()
}

# =============================================================================
# 4.  DRAW ENGINE — region-stratified GRTS with proportional-to-size allocation
#     (floor + ceiling), run once per band. Same allocation machinery as before,
#     just factored out so all three draws share it.
# =============================================================================
run_band_draw <- function(df_all, band_lab, n_total, priority_fn, purpose, k) {
  d <- dplyr::filter(df_all, band == band_lab)
  if (nrow(d) == 0) {
    message(sprintf('  [%s] no pixels in band "%s" - skipped.', purpose, band_lab))
    return(NULL)
  }
  d <- priority_fn(d, k)
  # normalize priority to mean 1 within region (grts aux_var wants relative wts)
  d$priority <- ave(d$priority, d$region, FUN = function(z) z / mean(z))
  d$priority <- pmax(d$priority, 1e-6)

  sframe <- sf::st_as_sf(d, coords = c('x', 'y'), crs = terra::crs(S_mean))

  avail <- table(droplevels(d$region))
  lvls  <- names(avail)
  n_reg <- length(lvls)

  budget_for_prop <- max(n_total - k$min_per_region * n_reg, 0)
  weight <- as.numeric(avail) ^ k$size_exponent
  prop   <- weight / sum(weight)

  n_base <- k$min_per_region + round(prop * budget_for_prop)
  n_base <- pmin(n_base, as.numeric(avail), k$max_per_region)
  n_base <- setNames(n_base, lvls)

  n_over <- pmin(round(n_base * k$oversample_mult), as.numeric(avail) - n_base)
  n_over <- setNames(pmax(n_over, 0), lvls)

  alloc <- data.frame(region = lvls, Site = site_lkp[lvls],
                      available = as.numeric(avail),
                      n_base = n_base[lvls], n_over = n_over[lvls])
  cat(sprintf('\n--- %s draw (band "%s") allocation ---\n', purpose, band_lab))
  print(alloc[order(-alloc$available), ], row.names = FALSE)
  message(sprintf('  [%s] %d base + %d oversample across %d Sites (target %d).',
                  purpose, sum(n_base), sum(n_over), n_reg, n_total))

  dsgn <- spsurvey::grts(
    sframe      = sframe,
    n_base      = n_base,
    n_over      = n_over,
    stratum_var = 'region',
    seltype     = 'proportional',
    aux_var     = 'priority',
    mindis      = k$mindis
  )

  dplyr::bind_rows(
    dsgn$sites_base |> dplyr::mutate(set = 'base'),
    dsgn$sites_over |> dplyr::mutate(set = 'over')
  ) |>
    dplyr::mutate(purpose = purpose, band = band_lab)
}

# =============================================================================
# 5.  THREE DRAWS
# =============================================================================
draw_eval <- run_band_draw(df, 'evaluation',    cfg$n_eval,
                           coverage_priority,  'evaluation',    cfg)
draw_near <- run_band_draw(df, 'discovery_near', cfg$n_disc_near,
                           discovery_priority, 'discovery_near', cfg)
draw_far  <- run_band_draw(df, 'discovery_far',  cfg$n_disc_far,
                           discovery_priority, 'discovery_far',  cfg)

sites_all <- dplyr::bind_rows(draw_eval, draw_near, draw_far)
message(sprintf('\nTotal drawn: %d sites (%d base, %d oversample).',
                nrow(sites_all),
                sum(sites_all$set == 'base'), sum(sites_all$set == 'over')))

# =============================================================================
# 6.  COVERAGE DIAGNOSTICS — GRTS only now. Two things worth eyeballing:
#     (a) evaluation sample spans S; (b) discovery sits high on S.
# =============================================================================
plot_df <- sites_all |> sf::st_drop_geometry() |> dplyr::filter(set == 'base')

p_S <- ggplot(plot_df, aes(S, fill = purpose)) +
  geom_histogram(bins = 15, alpha = .6, position = 'identity') +
  facet_wrap(~purpose, ncol = 1, scales = 'free_y') +
  labs(title = sprintf('Suitability coverage by purpose (%s)', res_tag),
       subtitle = 'evaluation should span S; discovery should sit high',
       x = 'Predicted suitability S', y = 'base sites') +
  theme_minimal() + theme(legend.position = 'none')

dir.create(file.path(p2res, 'GroundTruthSampling'), showWarnings = FALSE)
ggsave(file.path(p2res, 'GroundTruthSampling',
                 paste0(cfg$ver, '-Test-coverage_S_by_purpose.png')),
       p_S, width = 7, height = 6, dpi = 130)

# =============================================================================
# 7.  OUTPUTS for the field — one gpkg, a layer per purpose (never blended).
# =============================================================================
out_path <- file.path(p2res, 'GroundTruthSampling', paste0(cfg$ver, '-Test-plots.gpkg'))
### The ':' characters in cfg$ver make GDAL misread the path as a
### "driver:connection" string - name the driver explicitly.
first <- TRUE
for (pp in c('evaluation', 'discovery_near', 'discovery_far')) {
  layer_sf <- dplyr::filter(sites_all, purpose == pp)
  if (nrow(layer_sf) == 0) next
  sf::st_write(layer_sf, out_path, layer = pp, driver = 'GPKG',
               delete_dsn = first, append = !first, quiet = TRUE)
  first <- FALSE
}

message('\nDone. Inspect the coverage_S_by_purpose PNG and: ', out_path,
        '\nLayers: evaluation / discovery_near / discovery_far')
