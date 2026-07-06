# =============================================================================
# Stratified ground-truth design for SDM verification
#
# Stratify across:  S  = predicted suitability (ensemble mean)
#                   E  = Euclidean distance from known occurrences
#                   R  = WITHIN-REGION residual of cost-distance on E
#                        ("resistance anomaly" - cost beyond what proximity buys)
#                  [D] = WITHIN-suitability residual of ensemble SD on the mean
#                        (kept only if the diagnostic shows it carries signal)
#
# Draw:   GRTS, stratified by region, unequal probability ~ a priority surface
#         you build and can defend, with a 100% oversample (n_over = n_base)
#         and reverse-hierarchical ordering for inaccessible-draw replacement.
#
# Compare against: cost-constrained conditioned Latin hypercube (clhs) on the
#         same stack, so you can see realized joint coverage side by side.
#
# This file RUNS END TO END on a synthetic landscape shaped like yours, so you
# can watch the pipeline behave before swapping in real rasters. Every place you
# substitute real data is marked   ### SWAP-IN ###.
# =============================================================================

suppressPackageStartupMessages({
  library(terra)      # rasters
  library(sf)         # vector / sample frame
  library(dplyr)      # wrangling + group-wise residualization
  library(mgcv)       # gam smooths for the two residualizations
  library(spsurvey)   # grts()  (>= 5.0 interface)
  library(clhs)       # cost-constrained conditioned LHS (comparison)
  library(ggplot2)    # coverage diagnostics
})

set.seed(42)

# -----------------------------------------------------------------------------
# 0.  CONFIG - every design choice is a knob here, nothing buried downstream
# -----------------------------------------------------------------------------
cfg <- list(
  n_per_region   = 30,     # base sites per region (equal allocation)
  oversample_mult = 1.0,   # n_over = mult * n_base ; 

  # priority-surface weights (relative; rescaled within region). These ARE the
  # design - set them deliberately, report them, defend them.
  w_coverage   = 1.0,      # flat base: guarantees span of S for identifiability
  w_discovery  = 1.0,      # bump in the contested mid-S band (find new occ.)
  w_far        = 0.8,      # reward high-E tail (break the "near" truncation)
  w_costanom   = 1.0,      # reward |R| (where the cost surface earns its keep)
  w_disagree   = 0.0,      # set > 0 only if the D diagnostic says D survives

  discovery_center = 0.5,  # mid suitability, on the within-region quantile scale
  discovery_sigma  = 0.18,

  # R-POCKET CAP - the mitigation for the R/inaccessibility threat.
  # Stop upweighting cells above this within-region R quantile, so the
  # oversample doesn't burn itself re-drawing unreachable resistance pockets.
  r_cap_q = 0.90,

  # decision threshold for keeping the D (suitability-disagreement) axis:
  # fraction of SD variance left UNexplained by a smooth of the mean.
  d_keep_resid_frac = 0.15
)

# =============================================================================
# 1.  INPUTS   synthetic generator here; replace the four rasters + occurrences
# =============================================================================
### SWAP-IN ###  Replace this whole block with:
###   S_mean <- rast("suitability_mean.tif")     # ensemble mean, [0,1]
###   S_sd   <- rast("suitability_sd.tif")        # ensemble SD
###   E      <- rast("euclidean_dist.tif")        # from known occurrences
###   Cost   <- rast("cost_dist.tif")             # GLOBAL least-cost from occ.
###   region <- rast("region_id.tif")             # ~14 mask regions, integer/factor
###   occ    <- st_read("known_occurrences.gpkg") # points, for resolution check
### Make sure all five rasters share extent/res/CRS (resample/align if not).

make_synthetic <- function(n = 180, n_occ = 25, n_tiles = 4) {
  tmpl <- rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n,
               crs = "EPSG:32612")

  # --- suitability: a few smooth gaussian ridges + noise, squashed to (0,1)
  xy <- xyFromCell(tmpl, 1:ncell(tmpl))
  bump <- function(cx, cy, s) exp(-((xy[,1]-cx)^2 + (xy[,2]-cy)^2) / (2*s^2))
  z <- 1.6*bump(55,120,28) + 1.3*bump(130,60,24) + 1.1*bump(95,95,40)
  z <- z + rnorm(length(z), 0, 0.15)
  S_mean <- tmpl; values(S_mean) <- plogis(scale(z)[,1]*1.6)

  # --- ensemble SD: deliberately MOSTLY a parabola of the mean (peaks ~0.45),
  #     i.e. your "variance all centered at 0.4-0.5", + a little independent
  #     structure so the D diagnostic has something real to rule in or out.
  m <- values(S_mean)[,1]
  sd_det <- 0.28 * (1 - (2*(m - 0.45))^2); sd_det[sd_det < 0] <- 0
  sd_extra <- 0.05 * bump(130,130,30)                     # a genuine hot-spot
  S_sd <- tmpl; values(S_sd) <- pmax(0.01, sd_det + sd_extra + rnorm(length(m),0,0.02))

  # --- regions: coarse tiling -> ~n_tiles^2 regions (drop a couple via mask)
  reg <- tmpl
  tx <- cut(xy[,1], n_tiles, labels = FALSE)
  ty <- cut(xy[,2], n_tiles, labels = FALSE)
  values(reg) <- (ty - 1)*n_tiles + tx
  # complex-terrain mask: knock out a corner + a band so it's irregular & < n^2
  keep <- !( (xy[,1] > 150 & xy[,2] < 35) | (abs(xy[,1]-xy[,2]) < 6) )
  reg[!keep] <- NA
  region <- reg

  # --- known occurrences: high-suitability cells
  pcell <- sample(which(keep), n_occ, prob = m[keep]^3)
  occ <- vect(xyFromCell(tmpl, pcell), type = "points", crs = crs(tmpl))

  # --- Euclidean distance from occurrences
  E <- distance(tmpl, occ); E <- mask(E, region)

  # --- COST distance (synthetic stand-in for your real least-cost surface):
  #     resistance rises where suitability is low + ridge noise; cost tracks E
  #     but is inflated through resistant ground, so R will carry real structure.
  res_local <- 1 + 1.8*(1 - S_mean) + 0.6*bump(40,40,30)
  Cost <- E * res_local * exp(rnorm(ncell(tmpl), 0, 0.25))
  Cost <- mask(Cost, region)

  list(S_mean = S_mean, S_sd = S_sd, E = E, Cost = Cost,
       region = region, occ = occ, template = tmpl)
}

inp <- make_synthetic()

# =============================================================================
# 2.  Assemble the per-pixel table over the mask
# =============================================================================
stk <- c(inp$S_mean, inp$S_sd, inp$E, inp$Cost, inp$region)
names(stk) <- c("S", "SD", "E", "Cost", "region")

df <- as.data.frame(stk, cells = TRUE, xy = TRUE, na.rm = TRUE) %>%
  mutate(region = factor(region)) %>%
  filter(!is.na(region))

message(sprintf("Pixels in frame: %s across %d regions",
                format(nrow(df), big.mark = ","), nlevels(df$region)))

# =============================================================================
# 3.  DIAGNOSTIC FIRST - does the suitability-disagreement axis (D) survive?
#     Fit SD ~ s(mean). If the residual is nearly flat, D is redundant with S
#     and we drop it honestly. If it has spread, keep D and oversample its tail.
# =============================================================================
m_sd  <- gam(SD ~ s(S, k = 8), data = df, method = "REML")
df$D  <- residuals(m_sd)                       # contested-beyond-expectation
resid_frac <- var(df$D) / var(df$SD)

keep_D <- resid_frac >= cfg$d_keep_resid_frac
message(sprintf(
  "D diagnostic: %.1f%% of SD variance is independent of the mean -> %s D axis.",
  100*resid_frac, if (keep_D) "KEEP" else "DROP"))
if (keep_D && cfg$w_disagree == 0) cfg$w_disagree <- 0.7  # auto-enable if signal

# =============================================================================
# 4.  Region-aware cost orthogonalization:  R = within-region residual of C on E
#     Group-wise smooth = "this pixel is cheaper/costlier than its OWN region's
#     E-C relationship predicts." Orthogonal to a region intercept by design,
#     and it is the partial cost effect you actually want to defend.
#     (A single global s(E,by=region) would let R proxy region -> collinear with
#      the random intercept. Group-wise avoids that.)
# =============================================================================
df <- df %>%
  group_by(region) %>%
  group_modify(~{
    if (nrow(.x) >= 30 && length(unique(.x$E)) > 10) {
      mm <- tryCatch(gam(Cost ~ s(E, k = 6), data = .x, method = "REML"),
                     error = function(e) lm(Cost ~ poly(E, 2), data = .x))
      .x$R <- resid(mm)
    } else .x$R <- as.numeric(scale(.x$Cost - .x$E))   # tiny-region fallback
    .x
  }) %>%
  ungroup()

# NOTE (modelling-time hook, not used for design): the BETWEEN-region cost
# effect is gone from this within-region R by construction. If you want it back,
# reconstruct a Mundlak decomposition at fit time from region means of a GLOBAL
# C-on-E residual, and enter region-mean(R_global) + within-R as two terms.

# =============================================================================
# 5.  Within-region standardization (rank/quantile) so no massif dominates the
#     global strata. Everything below is built on these uniform [0,1] scales.
# =============================================================================
qn <- function(v) (rank(v, ties.method = "average") - 0.5) / length(v)
df <- df %>%
  group_by(region) %>%
  mutate(Sq = qn(S), Eq = qn(E), Rq = qn(R), Dq = qn(D)) %>%
  ungroup()

# =============================================================================
# 6.  PRIORITY SURFACE — the inclusion-probability blend (the knob-set).
#     priority ∝ w_coverage              (flat: keep full S span identifiable)
#             + w_discovery * mid-S bump  (contested frontier: find new occ.)
#             + w_far       * Eq          (break the historical near-truncation)
#             + w_costanom  * |Rq-0.5|*2  (where cost surface earns its keep) ...
#                ... but CAPPED above r_cap_q so we don't chase unreachable pockets
#             + w_disagree  * Dq          (model adjudicators, if D survived)
# =============================================================================
build_priority <- function(d, k) {
  disc <- exp(-((d$Sq - k$discovery_center)^2) / (2 * k$discovery_sigma^2))
  far  <- d$Eq
  # |R| signal, then flatten the extreme upper tail so the oversample isn't
  # spent re-drawing resistance pockets that keep coming back inaccessible:
  costanom <- abs(d$Rq - 0.5) * 2
  costanom[d$Rq > k$r_cap_q] <- abs(k$r_cap_q - 0.5) * 2
  dis  <- if (k$w_disagree > 0) d$Dq else 0

  p <- k$w_coverage +
       k$w_discovery * disc +
       k$w_far       * far +
       k$w_costanom  * costanom +
       k$w_disagree  * dis

  # rescale within region (grts aux_var draws prob ∝ value WITHIN stratum)
  d$priority <- ave(p, d$region, FUN = function(z) z / sum(z)) * length(z <- p)
  d$priority <- pmax(d$priority, 1e-6)         # strictly positive for aux_var
  d
}
df <- build_priority(df, cfg)

# =============================================================================
# 7.  GRTS DRAW — stratified by region, unequal prob ∝ priority, 100% oversample
# =============================================================================
sframe <- st_as_sf(df, coords = c("x", "y"), crs = crs(inp$template))

lvls    <- levels(df$region)
n_base  <- setNames(rep(cfg$n_per_region, length(lvls)), lvls)
n_over  <- setNames(round(cfg$n_per_region * cfg$oversample_mult), lvls)
n_over  <- setNames(rep(n_over[1], length(lvls)), lvls)

### SWAP-IN ###  legacy / already-visited sites to force into the design.
### sf POINTS, same CRS as sframe. grts includes them, accounts for their
### inclusion probabilities, and keeps the rest spatially balanced around them.
### Set to NULL if none. (Optionally add a legacy_var flag column + legacy_var=.)
legacy_sites <- NULL   # e.g. st_read("historical_visited.gpkg") |> st_transform(crs(sframe))

dsgn <- grts(
  sframe       = sframe,
  n_base       = n_base,
  n_over       = n_over,
  stratum_var  = "region",
  seltype      = "proportional",  # inclusion prob ∝ aux_var within stratum
  aux_var      = "priority",
  legacy_sites = legacy_sites,    # force in prior/known-visited points
  mindis       = 100              # min spacing between sites, in CRS UNITS
                                  # -> 100 only means 100 m if CRS is metric (yours is)
)

base_sites <- dsgn$sites_base %>% mutate(set = "base")
over_sites <- dsgn$sites_over %>% mutate(set = "over")  # already in RHO order
message(sprintf("Drew %d base + %d oversample sites.",
                nrow(base_sites), nrow(over_sites)))

# =============================================================================
# 8.  cLHS COMPARISON — cost-constrained conditioned LHS on the same axes.
#     cost = a travel-effort proxy (use real access cost when you have it).
# =============================================================================
### SWAP-IN ###  access_cost should be effort from ROAD/TRAILHEAD, not from
### occurrences — that distinction is exactly what decouples it from R.
df$access_cost <- df$Rq                 # placeholder proxy for the demo
clhs_cols <- c("S", "E", "R", if (keep_D) "D")
total_n   <- cfg$n_per_region * length(lvls)

clhs_res <- clhs(
  x      = as.data.frame(df[, c(clhs_cols, "access_cost")]),
  size   = total_n,
  cost   = "access_cost",
  iter   = 5000,
  simple = FALSE,
  progress = FALSE
)
clhs_pts <- df[clhs_res$index_samples, ]

# =============================================================================
# 9.  REALIZED COVERAGE + THE R-THREAT MISSINGNESS INSTRUMENT
# =============================================================================

## 9a. Simulate accessibility on the synthetic so you can SEE the threat fire.
### SWAP-IN ###  In the field this is just: did the crew reach the pixel? (1/0)
df$p_inaccess <- plogis(-1.2 + 2.5 * (df$Rq - 0.5))     # high-R -> often unreachable
sel <- df %>% inner_join(base_sites %>% st_drop_geometry() %>%
                           select(any_of("cell")), by = "cell")
# Map drawn sites back to the frame to pull Rq / inaccessibility:
key <- df %>% select(cell, region, Sq, Eq, Rq, p_inaccess)
base_keyed <- base_sites %>% st_drop_geometry() %>% left_join(key, by = "cell")
over_keyed <- over_sites %>% st_drop_geometry() %>% left_join(key, by = "cell")

base_keyed$inaccessible <- rbinom(nrow(base_keyed), 1, base_keyed$p_inaccess)

## 9b. Replace inaccessible base draws by walking down the RHO oversample list,
##     preserving spatial balance + the priority weighting (that's the point).
n_need <- sum(base_keyed$inaccessible)
replacements <- over_keyed %>%
  mutate(inaccessible = rbinom(n(), 1, p_inaccess)) %>%
  filter(inaccessible == 0) %>%        # reachable, in RHO order
  slice_head(n = n_need)

realized <- bind_rows(
  base_keyed %>% filter(inaccessible == 0),
  replacements
)
message(sprintf("Base inaccessible: %d/%d. Backfilled %d from oversample.",
                n_need, nrow(base_keyed), nrow(replacements)))

## 9c. THE diagnostic: does dropout pile up in the high-R stratum?
##     Flat across R-bins -> MCAR-ish, proceed. Pile-up in high R -> you have
##     DOCUMENTED informative missingness, and your estimand is honestly
##     "cost contribution among ACCESSIBLE habitat." Report this table.
miss_by_R <- base_keyed %>%
  mutate(Rbin = cut(Rq, c(0,.2,.4,.6,.8,1), include.lowest = TRUE)) %>%
  group_by(Rbin) %>%
  summarise(drawn = n(), inaccessible = sum(inaccessible),
            drop_rate = mean(inaccessible), .groups = "drop")
cat("\n--- Inaccessibility by within-region R quantile (the R-threat) ---\n")
print(miss_by_R)

cat("\n--- Realized R-coverage: drawn vs realized (after dropout) ---\n")
print(rbind(
  drawn    = summary(base_keyed$Rq),
  realized = summary(realized$Rq)
))

## 9d. Joint-coverage comparison: GRTS-realized vs cLHS on the E–R plane,
##     plus the S marginal (the identifiability axis).
cover_plot <- bind_rows(
  realized  %>% transmute(Eq, Rq, Sq, method = "GRTS (realized)"),
  clhs_pts  %>% transmute(Eq, Rq, Sq, method = "cLHS")
)
p_plane <- ggplot(cover_plot, aes(Eq, Rq)) +
  geom_point(data = df, color = "grey88", size = 0.3) +
  geom_point(aes(color = method), size = 1.6) +
  geom_hline(yintercept = cfg$r_cap_q, linetype = 3) +
  facet_wrap(~method) +
  labs(title = "Realized joint coverage on the E–R plane",
       subtitle = "grey = available pixels; dotted = R-pocket cap",
       x = "Euclidean (within-region quantile)",
       y = "Cost residual R (within-region quantile)") +
  theme_minimal() + theme(legend.position = "none")

p_marg <- ggplot(cover_plot, aes(Sq, fill = method)) +
  geom_histogram(bins = 12, alpha = .6, position = "identity") +
  labs(title = "Suitability marginal coverage (identifiability axis)",
       x = "Suitability (within-region quantile)", y = "sites") +
  theme_minimal()

ggsave("coverage_E_R_plane.png", p_plane, width = 9, height = 4.2, dpi = 130)
ggsave("coverage_S_marginal.png", p_marg,  width = 7, height = 3.6, dpi = 130)

# =============================================================================
# OUTPUTS for the field
# =============================================================================
### SWAP-IN ###  write realized + oversample (in RHO order) to GPKG for GPS:
# st_write(st_as_sf(realized,  ...), "field_sites_base.gpkg")
# st_write(over_sites,                "field_sites_oversample_RHO_order.gpkg")

message("\nDone. Inspect: miss_by_R, coverage_E_R_plane.png, coverage_S_marginal.png")
message("If miss_by_R shows high-R pile-up, your estimand is cost-among-accessible — state it.")
