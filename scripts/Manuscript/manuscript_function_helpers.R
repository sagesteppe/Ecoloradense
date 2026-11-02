## Small reusable helpers for the manuscript's inline statistics (sourced after
## theme_manuscript.R). Keep these generic - model-specific wrangling stays in the chunk.

## Format a p-value as markdown scientific notation, e.g. 2.6x10^-8^ (renders as 2.6x10^-8).
## Tiny p-values go scientific; p >= `plain_above` print as a plain decimal (e.g. 0.30).
fmt_p <- function(p, plain_above = 0.01, digits = 1) {
  if (is.na(p)) return('NA')
  if (p >= plain_above) return(formatC(p, format = 'f', digits = 2))
  e <- floor(log10(p))
  sprintf(paste0('%.', digits, 'fx10^%d^'), p / 10^e, e)
}

## Partial R^2 of one term in a fitted lm: the share of residual variance it uniquely explains,
## = SS(term dropped) / (RSS_full + SS(term dropped)), from a single-term F drop.
pr2_term <- function(model, term) {
  d  <- drop1(model, scope = term, test = 'F')
  ss <- d[term, 'Sum of Sq']
  ss / (stats::deviance(model) + ss)
}

## F-test p-value for dropping one term from an lm (whole factor, any df).
pval_term <- function(model, term) drop1(model, scope = term, test = 'F')[term, 'Pr(>F)']


cm_lab <- function(m) {
  m <- str_replace(m, 'Constant: (High|Low)-phase training mean.*',   'Arithmetic mean')
  m <- str_replace(m, 'Constant: (High|Low)-phase training median.*', 'Median')
  m <- str_replace(m, 'GLM Poisson \\((High|Low)-phase\\): suitability only', 'GLM Poisson ~ suitability')
  m <- str_replace(m, 'Spat[.]', '(Spatial CV)')
  m <- str_replace(m, 'Ordinary kriging.*', 'Ordinary Kriging')          # drop '(raw count)'
  m <- str_replace(m, 'IDW \\(gstat default, idp=2\\)', 'IDW (inverse-square)')
  m <- str_replace(m, 'brms Hurdlenegbinomial', 'Hurdle NegBinomial')
  m <- str_replace(m, 'brms Negbinomial', 'NegBinomial')
  m
}

cm_groups <- tibble::tribble(
  ~Group,                         ~Model,
  'Null models',                  'Arithmetic mean',
  'Null models',                  'Median',
  'Interpolation & GLM baseline', 'Ordinary Kriging',
  'Interpolation & GLM baseline', 'IDW (inverse-square)',
  'Interpolation & GLM baseline', 'kNN (within pop)',
  'Interpolation & GLM baseline', 'kNN (across pop)',
  'Interpolation & GLM baseline', 'GLM Poisson ~ suitability',
  'Gradient-boosted trees',       'XGB Tweedie',
  'Gradient-boosted trees',       'XGB Tweedie (Spatial CV)',
  'Gradient-boosted trees',       'XGB Poisson',
  'Gradient-boosted trees',       'XGB Poisson (Spatial CV)',
  'Gradient-boosted trees',       'LGBM Poisson',
  'Gradient-boosted trees',       'LGBM Poisson (Spatial CV)',
  'Bayesian (brms)',              'Hurdle NegBinomial',
  'Bayesian (brms)',              'NegBinomial'
)

cm_fix <- function(x) ifelse(is.finite(x) & x <= 1e6, round(x, 2), NA_real_)  

## One cell block: phase (High/Low) x domain (interp = plots within the sampled
## range; extrap = the Antero site, 86 km out). Antero is low-elevation only, so
## there is no high-elevation extrapolation block.
cm_slice <- function(phase, domain) {
  keep <- if (domain == 'extrap') function(s) s == 'Antero' else function(s) s != 'Antero'
  pb <- cm_pred |> filter(ScoredBy == 'Phase B (per-phase)', Phase == phase, keep(Site))
  pa <- cm_pred |> filter(ScoredBy == 'Phase A (pooled)', Phase == phase, keep(Site),
                          Tier == '3 Spatial interpolation' | Model %in% brms_keep)
  bind_rows(pb, pa) |>
    filter(is.finite(AbsError)) |>   # per-phase baselines are NA on their off-phase plots; drop those
    group_by(Model = cm_lab(Model)) |>
    summarise(all = mean(AbsError), occ = mean(AbsError[Occupied]), .groups = 'drop')
}

nn  <- function(ph, ant) sum(cm_n$Phase == ph & (cm_n$Site == 'Antero') == ant)
nno <- function(ph, ant) sum(cm_n$Phase == ph & (cm_n$Site == 'Antero') == ant & cm_n$Occupied)


get_p <- function(term_lbl) { v <- occ_lrt$p_value[occ_lrt$term == term_lbl]; if (length(v) == 1) signif(v, 2) else 'XX' }

get_band <- function(b) { v <- bands$n_presences[bands$band == b]; if (length(v) == 1) v else 'XX' }


## Occurrence/absence sample sizes per resolution, read directly from the
## presence files each PA-ratio search was built on (absences here are only
## each file's own baked-in absences - NOT the shared iter1 absence pool
## that adaptive_PAratio_search() also draws from, since that pool
## ('iter1-pa.gpkg') isn't available in this data checkout). Falls back to
## 'XX' per-resolution when a presence file is missing (e.g. the 1/3
## arc-second presence-iter1 file doesn't exist yet).
sample_sizes <- function(res_key) {
  f <- file.path(PROJ_ROOT,'data', 'Data4modelling', paste0(res_key, 'm-presence-iter1.gpkg'))
  if (!file.exists(f)) return(c(occ = NA_integer_, abs = NA_integer_))
  d <- sf::st_read(f, quiet = TRUE)
  c(occ = sum(d$Presenc == 1, na.rm = TRUE), abs = sum(d$Presenc == 0, na.rm = TRUE))
}

## min/max of metric m over all replicates at each group, looked up for the `keys` of interest

rng <- function(df, group_cols, m, keys) {
  k <- do.call(paste, df[group_cols])
  list(lo = as.numeric(tapply(df[[m]], k, min, na.rm = TRUE)[keys]),
       hi = as.numeric(tapply(df[[m]], k, max, na.rm = TRUE)[keys]))
}