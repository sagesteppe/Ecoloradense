## Raster-side entry points for Partitioner.R.
##
## Partitioner.R is logic-only (sf/FNN/twosamples). This file adds the Step 0 PCA on
## the production rastReader stack and two thin wrappers that compose Step 0 with
## build_partition():
##
##   build_partition_mapping(...)  -> inner CV = partition_folds       (spatial_mapping)
##   build_partition_lmpo(...)     -> inner CV = partition_folds_lopo  (LMPO; leave
##                                                                      multiple pops out)
##
## Both wrappers return the same shape (a build_partition() object); the only
## difference is the inner-CV strategy wired into `folds` and the index triple.
## Source Partitioner.R and functions.R (for rastReader/PROJ_ROOT) first.

if (!exists('build_partition'))
  stop('source scripts/final_sdm/Partitioner.R before partition_entrypoints.R')
if (!exists('rastReader'))
  stop('source scripts/functions.R before partition_entrypoints.R')

#' Step 0 - PCA on raster predictors
#'
#' Extracts the predictor stack at each presence point, drops rows with any NA
#' (same convention as modeller()), and runs centred/scaled PCA. Keeps the first
#' axes whose cumulative variance reaches `pca_var_target`. Returns the presences
#' as sf with an integer `.id` column and PC1..PCn score columns, plus the
#' retained PCA names.
#'
#' @param x_pres sf of presences (any CRS; env is extracted in rastReader's CRS)
#' @param rast_dat SpatRaster predictor stack (output of rastReader())
#' @param pca_var_target cumulative-variance target (default from partition_params())
partition_pca_scores <- function(x_pres, rast_dat, pca_var_target = 0.90) {
  stopifnot(inherits(x_pres, 'sf'), inherits(rast_dat, 'SpatRaster'))

  env <- terra::extract(rast_dat, x_pres, ID = FALSE)
  keep <- stats::complete.cases(env)
  if (!all(keep))
    message(sprintf('partition_pca_scores: dropped %d presences with NA predictors.',
                    sum(!keep)))
  x_pres <- x_pres[keep, , drop = FALSE]
  env    <- env[keep, , drop = FALSE]

  pca <- stats::prcomp(env, center = TRUE, scale. = TRUE)
  cumvar <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  n_keep <- max(1L, which(cumvar >= pca_var_target)[1])
  pcs <- paste0('PC', seq_len(n_keep))

  scores <- as.data.frame(pca$x[, seq_len(n_keep), drop = FALSE])
  names(scores) <- pcs
  x_pres$.id <- seq_len(nrow(x_pres))
  x_pres <- cbind(x_pres, scores)

  list(x = x_pres, pcs = pcs, pca = pca,
       cum_var = cumvar[n_keep], n_axes = n_keep)
}

#' spatial_mapping entry point - inner CV via partition_folds()
#'
#' @param x_pres sf of presences; CRS matches rast_dat (we transform if needed).
#' @param bg sf of background absences (Occurrence == 0); env isn't needed here.
#' @param rast_dat SpatRaster predictor stack (rastReader() output).
#' @param p partition_params() list; override fields as needed.
build_partition_mapping <- function(x_pres, bg, rast_dat, p = partition_params()) {
  if (sf::st_crs(x_pres) != sf::st_crs(rast_dat))
    x_pres <- sf::st_transform(x_pres, sf::st_crs(rast_dat))
  if (sf::st_crs(bg) != sf::st_crs(rast_dat))
    bg <- sf::st_transform(bg, sf::st_crs(rast_dat))
  s0 <- partition_pca_scores(x_pres, rast_dat, p$pca_var_target)
  part <- build_partition(s0$x, bg, s0$pcs, p, folds_fn = partition_folds)
  part$pca <- s0$pca; part$pcs <- s0$pcs
  part$cv_mode <- 'spatial_mapping'
  part
}

#' LMPO entry point - inner CV via partition_folds_lopo()
#'
#' Whole segment -> one fold, greedily load-balanced. "LMPO" because any single
#' fold may hold several connected components (populations), not just one.
#' Same input/output contract as build_partition_mapping().
build_partition_lmpo <- function(x_pres, bg, rast_dat, p = partition_params()) {
  if (sf::st_crs(x_pres) != sf::st_crs(rast_dat))
    x_pres <- sf::st_transform(x_pres, sf::st_crs(rast_dat))
  if (sf::st_crs(bg) != sf::st_crs(rast_dat))
    bg <- sf::st_transform(bg, sf::st_crs(rast_dat))
  s0 <- partition_pca_scores(x_pres, rast_dat, p$pca_var_target)
  part <- build_partition(s0$x, bg, s0$pcs, p, folds_fn = partition_folds_lopo)
  part$pca <- s0$pca; part$pcs <- s0$pcs
  part$cv_mode <- 'LMPO'
  part
}
