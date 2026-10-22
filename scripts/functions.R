## Thin loader - every helper used to live here. Split into modules:
##   core_functions.R          identifiers, paths, small utilities
##   metric_functions.R        scoring (presenceScores, DCA, Cox calibration)
##   cost_distance_functions.R connectivity / least-cost pipeline
##   groundtruth_functions.R   field / transect / point-ordering helpers
##   cv_functions.R            splits, replicates, PAratio/distOrder search
##   suitability_functions.R   SDM (modeller, AOA, case weights, Platt)
##   density_functions.R       densityModeller, gbs, kNN null, tune helpers
## Sourcing this file loads all of them, so existing `source('functions.R')`
## callers keep working with no change.

PROJ_ROOT <- '/media/steppe/hdd/EriogonumColoradenseTaxonomy'

local({
  # Resolve this file's directory from the live source() stack - works whether
  # functions.R is sourced from scripts/ or from a subdir with `source('../functions.R')`.
  this_dir <- (function(){
    for (fr in rev(sys.frames())) {
      of <- fr$ofile
      if (!is.null(of)) return(dirname(normalizePath(of, mustWork = FALSE)))
    }
    getwd()
  })()
  for (mod in c('core_functions','metric_functions','cost_distance_functions',
                'groundtruth_functions','cv_functions','suitability_functions',
                'density_functions')) {
    source(file.path(this_dir, paste0(mod, '.R')))
  }
})
