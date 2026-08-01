## Shared ggplot2 theme and scales for SDManuscript.Rmd figures.
## source('theme_manuscript.R') from within the Rmd's working directory.

theme_manuscript <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = 'grey90', linewidth = 0.2),
      axis.text         = ggplot2::element_text(size = ggplot2::rel(0.75)),
      axis.title        = ggplot2::element_text(size = ggplot2::rel(0.85)),
      strip.text        = ggplot2::element_text(size = ggplot2::rel(0.9), face = 'bold'),
      strip.background  = ggplot2::element_blank(),
      legend.position   = 'bottom',
      legend.title      = ggplot2::element_text(size = ggplot2::rel(0.8)),
      legend.text       = ggplot2::element_text(size = ggplot2::rel(0.7)),
      plot.title        = ggplot2::element_text(size = ggplot2::rel(1), face = 'bold', hjust = 0),
      plot.margin       = ggplot2::margin(4, 4, 4, 4)
    )
}

## Variant for map panels in multi-panel composites: axis text/ticks dropped
## (lat/lon labels just add clutter once you're tiling 3-4 maps together),
## a light border kept so each panel still reads as a distinct map "window".
theme_manuscript_map <- function(base_size = 11) {
  theme_manuscript(base_size) +
    ggplot2::theme(
      axis.text   = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      axis.title  = ggplot2::element_blank(),
      panel.grid  = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = 'grey40', fill = NA, linewidth = 0.3)
    )
}

## Shared fill scale for habitat suitability (0-1 probability) rasters, so
## every suitability panel across every figure uses the same palette/limits.
scale_fill_suitability <- function(...) {
  ggplot2::scale_fill_viridis_c(name = 'Predicted\nsuitability', limits = c(0, 1), ...)
}

## Shared fill scale for plant density/count rasters.
scale_fill_density <- function(...) {
  ggplot2::scale_fill_viridis_c(name = 'Predicted\ndensity', option = 'magma', ...)
}
