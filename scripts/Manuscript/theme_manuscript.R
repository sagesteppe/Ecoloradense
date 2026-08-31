## Shared ggplot2 theme and scales for SDManuscript.Rmd figures.
## source('theme_manuscript.R') from within the Rmd's working directory.

## Base-R-styled theme: white panel, full black box (bty = 'o'), black ticks
## and axis text, no gridlines - reads like base plot() rather than ggplot's
## usual grey-panel look.
theme_manuscript <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid        = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(colour = 'black', fill = NA, linewidth = 0.5),
      axis.ticks        = ggplot2::element_line(colour = 'black', linewidth = 0.4),
      axis.text         = ggplot2::element_text(colour = 'black', size = ggplot2::rel(0.75)),
      axis.title        = ggplot2::element_text(colour = 'black', size = ggplot2::rel(0.85)),
      strip.text        = ggplot2::element_text(size = ggplot2::rel(0.9), colour = 'black'),
      strip.background  = ggplot2::element_rect(fill = 'white', colour = 'black', linewidth = 0.4),
      legend.position    = 'bottom',
      legend.background = ggplot2::element_blank(),
      legend.key        = ggplot2::element_blank(),
      legend.title      = ggplot2::element_text(size = ggplot2::rel(0.8)),
      legend.text       = ggplot2::element_text(size = ggplot2::rel(0.7)),
      plot.title        = ggplot2::element_text(size = ggplot2::rel(1), hjust = 0.5),
      plot.margin       = ggplot2::margin(4, 4, 4, 4)
    )
}

## Variant for map panels in multi-panel composites: axis title dropped (no
## "Easting"/"Northing" needed), but axis ticks/text kept - mapPanel() only
## ever sets two breaks per axis (the panel's corners), so this reads as
## normal ggplot coordinate ticks rather than a full graticule. The two
## labels are justified inward (hjust/vjust recycle across breaks) so they
## sit right of the lower-left tick / left of the lower-right tick on x,
## and above the lower tick / below the upper tick on y, instead of
## centering on the tick and overhanging the edge of the figure. Extra top
## margin makes room for a flush top-left patchwork tag (theme_manuscript_
## tags()) so it doesn't collide with the upper-left latitude label.
theme_manuscript_map <- function(base_size = 11) {
  theme_manuscript(base_size) +
    ggplot2::theme(
      axis.title  = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(colour = 'black', size = ggplot2::rel(0.75), hjust = c(0, 1)),
      axis.text.y = ggplot2::element_text(colour = 'black', size = ggplot2::rel(0.75), vjust = c(0, 1)),
      plot.margin = ggplot2::margin(12, 0, 0, 0)
    )
}

## Enlarged legend swatches for composite map figures - the default
## legend.key.size (~1.2 lines) reads too small once a fill gradient is
## squeezed under a wide multi-panel composite. `width_mult`/`height_mult`
## scale relative to that default. Apply with `& theme_manuscript_legend()`
## alongside theme_manuscript_tags() after plot_annotation().
theme_manuscript_legend <- function(width_cm = 1.5, height_cm = 0.5) {
  ggplot2::theme(
    legend.key.width  = ggplot2::unit(width_cm, 'cm'),
    legend.key.height = ggplot2::unit(height_cm, 'cm')
  )
}

## Patchwork panel-tag styling ("A", "B", ... from plot_annotation(tag_levels
## = 'A')): flush to the top-left corner of each panel instead of
## patchwork's default inset margin. plot.tag.location = 'margin' draws it
## in the panel's margin gutter (the extra top plot.margin added by
## theme_manuscript_map()) rather than inside/over the panel, so it doesn't
## collide with the upper-left axis label. Apply with
## `& theme_manuscript_tags()` after plot_annotation().
theme_manuscript_tags <- function() {
  ggplot2::theme(
    plot.tag.position = 'topleft',
    plot.tag.location = 'margin',
    plot.tag = ggplot2::element_text(hjust = 0, vjust = 1, margin = ggplot2::margin(0, 0, 0, 0), face = 'bold')
  )
}

## Shared fill scale for habitat suitability (0-1 probability) rasters, so
## every suitability panel across every figure uses the same palette/limits.
## Diverging red -> yellow -> green through `midpoint`, the suitable/
## unsuitable cutoff - not necessarily the middle of the 0-1 range.
scale_fill_suitability <- function(midpoint = 0.45, ...) {
  ggplot2::scale_fill_gradient2(
    name = 'Predicted\nsuitability', low = '#D82216', mid = '#ffe68f', high = '#197A56',
    midpoint = midpoint, limits = c(0, 1), ...
  )
}

## Shared fill scale for plant density/count rasters. Pass a shared `limits`
## (e.g. range() across every site's cropped data) so multiple density
## panels stay on the same colour mapping and can share one legend.
scale_fill_density <- function(limits = NULL, ...) {
  ggplot2::scale_fill_gradient(
    name = 'Predicted\ndensity', high = '#0a477d', low = '#76B9F4',
    limits = limits, ...
  )
}

## Shared colour scale for figures comparing ground-verification phases:
## H0 (fit on pre-existing occurrence records only) vs H1 (refit after
## incorporating a round of ground-truthing). Reuses the same red/green
## semantic as scale_fill_suitability() (red = less trustworthy baseline,
## green = the ground-truth-informed result) for a consistent visual
## language across figures.
phase_colors <- c(H0 = '#D82216', H1 = '#197A56')
scale_color_phase <- function(...) ggplot2::scale_color_manual(values = phase_colors, ...)

## Shared colour/shape scales for figures comparing spatial vs non-spatial
## ("classic") train-test splits or cross-validation - reused across the
## PA-ratio search (H2) and, later, the density CV comparison (H4a).
split_colors <- c(Spatial = '#197A56', Classic = '#0A477D')
split_shapes <- c(Spatial = 16, Classic = 17)

scale_color_split <- function(...) ggplot2::scale_color_manual(values = split_colors, ...)
scale_fill_split  <- function(...) ggplot2::scale_fill_manual(values = split_colors, ...)
scale_shape_split <- function(...) ggplot2::scale_shape_manual(values = split_shapes, ...)

## crop (optionally) + tidy a raster for geom_raster.
rasterToDF <- function(r, bb = NULL) {
  if (!is.null(bb)) r <- terra::crop(r, terra::ext(bb[c('xmin', 'xmax', 'ymin', 'ymax')]))
  setNames(as.data.frame(r, xy = TRUE), c('x', 'y', 'value'))
}

## Null out predictions outside the model's area of applicability (AOA == 0)
## before plotting - those pixels are extrapolations, not real predictions.
mask_aoa <- function(r, aoa) {
  terra::mask(r, aoa, maskvalues = 0)
}

## Regular grid of points clipped to the outside-AOA area (AOA == 0), for
## drawing that area as a speckled/stippled fill instead of a blank NA gap.
## No ggpattern dependency: sample a grid over the region's bbox, then keep
## only points that actually fall inside the (possibly multi-part,
## irregularly shaped) outside-AOA polygon.
aoa_speckle <- function(aoa, spacing = NULL) {
  outside <- terra::as.polygons(aoa == 0, dissolve = TRUE)
  outside <- sf::st_as_sf(outside)
  outside <- outside[outside[[1]] == 1, ]
  if (nrow(outside) == 0) return(NULL)
  outside <- sf::st_union(outside)

  if (is.null(spacing)) spacing <- terra::res(aoa)[1] * 10
  grid_pts <- sf::st_make_grid(outside, cellsize = spacing, what = 'centers')
  pts <- sf::st_coordinates(grid_pts[sf::st_intersects(grid_pts, outside, sparse = FALSE)[, 1]])
  setNames(as.data.frame(pts), c('x', 'y'))
}

## Speckle points as a plot layer - minimal point size, drawn over whatever
## sits underneath (the raster's NA gap reads as blank/white there).
speckle_layer <- function(speckle, size = 0.1, color = 'grey30', alpha = 0.5) {
  ggplot2::geom_point(
    data = speckle, ggplot2::aes(x = x, y = y), inherit.aes = FALSE,
    size = size, shape = 16, color = color, alpha = alpha
  )
}

## Lon/lat tick breaks + labels for a projected map panel's x/y axes - just
## the panel's min/max on each axis (so ticks land at the lower-left,
## upper-left, and lower-right corners), labelled in degrees regardless of
## the raster's working CRS (UTM 13N here). Ordinary ggplot axis machinery
## does the rest (ticks/text drawn outside the panel, as normal).
lonlat_breaks <- function(x_range, y_range, crs) {
  x_pts <- sf::st_as_sf(data.frame(x = x_range, y = y_range[1]), coords = c('x', 'y'), crs = crs)
  y_pts <- sf::st_as_sf(data.frame(x = x_range[1], y = y_range), coords = c('x', 'y'), crs = crs)
  x_geo <- sf::st_coordinates(sf::st_transform(x_pts, 4326))
  y_geo <- sf::st_coordinates(sf::st_transform(y_pts, 4326))
  list(
    x_breaks = x_range,
    x_labels = sprintf('%.2f°%s', abs(x_geo[, 1]), ifelse(x_geo[, 1] >= 0, 'E', 'W')),
    y_breaks = y_range,
    y_labels = sprintf('%.2f°%s', abs(y_geo[, 2]), ifelse(y_geo[, 2] >= 0, 'N', 'S'))
  )
}

## Region bounding boxes drawn as outline rectangles - e.g. to show where
## the focal-region crops sit within a full-domain panel. If `tags` is
## given (named to match `bboxes`, e.g. c(cb = 'B', coche = 'C', antero =
## 'D') - the patchwork letter of each box's own inset panel), a grey-
## background label is stamped at each box's top-left corner so the
## full-domain panel cross-references its insets without a separate key.
bbox_layer <- function(bboxes, tags = NULL, color = 'black', linewidth = 0.6, linetype = '3313',
                        tag_fill = 'grey85', tag_size = 3) {
  bbox_df <- do.call(rbind, lapply(names(bboxes), function(nm) {
    b <- bboxes[[nm]]
    data.frame(name = nm, xmin = b[['xmin']], xmax = b[['xmax']], ymin = b[['ymin']], ymax = b[['ymax']])
  }))
  layers <- list(ggplot2::geom_rect(
    data = bbox_df, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, fill = NA, color = color, linewidth = linewidth, linetype = linetype
  ))
  if (!is.null(tags)) {
    bbox_df$tag <- tags[bbox_df$name]
    layers <- c(layers, list(ggplot2::geom_label(
      data = bbox_df, ggplot2::aes(x = xmin, y = ymax, label = tag),
      inherit.aes = FALSE, hjust = 0, vjust = 1, size = tag_size, fontface = 'bold',
      label.size = 0, label.padding = ggplot2::unit(0.1, 'lines'), fill = tag_fill
    )))
  }
  layers
}

## Centered "zoom" bounding box: shrinks bb to `factor` of its width/height
## around the same center - a nested inset box for a further zoomed-in
## panel (e.g. a close-up crop within one of the focal regions).
shrink_bbox <- function(bb, factor = 0.35) {
  cx <- mean(c(bb[['xmin']], bb[['xmax']]))
  cy <- mean(c(bb[['ymin']], bb[['ymax']]))
  hw <- (bb[['xmax']] - bb[['xmin']]) * factor / 2
  hh <- (bb[['ymax']] - bb[['ymin']]) * factor / 2
  c(xmin = cx - hw, xmax = cx + hw, ymin = cy - hh, ymax = cy + hh)
}

## Shared chrome for a map panel: in-panel upper-left title, lon/lat axis
## ticks, a scale bar, and (optionally) a north arrow.
decorate_map <- function(p, x_range, y_range, label = NULL, crs = NULL,
                          scale_bar = TRUE, north_arrow = FALSE, scale_width_hint = 0.25) {
  if (!is.null(label)) {
    p <- p + ggplot2::annotate(
      'label', x = x_range[1], y = y_range[2], label = label,
      hjust = 0, vjust = 1, size = 3.2, fontface = 'bold',
      linewidth = 0, fill = ggplot2::alpha('white', 0.75)
    )
  }
  if (!is.null(crs)) {
    lb <- lonlat_breaks(x_range, y_range, crs)
    p <- p +
      ggplot2::scale_x_continuous(breaks = lb$x_breaks, labels = lb$x_labels) +
      ggplot2::scale_y_continuous(breaks = lb$y_breaks, labels = lb$y_labels)
  }
  if (scale_bar) {
    ## width_hint controls the scale bar's target length as a fraction of
    ## panel width, which in turn decides ggspatial's m/km (or ft/mi)
    ## bucket: it picks km only once that target length exceeds 1600m
    ## (scalebar_params()'s internal threshold). Tightly zoomed panels can
    ## land just under that with the default 0.25 hint and render as an
    ## awkward '1000 m' instead of '1 km' - bump scale_width_hint on those
    ## panels to push the target comfortably past the boundary.
    p <- p + ggspatial::annotation_scale(
      location = 'tr', style = 'ticks', text_cex = 0.6,
      pad_y = grid::unit(if (north_arrow) 1.1 else 0.3, 'cm'),
      width_hint = scale_width_hint
    )
  }
  if (north_arrow) {
    p <- p + ggspatial::annotation_north_arrow(
      location = 'tr', which_north = 'grid',
      height = grid::unit(0.8, 'cm'), width = grid::unit(0.8, 'cm'),
      style = ggspatial::north_arrow_minimal()
    )
  }
  p
}

## One styled map panel, built from a tidied raster data frame. `crs` (the
## raster's CRS) is only needed to compute lon/lat axis breaks; pass
## crs = NULL to skip them (e.g. for the small bivariate legend). `bboxes`
## draws region bounding boxes as outline rectangles (e.g. on a full-domain
## panel, to show where the focal-region crops sit). `speckle` (from
## aoa_speckle()) stipples the outside-AOA gap instead of leaving it blank.
mapPanel <- function(df, scale_fn, label = NULL, crs = NULL, bboxes = NULL, bbox_tags = NULL,
                      speckle = NULL, scale_bar = TRUE, north_arrow = FALSE, scale_width_hint = 0.25) {
  x_range <- range(df$x, na.rm = TRUE)
  y_range <- range(df$y, na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_raster() +
    scale_fn() +
    ggplot2::coord_equal(xlim = x_range, ylim = y_range, expand = FALSE) +
    theme_manuscript_map()

  if (!is.null(speckle)) p <- p + speckle_layer(speckle)
  if (!is.null(bboxes)) p <- p + bbox_layer(bboxes, tags = bbox_tags)

  decorate_map(p, x_range, y_range, label, crs, scale_bar, north_arrow, scale_width_hint)
}

## Bivariate panel: suitability as fill (red = low, green = high) with
## prediction SE mapped to alpha, so more uncertain pixels fade toward the
## white background instead of needing a second, harder-to-read legend grid.
mapPanelBivariate <- function(df, label = NULL, crs = NULL, bboxes = NULL, speckle = NULL,
                               scale_bar = TRUE, north_arrow = FALSE) {
  x_range <- range(df$x, na.rm = TRUE)
  y_range <- range(df$y, na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = value, alpha = se)) +
    ggplot2::geom_raster() +
    scale_fill_suitability() +
    ggplot2::scale_alpha_continuous(name = 'Prediction\nSE', range = c(1, 0.15)) +
    ggplot2::coord_equal(xlim = x_range, ylim = y_range, expand = FALSE) +
    theme_manuscript_map()

  if (!is.null(speckle)) p <- p + speckle_layer(speckle)
  if (!is.null(bboxes)) p <- p + bbox_layer(bboxes)

  decorate_map(p, x_range, y_range, label, crs, scale_bar, north_arrow)
}

## Smoothed random-noise raster over a bounding box, rescaled to `range`.
## Stand-in for the final suitability/density rasters (still being
## produced) so map-panel styling can be previewed against something
## spatially plausible instead of static or an empty plot.
example_surface <- function(bb, res = 0.003, range = c(0, 1), smooth_iter = 8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  r <- terra::rast(
    xmin = bb[['xmin']], xmax = bb[['xmax']],
    ymin = bb[['ymin']], ymax = bb[['ymax']],
    resolution = res, crs = 'EPSG:4326'
  )
  terra::values(r) <- stats::runif(terra::ncell(r))
  for (i in seq_len(smooth_iter)) {
    r <- terra::focal(r, w = 5, fun = 'mean', na.policy = 'omit', na.rm = TRUE)
  }
  v <- terra::values(r)
  v <- (v - min(v, na.rm = TRUE)) / (max(v, na.rm = TRUE) - min(v, na.rm = TRUE))
  terra::values(r) <- range[1] + v * (range[2] - range[1])
  r
}
