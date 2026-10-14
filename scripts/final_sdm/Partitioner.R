#' unified spatial-environmental partitioner for clustered, bimodal SDM data
#'
#' Sourceable extract of scripts/final_sdm/unified_partitioner.qmd (which carries the
#' full derivation, diagnostics and plots). This file is the *logic only*: the
#' partition_* helpers compose into build_partition(). Base R + sf + FNN + twosamples
#' + stats (CAST is optional, only for print()/plot() on the knndm-classed fold objects).
#'
#' INPUT CONTRACT
#'   x   : sf of PRESENCES carrying an integer .id column AND the environmental PCA
#'         score columns (named by `pcs`, e.g. PC1..PCn). On the data machine build these
#'         by extracting predictors (terra::extract on the rastReader stack) and running
#'         prcomp(center, scale), keeping axes to a cumulative-variance target - the Step 0
#'         code in the .qmd. Occurrence == 1.
#'   bg  : sf of background ABSENCES (Occurrence == 0); no env needed here (attached at fit).
#'   pcs : character vector naming the PCA columns on `x`.
#'   p   : parameter list; use partition_params() and override as needed.
#'
#' OUTPUT (build_partition): a list with the presence-level structure (train_id/test_id,
#' components, geo_d/env_d, roles, split_coverage, fold_meta, W) and the caret-ready
#' class-aware frame: combined, train_frame, and the nndm-style index triple
#' indx_train / indx_test / indx_exclude (positional into train_frame).
#'
#' HANDOFF TO THE FIT (so CAST::aoaSurfaces() derives its AOA threshold from THIS partition):
#'   part <- build_partition(x, bg, pcs, partition_params())
#'   ctrl <- caret::trainControl(index = part$indx_train, indexOut = part$indx_test, ...)
#'   # fit with ctrl; then CAST::trainDI(model) inside aoaSurfaces() picks up these folds.
#' indx_exclude (extrap, validation-only) is already removed from indx_train.

if (!requireNamespace("sf",         quietly = TRUE)) stop("Partitioner.R requires sf")
if (!requireNamespace("FNN",        quietly = TRUE)) stop("Partitioner.R requires FNN")
if (!requireNamespace("twosamples", quietly = TRUE)) stop("Partitioner.R requires twosamples")

#' Default partitioner parameters (mirrors params_part in the .qmd). Override fields as needed.
partition_params <- function(
    resolution     = 90,
    connect_m      = 500,
    superclust_m   = 40000,
    train_target   = c(0.70, 0.80),
    k              = 10,
    purpose        = "mapping",       # discovery | mapping | projection
    pca_var_target = 0.90,
    min_split_n    = 12,
    extrap_frac    = 0.8,
    seed           = 1021) {
  list(resolution = resolution, connect_m = connect_m, superclust_m = superclust_m,
       train_target = train_target, k = k, purpose = purpose, pca_var_target = pca_var_target,
       min_split_n = min_split_n, extrap_frac = extrap_frac, seed = seed)
}

# ----------------------------------------------------------------------------
# Step 1 - connected components (single-linkage graph at connect_m)
# ----------------------------------------------------------------------------
partition_components <- function(x, connect_m) {
  # Connected components of the within-connect_m graph
  within <- sf::st_is_within_distance(x, dist = connect_m)  # sparse list kept for the edge plot
  d      <- sf::st_distance(x); units(d) <- NULL
  memb   <- stats::cutree(stats::hclust(stats::as.dist(d), method = "single"), h = connect_m)
  x$component <- factor(memb)
  attr(x, "adj") <- within  
  x
}

# ----------------------------------------------------------------------------
# Step 2 - inter-component minimum geographic distance
# ----------------------------------------------------------------------------
partition_geo_dist <- function(x, hull = FALSE) {
  lev   <- levels(x$component)
  geoms <- do.call(c, lapply(lev, function(l)                      # one (MULTI)geometry per component
             sf::st_union(sf::st_geometry(x)[x$component == l])))
  parts <- sf::st_sf(component = lev, geometry = geoms)
  if (hull) parts <- sf::st_convex_hull(parts)
  d <- sf::st_distance(parts)                                      # min distance between geometries
  units(d) <- NULL                                                 # strip units (CAST idiom)
  dimnames(d) <- list(lev, lev)
  d
}

# ----------------------------------------------------------------------------
# Step 2a - inter-component environmental (PCA-medoid) distance
# ----------------------------------------------------------------------------
partition_env_dist <- function(x, pcs) {
  sc  <- sf::st_drop_geometry(x)
  lev <- levels(x$component)
  m   <- t(vapply(lev, function(l)                                 # per-component PCA medoid (column medians)
            apply(sc[sc$component == l, pcs, drop = FALSE], 2, stats::median, na.rm = TRUE),
            numeric(length(pcs))))
  rownames(m) <- lev; colnames(m) <- pcs
  as.matrix(stats::dist(m))
}

# ----------------------------------------------------------------------------
# Step 3 - component roles (extrapolation / splittable / atomic) via per-component DI/AOA
# ----------------------------------------------------------------------------
partition_splittable <- function(x, pcs, purpose, min_split_n = 12, extrap_frac = 0.8) {
  xd    <- sf::st_drop_geometry(x)
  comps <- levels(x$component)
  P     <- as.matrix(xd[pcs])

  # REALIZED extrapolation flag (shared by Steps 4-5): CAST's DI/AOA per component. Hold each
  # component out WHOLE; the env-NN p95 of the rest-of-data is the AOA threshold (`thr`), each
  # held-out point's env-NN distance is its DI (`enn`), and frac_beyond is the share outside
  # the AOA (DI > thr). frac_beyond ~1 ⇒ an env-distinct segment (an outlier population or a
  # resolved mode); a fraction only means a self-exclusion gap. Direct and tuning-light, and it
  # drives allocation in Steps 4-5: extrap segments are kept out of training and spread across
  # the holdouts. (trainDI()/aoa() give the same DI/threshold for a global cross-check.)
  seg_extrap <- do.call(rbind, lapply(comps, function(cmp) {
    te <- P[xd$component == cmp, , drop = FALSE]
    tr <- P[xd$component != cmp, , drop = FALSE]
    tnn <- c(FNN::knn.dist(tr, k = 1)) # rest-of-data sample-to-sample NN
    thr <- stats::quantile(tnn, 0.95) # AOA threshold
    enn <- c(FNN::knnx.dist(data = tr, query = te, k = 1))         # held-out component's DI
    data.frame(component = cmp, n_beyond = sum(enn > thr), frac_beyond = mean(enn > thr),
               stringsAsFactors = FALSE)
  }))

  # one row per component, in factor-level order
  out <- data.frame(
    component    = comps,
    supercluster = xd$supercluster[match(comps, as.character(xd$component))],
    n            = as.integer(table(xd$component)[comps]),
    stringsAsFactors = FALSE)
  out <- out[, c("component", "supercluster", "n")]
  mi  <- match(out$component, seg_extrap$component)
  out$n_beyond    <- seg_extrap$n_beyond[mi]
  out$frac_beyond <- seg_extrap$frac_beyond[mi]
  out$extrap     <- purpose == "mapping" & out$frac_beyond >= extrap_frac
  out$splittable <- purpose == "mapping" & !out$extrap & out$n >= min_split_n
  # `purpose` is scalar, so branch it with if/else (a scalar test in ifelse() would collapse the
  # result to length 1 and recycle); only the per-component flags go through vectorized ifelse().
  out$role <-
    if (purpose == "discovery")  rep("atomic (supercluster = unit)", nrow(out))
    else if (purpose == "projection") rep("defer to env thinning", nrow(out))
    else ifelse(out$extrap,     "extrapolation (spread, held from training)",
         ifelse(out$splittable, "splittable (interpolation)",
                                "atomic (interpolation)"))
  out
}

# ----------------------------------------------------------------------------
# knndm distance quantities: Gj / Gjstar / Gij (proxy) / W (twosamples::wass_stat)
# ----------------------------------------------------------------------------
partition_knndm_dists <- function(P, clusters, extrap) {
  Gj     <- c(FNN::knn.dist(P, k = 1))                 # sample-to-sample LOO NN
  Gjstar <- numeric(nrow(P))                           # CV NN: point -> its fold's training set
  for (f in sort(unique(clusters))) {
    te <- which(clusters == f); tr <- which(clusters != f & !extrap)
    Gjstar[te] <- c(FNN::knnx.dist(data = P[tr, , drop = FALSE], query = P[te, , drop = FALSE], k = 1))
  }
  Gij <- Gj
  list(Gj = Gj, Gjstar = Gjstar, Gij = Gij, W = twosamples::wass_stat(Gjstar, Gij))
}

# ----------------------------------------------------------------------------
# Step 4 - outer test/train split (soft train window, NNDM coverage)
# ----------------------------------------------------------------------------
partition_split <- function(x, spl, pcs, train_target, seed = 1) {
  set.seed(seed)
  xd <- sf::st_drop_geometry(x)
  mi <- match(as.character(xd$component), spl$component)
  xd$role <- spl$role[mi]; xd$extrap <- spl$extrap[mi]; xd$n_beyond <- spl$n_beyond[mi]
  N  <- nrow(xd)
  test_frac   <- 1 - mean(train_target)
  target_test <- round(N * test_frac)

  # k-means-medoid sampler: representative (not extreme) spread over all env axes
  kmed_pick <- function(M, n_pick) {
    n <- nrow(M)
    if (n_pick >= n) return(seq_len(n))
    if (n_pick < 1)  return(integer(0))
    n_pick <- min(n_pick, nrow(unique(M)))     # kmeans needs centers <= distinct points
    km <- suppressWarnings(stats::kmeans(M, centers = n_pick, nstart = 5, iter.max = 50))
    vapply(seq_len(n_pick), function(k) {
      idx <- which(km$cluster == k)
      idx[which.min(rowSums((M[idx, , drop = FALSE] -
                             matrix(km$centers[k, ], length(idx), ncol(M), byrow = TRUE))^2))]
    }, integer(1))
  }

  test_id <- integer(0)

  # 1. extrapolation segments: SPREAD their points — a representative test_frac slice to the
  #    stringent test, the rest stay in TRAIN so Step 5's inner folds have extrapolation to
  #    spread too. Env-distinct points are never dumped wholesale to either side.
  ex_ids <- xd$.id[xd$extrap]
  if (length(ex_ids)) {
    M <- as.matrix(xd[match(ex_ids, xd$.id), pcs])
    test_id <- c(test_id, ex_ids[kmed_pick(M, round(length(ex_ids) * test_frac))])
  }

  # 2. interpolation splittable: within-component representative test subset
  for (cmp in spl$component[spl$role == "splittable (interpolation)"]) {
    idx  <- xd$.id[xd$component == cmp]
    M    <- as.matrix(xd[match(idx, xd$.id), pcs])
    pick <- idx[kmed_pick(M, round(length(idx) * test_frac))]
    test_id <- c(test_id, pick)
  }

  # 3. interpolation atomic: add whole units until the budget is met, HARDEST first (most
  #    points beyond rest-of-data p95), so the residual budget is spent deliberately rather
  #    than by luck of a random draw.
  atomic <- spl[spl$role == "atomic (interpolation)", , drop = FALSE]
  atomic_order <- atomic$component[order(-atomic$n_beyond)]      # hardest (most beyond p95) first
  for (cmp in atomic_order) {
    if (length(test_id) >= target_test) break
    test_id <- c(test_id, xd$.id[xd$component == cmp])
  }

  train_id <- setdiff(xd$.id, test_id)

  # NNDM-style coverage: test->train NN env distances vs train->train NN distances
  Ptr <- as.matrix(xd[match(train_id, xd$.id), pcs])
  Pte <- as.matrix(xd[match(test_id,  xd$.id), pcs])
  test_nn  <- c(FNN::knnx.dist(data = Ptr, query = Pte, k = 1))   # test-to-sample NN
  train_nn <- c(FNN::knn.dist(Ptr, k = 1))                        # sample-to-sample LOO NN

  set     <- ifelse(xd$.id %in% test_id, "test", "train")
  comps_u <- levels(xd$component)
  tr_n    <- table(xd$component[set == "train"]); te_n <- table(xd$component[set == "test"])
  by_unit <- data.frame(
    component = comps_u,
    role      = xd$role[match(comps_u, as.character(xd$component))],
    train     = as.integer(tr_n[comps_u]),
    test      = as.integer(te_n[comps_u]),
    stringsAsFactors = FALSE)
  by_unit$train[is.na(by_unit$train)] <- 0L
  by_unit$test[is.na(by_unit$test)]   <- 0L

  list(
    train_id = train_id, test_id = test_id, by_unit = by_unit,
    achieved_train_prop = length(train_id) / N, target_test = target_test,
    nn = list(train = train_nn, test = test_nn),
    coverage = list(
      frac_beyond_p95 = mean(test_nn > stats::quantile(train_nn, 0.95)),
      frac_beyond_max = mean(test_nn > max(train_nn)))
  )
}

# ----------------------------------------------------------------------------
# Step 5 - inner CV folds (cut & re-represent); returns a classed knndm object
# ----------------------------------------------------------------------------
partition_folds <- function(x, spl, pcs, train_id, k = 10, seed = 1) {
  set.seed(seed)
  xt <- sf::st_drop_geometry(x)
  xt <- xt[xt$.id %in% train_id, , drop = FALSE]
  mi <- match(as.character(xt$component), spl$component)
  xt$role <- spl$role[mi]; xt$extrap <- spl$extrap[mi]
  xt$row <- seq_len(nrow(xt))                       # positional index into the training frame
  k <- as.integer(k)
  fold <- integer(nrow(xt))

  # interpolation splittable: spread each component's points across the k folds (re-represent)
  for (cmp in unique(xt$component[!xt$extrap & xt$role == "splittable (interpolation)"])) {
    idx <- which(xt$component == cmp & !xt$extrap)
    M   <- as.matrix(xt[idx, pcs])
    kk  <- min(k, nrow(unique(M)))            # cap at distinct env points
    cl  <- if (kk <= 1)          rep(1L, length(idx))
           else if (kk >= nrow(M)) seq_len(nrow(M))   # each point its own cluster (avoids kmeans centers==nrow edge)
           else suppressWarnings(stats::kmeans(M, centers = kk, nstart = 5, iter.max = 50)$cluster)
    fold[idx] <- ((cl - 1) %% k) + 1
  }

  # interpolation atomic: whole unit -> one fold, greedily load-balanced, largest unit first
  sel   <- !xt$extrap & xt$role != "splittable (interpolation)"
  wcnt  <- table(droplevels(xt$component[sel]))
  whole <- names(wcnt)[order(-as.integer(wcnt))]
  load  <- tabulate(fold[fold > 0], nbins = k)
  for (cmp in whole) {
    f <- which.min(load); idx <- which(xt$component == cmp & !xt$extrap)
    fold[idx] <- f; load[f] <- load[f] + length(idx)
  }

  # extrapolation points: round-robin across folds as VALIDATION-ONLY. They are given a
  # validation fold but are excluded from every fold's training (below), so each fold gets a
  # comparable slice of genuinely-unseen points instead of one fold absorbing a whole segment.
  eload <- integer(k)
  for (i in which(xt$extrap)) { f <- which.min(eload); fold[i] <- f; eload[f] <- eload[f] + 1 }

  # = extrap rows assigned to OTHER folds — neither trained nor validated in fold f.
  indx_test    <- lapply(seq_len(k), function(f) xt$row[fold == f])
  indx_train   <- lapply(seq_len(k), function(f) xt$row[fold != f & !xt$extrap])
  indx_exclude <- lapply(seq_len(k), function(f) xt$row[fold != f &  xt$extrap])

  fold_meta <- do.call(rbind, lapply(seq_len(k), function(f) {
    in_val <- fold == f
    tr <- as.matrix(xt[fold != f & !xt$extrap, pcs]); te <- as.matrix(xt[in_val, pcs])
    tnn <- c(FNN::knn.dist(tr, k = 1))                          # fold-train sample-to-sample NN
    enn <- c(FNN::knnx.dist(data = tr, query = te, k = 1))      # CV-distances (validation -> fold-train)
    data.frame(fold = f, n_val = nrow(te), val_prop = nrow(te) / nrow(xt),
               n_extrap_val = sum(xt$extrap[in_val]),
               beyond_p95 = mean(enn > stats::quantile(tnn, 0.95)),
               beyond_max = mean(enn > max(tnn)))
  }))

  d <- partition_knndm_dists(as.matrix(xt[pcs]), fold, xt$extrap)
  res <- list(
    clusters   = fold,                                  # per-training-row fold id (knndm$clusters)
    indx_train = indx_train, indx_test = indx_test, indx_exclude = indx_exclude,
    Gij = d$Gij, Gj = d$Gj, Gjstar = d$Gjstar, W = d$W,
    method = "single-linkage graph", q = "role-based (cut/re-represent/spread)",
    space = "feature", dist_space = "feature",
    k = k, fold = setNames(fold, xt$.id), extrap_id = xt$.id[xt$extrap], fold_meta = fold_meta)
  class(res) <- c("knndm", "list")
  res
}

# ----------------------------------------------------------------------------
# Step 5b - leave-whole-segment-out folds (classed knndm object)
# ----------------------------------------------------------------------------
partition_folds_lopo <- function(x, spl, pcs, train_id, k = 10, seed = 1) {
  set.seed(seed)
  xt <- sf::st_drop_geometry(x)
  xt <- xt[xt$.id %in% train_id, , drop = FALSE]
  mi <- match(as.character(xt$component), spl$component)
  xt$role <- spl$role[mi]; xt$extrap <- spl$extrap[mi]
  xt$row <- seq_len(nrow(xt))                    # positional index into the training frame
  k <- as.integer(k)
  fold <- integer(nrow(xt))

  # interpolation segments: whole segment -> one fold, greedy size-balance (true leave-segment-out)
  cnt  <- table(droplevels(xt$component[!xt$extrap]))
  csz  <- data.frame(component = names(cnt), n = as.integer(cnt), stringsAsFactors = FALSE)
  csz  <- csz[order(-csz$n), , drop = FALSE]                     # largest segment first
  load <- integer(k)
  for (i in seq_len(nrow(csz))) {
    f <- which.min(load); idx <- which(xt$component == csz$component[i] & !xt$extrap)
    fold[idx] <- f; load[f] <- load[f] + csz$n[i]
  }

  # extrapolation points: round-robin across folds as VALIDATION-ONLY (excluded from training
  # below), so every fold gets a comparable slice of genuinely-unseen points.
  eload <- integer(k)
  for (i in which(xt$extrap)) { f <- which.min(eload); fold[i] <- f; eload[f] <- eload[f] + 1 }

  # CAST triple (nndm vocabulary): training EXCLUDES all extrap rows; validation = fold members
  # (whole interp segment + extrap dealt there); indx_exclude[[f]] = extrap rows of OTHER folds.
  indx_test    <- lapply(seq_len(k), function(f) xt$row[fold == f])
  indx_train   <- lapply(seq_len(k), function(f) xt$row[fold != f & !xt$extrap])
  indx_exclude <- lapply(seq_len(k), function(f) xt$row[fold != f &  xt$extrap])

  fold_meta <- do.call(rbind, lapply(seq_len(k), function(f) {
    in_val <- fold == f
    tr <- as.matrix(xt[fold != f & !xt$extrap, pcs]); te <- as.matrix(xt[in_val, pcs])
    tnn <- c(FNN::knn.dist(tr, k = 1))                          # fold-train sample-to-sample NN
    enn <- c(FNN::knnx.dist(data = tr, query = te, k = 1))      # CV-distances (validation -> fold-train)
    data.frame(fold = f, n_val = nrow(te), val_prop = nrow(te) / nrow(xt),
           n_comp = length(unique(xt$component[in_val & !xt$extrap])),
           n_superc = length(unique(xt$supercluster[in_val])),
           n_extrap_val = sum(xt$extrap[in_val]),
           beyond_p95 = mean(enn > stats::quantile(tnn, 0.95)),
           beyond_max = mean(enn > max(tnn)))
  }))

  # knndm-canonical classed object (same as partition_folds), so CAST's print/plot dispatch here too
  d <- partition_knndm_dists(as.matrix(xt[pcs]), fold, xt$extrap)
  res <- list(
    clusters   = fold,
    indx_train = indx_train, indx_test = indx_test, indx_exclude = indx_exclude,
    Gij = d$Gij, Gj = d$Gj, Gjstar = d$Gjstar, W = d$W,
    method = "segment-LOPO (whole segment)", q = "leave-segment-out",
    space = "feature", dist_space = "feature",
    k = k, fold = setNames(fold, xt$.id), extrap_id = xt$.id[xt$extrap], fold_meta = fold_meta)
  class(res) <- c("knndm", "list")
  res
}

# ----------------------------------------------------------------------------
# Step 6 - class-aware assembly (presences + absences -> caret-ready frame)
# ----------------------------------------------------------------------------
partition_classes <- function(pts, bg, split, folds, train_target, k, seed = 1) {
  set.seed(seed)
  test_prop <- 1 - mean(train_target)

  # presences: structured assignment from Steps 4-5 (test presences get fold = NA).
  # extrap presences carry their validation fold but are flagged so they are NEVER put in
  # any fold's training (they are validation-only, matching the Step 5 spread).
  pres <- data.frame(
    id     = paste0("pres_", pts$.id),
    class  = 1L,
    set    = ifelse(pts$.id %in% split$test_id, "test", "train"),
    fold   = unname(folds$fold[as.character(pts$.id)]),
    extrap = pts$.id %in% folds$extrap_id,
    stringsAsFactors = FALSE)

  # absences: random outer split, then balanced random folds within the train side
  nb  <- nrow(bg)
  aid <- paste0("abs_", seq_len(nb))
  a_set <- rep("train", nb); a_set[sample.int(nb, round(nb * test_prop))] <- "test"
  a_fold <- rep(NA_integer_, nb)
  tr <- which(a_set == "train")
  a_fold[tr] <- sample(rep_len(seq_len(k), length(tr)))      # balanced across folds
  abs <- data.frame(id = aid, class = 0L, set = a_set, fold = a_fold, extrap = FALSE,
                    stringsAsFactors = FALSE)

  combined <- rbind(pres, abs)
  tf <- combined[combined$set == "train", ]; tf$row <- seq_len(nrow(tf))
  # CAST triple (nndm vocabulary) carried through to the caret-ready frame.
  indx_test    <- lapply(seq_len(k), function(f) tf$row[which(tf$fold == f)])
  indx_train   <- lapply(seq_len(k), function(f) tf$row[which(tf$fold != f & !tf$extrap)])
  indx_exclude <- lapply(seq_len(k), function(f) tf$row[which(tf$fold != f &  tf$extrap)])

  list(combined = combined, train_frame = tf,
       indx_train = indx_train, indx_test = indx_test, indx_exclude = indx_exclude,
       test_ids = combined$id[combined$set == "test"])
}

# ----------------------------------------------------------------------------
# Assemble - one call composing Steps 1-6 into the shippable artifact
# ----------------------------------------------------------------------------
build_partition <- function(x, bg, pcs, p) {   # p = params_part
  x     <- partition_components(x, p$connect_m)
  geo_d <- partition_geo_dist(x)
  x$supercluster <- factor(
    cutree(hclust(as.dist(geo_d), method = "single"), h = p$superclust_m)[as.character(x$component)])
  env_d   <- partition_env_dist(x, pcs)
  spl     <- partition_splittable(x, pcs, p$purpose, p$min_split_n, p$extrap_frac)
  split   <- partition_split(x, spl, pcs, p$train_target, p$seed)
  folds   <- partition_folds(x, spl, pcs, split$train_id, p$k, p$seed)
  classes <- partition_classes(x, bg, split, folds, p$train_target, p$k, p$seed)

  list(
    # presence-level structure
    train_id   = split$train_id,
    test_id    = split$test_id,
    components = sf::st_drop_geometry(x)[, c(".id", "component", "supercluster")],
    geo_d = geo_d, env_d = env_d, splittable = spl,
    split_coverage = split$coverage, fold_meta = folds$fold_meta, W = folds$W,
    folds = folds,                        # the classed knndm object (CAST print/plot dispatch)
    # class-aware frame for fitting (caret-ready)
    combined   = classes$combined,        # every point: id, class, set, fold
    train_frame = classes$train_frame,
    indx_train   = classes$indx_train,    # positional into train_frame, for caret::trainControl(index=)
    indx_test    = classes$indx_test,     # ...and indexOut=
    indx_exclude = classes$indx_exclude,  # nndm-style: held from training, not validated (extrap)
    test_ids     = classes$test_ids,
    meta = p[c("connect_m", "superclust_m", "train_target", "k", "purpose",
               "pca_var_target", "min_split_n", "extrap_frac", "seed")]
  )
}
