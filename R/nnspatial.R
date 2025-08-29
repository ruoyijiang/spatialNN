#' Precompute pairwise distances between two cluster sets (Visium)
#'
#' @param obj_input Seurat object with VISIUM image
#' @param cl_a character vector of cluster labels for set A
#' @param cl_b character vector of cluster labels for set B
#' @param px_to_um numeric, pixel-to-micron conversion
#' @param cluster_col, pick the ident to subset on, if NULL, will use Idents(obj_input)
#' @return list with D_um (matrix), a_ids (cell ids for A), b_ids (cell ids for B)
#' @export
precompute_D <- function(obj_input, cl_a, cl_b,
                         px_to_um = 0.12028,
                         cluster_col = "SUBSET") {
  fov  <- Seurat::Images(obj_input)[1]
  cent <- SeuratObject::GetTissueCoordinates(obj_input[[fov]], boundary = "centroids")

  # pick cluster labels
  clusters <- if (is.null(cluster_col)) {
    as.character(Seurat::Idents(obj_input))
  } else {
    if (!cluster_col %in% colnames(obj_input@meta.data)) {
      stop(sprintf("cluster_col '%s' not found in meta.data. Available: %s",
                   cluster_col, paste(colnames(obj_input@meta.data), collapse = ", ")))
    }
    as.character(obj_input@meta.data[[cluster_col]])
  }

  cent$cluster <- clusters[match(cent$cell, rownames(obj_input@meta.data))]

  pts <- sf::st_as_sf(cent, coords = c("x","y"), remove = FALSE); sf::st_crs(pts) <- NA
  labs <- as.character(pts$cluster)
  a <- pts[labs %in% cl_a, ]
  b <- pts[labs %in% cl_b, ]
  stopifnot(nrow(a) > 0, nrow(b) > 0)

  D_px <- sf::st_distance(a, b)
  D_um <- as.matrix(D_px); mode(D_um) <- "numeric"
  D_um <- D_um * px_to_um

  list(D_um = D_um, a_ids = a$cell, b_ids = b$cell)
}

#' Nearest neighbors between A and B (exact or bootstrap)
#'
#' @param pre list returned by precompute_D()
#' @param n_a,n_b optional sizes when bootstrap = TRUE
#' @param bootstrap logical; if TRUE, sample rows/cols with replacement
#' @param replace_a,replace_b logical; sampling behavior
#' @return list with cell_id_a, nearest_b, nearest_dist_um (and indices in bootstrap)
#' @export
compute_nearest <- function(input_D_um, n_a = NULL, n_b = NULL, bootstrap = FALSE) {
 if(!bootstrap){
      a_ids <- input_D_um$a_ids
      b_ids <- input_D_um$b_ids
      D <- input_D_um$D_um
  } else {
      a_ids  <- sample.int(nrow(input_D_um$D_um), n_a, replace = bootstrap)
      b_ids  <- sample.int(ncol(input_D_um$D_um), n_b, replace = bootstrap)
      D <- input_D_um$D_um[a_ids, b_ids, drop = FALSE]
  }

  # check
  stopifnot(is.matrix(D), length(a_ids) == nrow(D), length(b_ids) == ncol(D))

  nn_idx <- max.col(-D, ties.method = "first")
  nn_dist_um  <- rowMins(D)

  list(
    cell_id_a       = a_ids,
    nearest_b       = b_ids[nn_idx],
    nearest_dist_um = as.numeric(nn_dist_um)
  )
}

#' Bootstrap summary of nearest-neighbor distances
#'
#' @param pre list from precompute_D()
#' @param n_a,n_b sample sizes for rows/cols
#' @param R number of bootstrap replicates
#' @param seed random seed
#' @param stat "median" or "mean"
#' @return list with bootstrap vector and 90% CI (5th/95th percentiles)
#' @export
boot_nn_stat <- function(pre, n_a, n_b, R = 1000, seed = 1, stat = c("median","mean")) {
  stat <- match.arg(stat)
  set.seed(seed)
  vals <- replicate(R, {
    out <- compute_nearest(pre, n_a = n_a, n_b = n_b, bootstrap = TRUE)
    if (stat == "median") stats::median(out$nearest_dist_um) else stats::mean(out$nearest_dist_um)
  })
  list(boot = vals, ci90 = stats::quantile(vals, c(0.05, 0.95)))
}

#' Fraction below distance thresholds
#'
#' @param out list from compute_nearest() (non-bootstrap), must contain nearest_dist_um
#' @param obj Seurat object to compute denominator from (SUBSET matches subset_id)
#' @param cluster_col, pick the ident to subset on, if NULL, will use Idents(obj_input)
#' @param subset_id character label for denominator
#' @param step_um,max_um numeric
#' @return tibble with thresholds and fractions
#' @export
compute_frac_vs_threshold <- function(out, obj, subset_id = "T cell", cluster_col = "SUBSET",
                                      step_um = 1,
                                      max_um = 20) {
  stopifnot("nearest_dist_um" %in% names(out))

  # choose denominator labels
  clusters <- if (is.null(cluster_col)) {
    as.character(Seurat::Idents(obj))
  } else {
    if (!cluster_col %in% colnames(obj@meta.data)) {
      stop(sprintf("cluster_col '%s' not found in meta.data.", cluster_col))
    }
    as.character(obj@meta.data[[cluster_col]])
  }
  denom <- sum(clusters == subset_id, na.rm = TRUE)
  if (denom == 0) {
    stop(sprintf(
      "No cells found with label '%s' in '%s'.",
      subset_id,
      ifelse(is.null(cluster_col), "Idents(obj)", cluster_col)
    ))
  }

  if (is.null(max_um)) max_um <- ceiling(max(out$nearest_dist_um, na.rm = TRUE))
  thresholds <- seq(0, max_um, by = step_um)

  tibble::tibble(distance_threshold_um = thresholds) |>
    dplyr::mutate(
      n_below = sapply(distance_threshold_um,
                       function(t) sum(out$nearest_dist_um < t, na.rm = TRUE)),
      denom = denom,
      fraction = pmin(n_below, denom) / denom,
      percent = 100 * fraction
    )
}
