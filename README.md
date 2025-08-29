# spatialNN
This R package computes distances with/without bootstrapping between cells identified from spatial transcriptome Seurat object

# spatialNN — nearest neighbors between spatial clusters

Utilities to (1) precompute pairwise distances between two Seurat cluster sets, (2) extract exact nearest neighbors, (3) run bootstrap summaries of nearest-neighbor distances, and (4) summarize fractions of cells within distance thresholds.

## Functions

- `precompute_D(obj_input, cl_a, cl_b, px_to_um = 0.12028, cluster_col = "SUBSET")`  
  Build a microns distance matrix `D_um` between clusters `cl_a` and `cl_b` in the same FOV, and store the corresponding cell IDs for A (rows) and B (columns). `cluster_col` chooses which metadata column defines clusters (use `NULL` to fall back to `Seurat::Idents(obj_input)`).

- `compute_nearest(pre, bootstrap = FALSE, n_a = NULL, n_b = NULL, replace_a = TRUE, replace_b = TRUE)`  
  **Exact mode (default):** returns `cell_id_a`, `nearest_b`, and `nearest_dist_um` for every A cell (nearest over all B).  
  **Bootstrap mode:** samples rows/cols (with replacement by default) and returns distances and chosen IDs computed on the sampled submatrix. Uses `matrixStats::rowMins` if available, otherwise `apply(..., 1, min)`.

- `boot_nn_stat(pre, n_a, n_b, R = 1000, seed = 1, stat = c("median","mean"))`  
  Run `R` bootstrap replicates of the summary statistic (median/mean) of nearest-neighbor distances via `compute_nearest(pre, bootstrap=TRUE)`; returns the bootstrap vector and a simple 90% CI (5th/95th percentiles).

- `compute_frac_vs_threshold(out, obj, subset_id = "T cell", step_um = 1, max_um = 20, cluster_col = "SUBSET")`  
  Return a tibble of fraction/percent of A cells whose nearest B is within each distance threshold, using the number of `subset_id` cells in `obj` as the denominator. `cluster_col` controls which metadata column defines the denominator group; set to `NULL` to use `Idents(obj)`.

---

## Install / load (for now)

If you have the standalone R file (e.g., `R/nnspatial.R`), source it directly:

```r
# one-time: install dependencies
install.packages(c("sf","Seurat","SeuratObject","matrixStats","tibble","dplyr"))

# load functions
source("path/to/R/nnspatial.R")
```

If you have the package folder (e.g., `~/dev/spatialNN`), you can hot-load with devtools:

```r
install.packages("devtools")
devtools::load_all("~/dev/spatialNN")
```

---

## Quick start

```r
library(Seurat)
library(sf)
library(tibble)
library(dplyr)

# Assume: Seurat Spatial object `obj` with the cluster labels in meta.data$SUBSET
# Choose two cluster sets
cl_a <- c("T cell")
cl_b <- c("B cell")

# 1) Precompute distances (microns) and capture cell IDs
pre <- precompute_D(
  obj_input = obj,
  cl_a = cl_a,
  cl_b = cl_b,
  px_to_um = 0.12028,        # adjust to your slide calibration
  cluster_col = "SUBSET"     # or NULL to use Idents(obj)
)

# 2) Exact nearest neighbors for all A cells
res_exact <- compute_nearest(pre)  # list: cell_id_a, nearest_b, nearest_dist_um
# Convert to tibble for easy viewing
res_exact_tbl <- tibble(
  cell_id_a       = res_exact$cell_id_a,
  nearest_b       = res_exact$nearest_b,
  nearest_dist_um = res_exact$nearest_dist_um
)
res_exact_tbl %>% slice_head(n = 6)

# 3) Bootstrap: sample subsets of A and B, compute distances on the sampled submatrix
set.seed(1)
res_boot <- compute_nearest(
  pre,
  bootstrap = TRUE,
  n_a = 500,            # sample 500 A cells (with replacement by default)
  n_b = 800
)
summary(res_boot$nearest_dist_um)

# 4) Bootstrap distribution of the median distance
bstat <- boot_nn_stat(
  pre,
  n_a = 500,
  n_b = 800,
  R = 1000,
  seed = 1,
  stat = "median"
)
str(bstat)
quantile(bstat$boot, c(0.05, 0.5, 0.95))

# 5) Fraction vs threshold curve (denominator = count of 'T cell' in obj)
frac_tbl <- compute_frac_vs_threshold(
  out = res_exact,
  obj = obj,
  subset_id = "T cell",
  step_um = 2,
  max_um = 50,
  cluster_col = "SUBSET"
)
frac_tbl %>% slice_head(n = 10)
```

---

## Examples: choosing different cluster columns

Use Seurat's default identities:
```r
# Make sure Idents(obj) holds the labels you want
Idents(obj) <- obj$seurat_clusters
pre <- precompute_D(obj, cl_a = c("0","1"), cl_b = c("2","3"), cluster_col = NULL)
res <- compute_nearest(pre)  # exact NN on these groups
```

Use a custom column, e.g., `meta.data$celltype_major`:
```r
pre <- precompute_D(
  obj, cl_a = c("CD8 T","CD4 T"), cl_b = c("Myeloid"),
  cluster_col = "celltype_major"
)
res <- compute_nearest(pre)
```

---

## Notes / tips

- Distances: `sf::st_distance` returns units; the code coerces to numeric and multiplies by `px_to_um` to produce microns. Verify your calibration constant per slide/scanner.
- Single FOV: `precompute_D` uses the first image from `Seurat::Images(obj_input)`. If your object has multiple images, consider looping over them.
- Ties: nearest-neighbor indices are computed with `max.col(-D, ties.method = "first")`.
- Performance: For very large A/B sets, bootstrap on submatrices is often faster for inference tasks.

