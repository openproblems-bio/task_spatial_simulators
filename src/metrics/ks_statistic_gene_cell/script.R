requireNamespace("anndata", quietly = TRUE)
requireNamespace("edgeR", quietly = TRUE)
requireNamespace("ks", quietly = TRUE)
library(Matrix)
library(matrixStats)

## VIASH START
dir <- "work/66/65ddd5a686b4c712172ff7aa70cce8/_viash_par"
par <- list(
  input_spatial_dataset = paste0(
    dir,
    "/input_spatial_dataset_1/dataset_sp.h5ad"
  ),
  input_singlecell_dataset = paste0(
    dir,
    "/input_singlecell_dataset_1/dataset_sc.h5ad"
  ),
  input_simulated_dataset = paste0(
    dir,
    "/input_simulated_dataset_1/spatialsimbench_mobnew.negative_normal.generate_sim_spatialcluster.output_sp.h5ad"
  ),
  output = "output.h5ad"
)
meta <- list(
  name = "ks_statistic"
)
## VIASH END

cat("Reading input files\n")
input_spatial_dataset <- anndataR::read_h5ad(par[["input_spatial_dataset"]])
input_singlecell_dataset <- anndataR::read_h5ad(par[[
  "input_singlecell_dataset"
]])
input_simulated_dataset <- anndataR::read_h5ad(par[["input_simulated_dataset"]])

real_counts <- input_spatial_dataset$layers[["counts"]]
sim_counts <- input_simulated_dataset$layers[["counts"]]

as_finite_kde_input <- function(x) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (is.matrix(x)) {
    x <- matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x))
    return(x[apply(is.finite(x), 1, all), , drop = FALSE])
  }

  x <- as.numeric(x)
  x[is.finite(x)]
}

ks_penalty_result <- function(reason) {
  warning(reason, " Returning worst-case KS statistic of 1.")
  list(zstat = 1, tstat = 1)
}

empirical_ks_statistic <- function(x1, x2) {
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)

  if (length(x1) < 2 || length(x2) < 2) {
    return(NA_real_)
  }

  result <- tryCatch(
    suppressWarnings(stats::ks.test(x1 = x1, x2 = x2, exact = FALSE)),
    error = function(e) NULL
  )

  if (is.null(result) || !is.finite(as.numeric(result$statistic))) {
    return(NA_real_)
  }

  as.numeric(result$statistic)
}

sliced_ks_statistic <- function(x1, x2) {
  n_dim <- ncol(x1)
  directions <- diag(n_dim)

  if (n_dim > 1) {
    directions <- rbind(
      directions,
      rep(1, n_dim),
      c(rep(1, n_dim - 1), -1)
    )
  }

  stats <- apply(directions, 1, function(direction) {
    direction <- direction / sqrt(sum(direction^2))
    empirical_ks_statistic(
      as.numeric(x1 %*% direction),
      as.numeric(x2 %*% direction)
    )
  })

  if (all(!is.finite(stats))) {
    return(NA_real_)
  }

  max(stats, na.rm = TRUE)
}

fallback_ks_result <- function(statistic, reason) {
  statistic <- as.numeric(statistic)
  warning(reason, " Falling back to empirical KS statistic.")
  list(zstat = statistic, tstat = statistic)
}

try_kde_test <- function(x1, x2) {
  x1 <- as_finite_kde_input(x1)
  x2 <- as_finite_kde_input(x2)
  is_2d <- is.matrix(x1) || is.matrix(x2)

  if (length(x1) == 0 || length(x2) == 0) {
    return(ks_penalty_result("No finite values available for KS statistic."))
  }

  if (is_2d) {
    if (!is.matrix(x1) || !is.matrix(x2) || ncol(x1) != ncol(x2)) {
      return(ks_penalty_result("KS statistic inputs have incompatible dimensions."))
    }

    if (nrow(x1) < 2 || nrow(x2) < 2) {
      return(ks_penalty_result("Not enough finite rows available for KS statistic."))
    }

    fallback_statistic <- sliced_ks_statistic(x1, x2)
  } else {
    if (length(x1) < 2 || length(x2) < 2) {
      return(ks_penalty_result("Not enough finite values available for KS statistic."))
    }

    fallback_statistic <- empirical_ks_statistic(x1, x2)
  }

  result <- tryCatch(
    ks::kde.test(x1 = x1, x2 = x2),
    error = function(e) e
  )

  if (!inherits(result, "error")) {
    return(result)
  }

  reason <- paste0("ks::kde.test failed: ", result$message, ".")
  if (!is.finite(fallback_statistic)) {
    return(ks_penalty_result(reason))
  }

  fallback_ks_result(fallback_statistic, reason)
}

subsample_correlations <- function(x, max_size = 10000L, seed = 1L) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) <= max_size) {
    return(x)
  }

  set.seed(seed)
  sample(x, size = max_size, replace = FALSE)
}

cat("Computing ks statistic of fraction of zeros per gene\n")
frac_zero_real_genes <- colMeans(real_counts == 0)
frac_zero_sim_genes <- colMeans(sim_counts == 0)
ks_statistic_frac_zero_genes <- try_kde_test(
  x1 = frac_zero_real_genes,
  x2 = frac_zero_sim_genes
)

cat("Computing ks statistic of fraction of zeros per cell\n")
frac_zero_real_cells <- rowMeans(real_counts == 0)
frac_zero_sim_cells <- rowMeans(sim_counts == 0)
ks_statistic_frac_zero_cells <- try_kde_test(
  x1 = frac_zero_real_cells,
  x2 = frac_zero_sim_cells
)

cat("Computing ks statistic of the library size\n")
lib_size_real_cells <- log1p(rowSums(real_counts))
lib_size_sim_cells <- log1p(rowSums(sim_counts))
ks_statistic_lib_size_cells <- try_kde_test(
  x1 = lib_size_real_cells,
  x2 = lib_size_sim_cells
)

real_dge <- edgeR::calcNormFactors(
  edgeR::DGEList(counts = Matrix::t(real_counts)),
  method = "TMM"
)
sim_dge <- edgeR::calcNormFactors(
  edgeR::DGEList(counts = Matrix::t(sim_counts)),
  method = "TMM"
)

cat("Computing ks statistic of the effective library size\n")
efflib_size_real_cells <- log1p(
  real_dge$samples$lib.size * real_dge$samples$norm.factors
)
efflib_size_sim_cells <- log1p(
  sim_dge$samples$lib.size * sim_dge$samples$norm.factors
)
ks_statistic_efflib_size_cells <- try_kde_test(
  x1 = efflib_size_real_cells,
  x2 = efflib_size_sim_cells
)

cat("Computing ks statistic of TMM\n")
tmm_real_cells <- real_dge$samples$norm.factors
tmm_sim_cells <- sim_dge$samples$norm.factors
ks_statistic_tmm_cells <- try_kde_test(x1 = tmm_real_cells, x2 = tmm_sim_cells)

cat("Computing ks statistic of the cell-level scaled variance\n")
scaled_var_real_cells <- scale(sparseMatrixStats::colVars(Matrix::t(
  real_counts
)))
scaled_var_sim_cells <- scale(sparseMatrixStats::colVars(Matrix::t(sim_counts)))
ks_statistic_scaled_var_cells <- try_kde_test(
  x1 = as.numeric(scaled_var_real_cells),
  x2 = as.numeric(scaled_var_sim_cells)
)

cat("Computing ks statistic of the cell-level scaled mean\n")
scaled_mean_real_cells <- scale(colMeans(Matrix::t(real_counts)))
scaled_mean_sim_cells <- scale(colMeans(Matrix::t(sim_counts)))
ks_statistic_scaled_mean_cells <- try_kde_test(
  x1 = as.numeric(scaled_mean_real_cells),
  x2 = as.numeric(scaled_mean_sim_cells)
)

cat(
  "Computing ks statistic of the library size vs fraction of zeros per cell\n"
)
lib_fraczero_real_cells <- data.frame(
  lib = lib_size_real_cells,
  fraczero = frac_zero_real_cells
)
lib_fraczero_sim_cells <- data.frame(
  lib = lib_size_sim_cells,
  fraczero = frac_zero_sim_cells
)
ks_statistic_lib_fraczero_cells <- try_kde_test(
  x1 = lib_fraczero_real_cells,
  x2 = lib_fraczero_sim_cells
)

cat("Computing ks statistic of the sample Pearson correlation\n")
# pearson_real_cells <- reshape2::melt(cor(as.matrix(Matrix::t(real_counts)), method = "pearson"))
pearson_real_cells <- proxyC::simil(
  real_counts,
  margin = 1,
  method = "correlation"
)
# pearson_sim_cells <- reshape2::melt(cor(as.matrix(Matrix::t(sim_counts)), method = "pearson"))
pearson_sim_cells <- proxyC::simil(
  sim_counts,
  margin = 1,
  method = "correlation"
)

ks_statistic_pearson_cells <- try_kde_test(
  x1 = subsample_correlations(pearson_real_cells),
  x2 = subsample_correlations(pearson_sim_cells)
)

cat("Computing ks statistic of the gene-level scaled variance\n")
scaled_var_real_genes <- scale(sparseMatrixStats::colVars(real_counts))
scaled_var_sim_genes <- scale(sparseMatrixStats::colVars(sim_counts))
ks_statistic_scaled_var_genes <- try_kde_test(
  x1 = as.numeric(scaled_var_real_genes),
  x2 = as.numeric(scaled_var_sim_genes)
)

cat("Computing ks statistic of the gene-level scaled mean\n")
scaled_mean_real_genes <- scale(colMeans(real_counts))
scaled_mean_sim_genes <- scale(colMeans(sim_counts))
ks_statistic_scaled_mean_genes <- try_kde_test(
  x1 = as.numeric(scaled_mean_real_genes),
  x2 = as.numeric(scaled_mean_sim_genes)
)

cat("Computing ks statistic of the gene Pearson correlation\n")
# pearson_real_genes <- reshape2::melt(cor(as.matrix(real_counts), method = "pearson"))
pearson_real_genes <- proxyC::simil(
  real_counts,
  margin = 2,
  method = "correlation"
)
# pearson_sim_genes <- reshape2::melt(cor(as.matrix(sim_counts), method = "pearson"))
pearson_sim_genes <- proxyC::simil(
  sim_counts,
  margin = 2,
  method = "correlation"
)
ks_statistic_pearson_genes <- try_kde_test(
  x1 = subsample_correlations(pearson_real_genes),
  x2 = subsample_correlations(pearson_sim_genes)
)

cat("Computing ks statistic of the mean expression vs variance expression\n")
mean_var_real_genes <- data.frame(
  mean = colMeans(real_counts),
  var = sparseMatrixStats::colVars(real_counts)
)
mean_var_sim_genes <- data.frame(
  mean = colMeans(sim_counts),
  var = sparseMatrixStats::colVars(sim_counts)
)
ks_statistic_mean_var_genes <- try_kde_test(
  x1 = mean_var_real_genes,
  x2 = mean_var_sim_genes
)

cat(
  "Computing ks statistic of the mean expression vs fraction of zeros per gene\n"
)
mean_fraczero_real_genes <- data.frame(
  mean = colMeans(real_counts),
  fraczero = frac_zero_real_genes
)
mean_fraczero_sim_genes <- data.frame(
  mean = colMeans(sim_counts),
  fraczero = frac_zero_sim_genes
)
ks_statistic_mean_fraczero_genes <- try_kde_test(
  x1 = mean_fraczero_real_genes,
  x2 = mean_fraczero_sim_genes
)

cat("Combining metric values\n")
uns_metric_ids <- c(
  "ks_statistic_frac_zero_genes_zstat",
  "ks_statistic_frac_zero_cells_zstat",
  "ks_statistic_lib_size_cells_zstat",
  "ks_statistic_efflib_size_cells_zstat",
  "ks_statistic_tmm_cells_zstat",
  "ks_statistic_scaled_var_cells_zstat",
  "ks_statistic_scaled_mean_cells_zstat",
  "ks_statistic_lib_fraczero_cells_zstat",
  "ks_statistic_pearson_cells_zstat",
  "ks_statistic_scaled_var_genes_zstat",
  "ks_statistic_scaled_mean_genes_zstat",
  "ks_statistic_pearson_genes_zstat",
  "ks_statistic_mean_var_genes_zstat",
  "ks_statistic_mean_fraczero_genes_zstat",

  "ks_statistic_frac_zero_genes_tstat",
  "ks_statistic_frac_zero_cells_tstat",
  "ks_statistic_lib_size_cells_tstat",
  "ks_statistic_efflib_size_cells_tstat",
  "ks_statistic_tmm_cells_tstat",
  "ks_statistic_scaled_var_cells_tstat",
  "ks_statistic_scaled_mean_cells_tstat",
  "ks_statistic_lib_fraczero_cells_tstat",
  "ks_statistic_pearson_cells_tstat",
  "ks_statistic_scaled_var_genes_tstat",
  "ks_statistic_scaled_mean_genes_tstat",
  "ks_statistic_pearson_genes_tstat",
  "ks_statistic_mean_var_genes_tstat",
  "ks_statistic_mean_fraczero_genes_tstat"
)
uns_metric_values <- c(
  ks_statistic_frac_zero_genes$zstat,
  ks_statistic_frac_zero_cells$zstat,
  ks_statistic_lib_size_cells$zstat,
  ks_statistic_efflib_size_cells$zstat,
  ks_statistic_tmm_cells$zstat,
  ks_statistic_scaled_var_cells$zstat,
  ks_statistic_scaled_mean_cells$zstat,
  ks_statistic_lib_fraczero_cells$zstat,
  ks_statistic_pearson_cells$zstat,
  ks_statistic_scaled_var_genes$zstat,
  ks_statistic_scaled_mean_genes$zstat,
  ks_statistic_pearson_genes$zstat,
  ks_statistic_mean_var_genes$zstat,
  ks_statistic_mean_fraczero_genes$zstat,

  ks_statistic_frac_zero_genes$tstat,
  ks_statistic_frac_zero_cells$tstat,
  ks_statistic_lib_size_cells$tstat,
  ks_statistic_efflib_size_cells$tstat,
  ks_statistic_tmm_cells$tstat,
  ks_statistic_scaled_var_cells$tstat,
  ks_statistic_scaled_mean_cells$tstat,
  ks_statistic_lib_fraczero_cells$tstat,
  ks_statistic_pearson_cells$tstat,
  ks_statistic_scaled_var_genes$tstat,
  ks_statistic_scaled_mean_genes$tstat,
  ks_statistic_pearson_genes$tstat,
  ks_statistic_mean_var_genes$tstat,
  ks_statistic_mean_fraczero_genes$tstat
)

cat("Writing output AnnData to file\n")
output <- anndataR::AnnData(
  uns = list(
    dataset_id = input_simulated_dataset$uns[["dataset_id"]],
    method_id = input_simulated_dataset$uns[["method_id"]],
    metric_ids = uns_metric_ids,
    metric_values = uns_metric_values
  ),
  shape = c(0L, 0L)
)
output$write_h5ad(par[["output"]], compression = "gzip")
