requireNamespace("anndata", quietly = TRUE)
requireNamespace("edgeR", quietly = TRUE)
requireNamespace("ks", quietly = TRUE)
requireNamespace("resample", quietly = TRUE)
requireNamespace("reshape2", quietly = TRUE)
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

add_kde_jitter <- function(x, amount) {
  x_range <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  x_scale <- max(x_range, stats::sd(as.numeric(x), na.rm = TRUE), 1)
  noise <- seq(-amount, amount, length.out = length(x)) * x_scale

  if (is.matrix(x)) {
    matrix(as.numeric(x) + noise, nrow = nrow(x), ncol = ncol(x))
  } else {
    as.numeric(x) + noise
  }
}

try_kde_test <- function(x1, x2) {
  x1 <- as_finite_kde_input(x1)
  x2 <- as_finite_kde_input(x2)

  if (length(x1) == 0 || length(x2) == 0) {
    warning("No finite values available for ks::kde.test; returning NA.")
    return(list(zstat = NA_real_, tstat = NA_real_))
  }

  if (is.matrix(x1) && (nrow(x1) < 2 || nrow(x2) < 2)) {
    warning("Not enough finite rows available for ks::kde.test; returning NA.")
    return(list(zstat = NA_real_, tstat = NA_real_))
  }

  if (!is.matrix(x1) && (length(x1) < 2 || length(x2) < 2)) {
    warning("Not enough finite values available for ks::kde.test; returning NA.")
    return(list(zstat = NA_real_, tstat = NA_real_))
  }

  last_error <- NULL
  for (jitter in c(0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4)) {
    x1_try <- if (jitter == 0) x1 else add_kde_jitter(x1, jitter)
    x2_try <- if (jitter == 0) x2 else add_kde_jitter(x2, jitter)

    result <- tryCatch(
      ks::kde.test(x1 = x1_try, x2 = x2_try),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )

    if (!is.null(result)) {
      if (jitter > 0) {
        warning(
          "Caught error in ks::kde.test: ",
          last_error$message,
          "\n\nSucceeded after adding deterministic jitter of size ",
          jitter,
          "."
        )
      }
      return(result)
    }
  }

  warning(
    "ks::kde.test failed after deterministic jitter retries: ",
    last_error$message,
    "\n\nReturning NA for this metric."
  )
  list(zstat = NA_real_, tstat = NA_real_)
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

cat("Computing ks statistic of the effective library size\n")
efflib_size_real_cells <- log1p(rowSums(real_counts))
efflib_size_sim_cells <- log1p(rowSums(sim_counts))
ks_statistic_efflib_size_cells <- try_kde_test(
  x1 = efflib_size_real_cells,
  x2 = efflib_size_sim_cells
)

cat("Computing ks statistic of TMM\n")
real_dge <- edgeR::DGEList(counts = Matrix::t(real_counts))
sim_dge <- edgeR::DGEList(counts = Matrix::t(sim_counts))
tmm_real_cells <- edgeR::calcNormFactors(
  real_dge,
  method = "TMM"
)$samples$norm.factors
tmm_sim_cells <- edgeR::calcNormFactors(
  sim_dge,
  method = "TMM"
)$samples$norm.factors
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
pearson_real_cells <- proxyC::simil(real_counts, method = "correlation")
# pearson_sim_cells <- reshape2::melt(cor(as.matrix(Matrix::t(sim_counts)), method = "pearson"))
pearson_sim_cells <- proxyC::simil(sim_counts, method = "correlation")

ks_statistic_pearson_cells <- try_kde_test(
  x1 = sample(as.numeric(pearson_real_cells), 10000),
  x2 = sample(as.numeric(pearson_sim_cells), 10000)
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
pearson_real_genes <- proxyC::simil(real_counts, method = "correlation")
# pearson_sim_genes <- reshape2::melt(cor(as.matrix(sim_counts), method = "pearson"))
pearson_sim_genes <- proxyC::simil(sim_counts, method = "correlation")
ks_statistic_pearson_genes <- try_kde_test(
  x1 = sample(as.numeric(pearson_real_genes), 10000),
  x2 = sample(as.numeric(pearson_sim_genes), 10000)
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
