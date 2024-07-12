requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_spatial_dataset = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  input_singlecell_dataset = "resources_test/datasets/MOBNEW/dataset_sc.h5ad",
  input_simulated_dataset = "resources_test/datasets/MOBNEW/simulated_dataset.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "ks_statistic"
)
## VIASH END

cat("Reading input files\n")
input_spatial_dataset <- anndata::read_h5ad(par[["input_spatial_dataset"]])
input_singlecell_dataset <- anndata::read_h5ad(par[["input_singlecell_dataset"]])
input_simulated_dataset <- anndata::read_h5ad(par[["input_simulated_dataset"]])

real_counts <- input_spatial_dataset$layers[["counts"]]
sim_counts <- input_simulated_dataset$layers[["counts"]]

cat("Compute ks statistic\n")
# ks_statistic_frac_zero_genes
frac_zero_real_genes <- colMeans(real_counts == 0)
frac_zero_sim_genes <- colMeans(sim_counts == 0)
ks_statistic_frac_zero_genes <- ks::kde.test(x1 = frac_zero_real_genes, x2 = frac_zero_sim_genes)

# ks_statistic_frac_zero_cells
frac_zero_real_cells <- rowMeans(real_counts == 0)
frac_zero_sim_cells <- rowMeans(sim_counts == 0)
ks_statistic_frac_zero_cells <- ks::kde.test(x1 = frac_zero_real_cells, x2 = frac_zero_sim_cells)

# ks_statistics_lib_size_cells
lib_size_real_cells <- log1p(rowSums(real_counts))
lib_size_sim_cells <- log1p(rowSums(sim_counts))
ks_statistics_lib_size_cells <- ks::kde.test(x1 = lib_size_real_cells, x2 = lib_size_sim_cells)

# ks_statistic_efflib_size_cells
# ks_statistic_tmm_cells
# ks_statistic_scaled_var_cells
# ks_statistic_scaled_mean_cells
# ks_statistic_lib_fraczero_cells
# ks_statistic_pearson_cells

# ks_statistic_scaled_var_genes
# ks_statistic_scaled_mean_genes
# ks_statistic_pearson_genes
# ks_statistic_mean_var_genes
# ks_statistic_mean_fraczeron_genes





cat("Combine metric values\n")
uns_metric_ids <- c(
  "ks_statistic_frac_zero_genes",
  "ks_statistic_frac_zero_cells",
  "ks_statistics_lib_size_cells"
)
uns_metric_values <- c(
  ks_statistic_frac_zero_genes$zstat,
  ks_statistic_frac_zero_cells$zstat,
  ks_statistics_lib_size_cells$zstat
)

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  uns = list(
    dataset_id = input_simulated_dataset$uns[["dataset_id"]],
    method_id = input_simulated_dataset$uns[["method_id"]],
    metric_ids = uns_metric_ids,
    metric_values = uns_metric_values
  ),
  shape = c(0L, 0L)
)
output$write_h5ad(par[["output"]], compression = "gzip")
