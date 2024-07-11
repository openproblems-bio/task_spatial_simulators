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
frac_zero_real_genes <- colMeans(real_counts == 0)
frac_zero_sim_genes <- colMeans(sim_counts == 0)
ks_statistic_frac_zero_genes <- ks::kde.test(x1 = frac_zero_real_genes, x2 = frac_zero_sim_genes)

frac_zero_real_cells <- rowMeans(real_counts == 0)
frac_zero_sim_cells <- rowMeans(sim_counts == 0)
ks_statistic_frac_zero_cells <- ks::kde.test(x1 = frac_zero_real_cells, x2 = frac_zero_sim_cells)

cat("Combine metric values\n")
uns_metric_ids <- c(
  "ks_statistic_frac_zero_genes",
  "ks_statistic_frac_zero_cells"
)
uns_metric_values <- c(
  ks_statistic_frac_zero_genes$zstat,
  ks_statistic_frac_zero_cells$zstat
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
