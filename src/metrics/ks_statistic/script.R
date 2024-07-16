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
ks_statistic_lib_size_cells <- ks::kde.test(x1 = lib_size_real_cells, x2 = lib_size_sim_cells)

# ks_statistic_efflib_size_cells
efflib_size_real_cells <- log1p(rowSums(real_counts))
efflib_size_sim_cells <- log1p(rowSums(sim_counts))
ks_statistic_efflib_size_cells <- ks::kde.test(x1 = efflib_size_real_cells, x2 = efflib_size_sim_cells)

# ks_statistic_tmm_cells (works only for variability in the data)
tmm_real_cells <- edgeR::DGEList(counts = Matrix::t(real_counts))$samples$norm.factors
tmm_sim_cells <- edgeR::DGEList(counts = Matrix::t(sim_counts))$samples$norm.factors
ks_statistic_tmm_cells <- ks::kde.test(x1 = tmm_real_cells, x2 = tmm_sim_cells)

# ks_statistic_scaled_var_cells
scaled_var_real_cells <- scale(resample::colVars(Matrix::t(real_counts)))
scaled_var_sim_cells <- scale(resample::colVars(Matrix::t(sim_counts)))
ks_statistic_scaled_var_cells <- ks::kde.test(x1 = as.numeric(scaled_var_sim_cells), x2 = as.numeric(scaled_var_sim_cells))

# ks_statistic_scaled_mean_cells
scaled_mean_real_cells <- scale(colMeans(Matrix::t(real_counts)))
scaled_mean_sim_cells <- scale(colMeans(Matrix::t(sim_counts)))
ks_statistic_scaled_mean_cells <- ks::kde.test(x1 = as.numeric(scaled_mean_sim_cells), x2 = as.numeric(scaled_mean_sim_cells))

# ks_statistic_lib_fraczero_cells
lib_fraczero_real_cells <- data.frame(lib = lib_size_real_cells, fraczero = frac_zero_real_cells)
lib_fraczero_sim_cells <- data.frame(lib = lib_size_sim_cells, fraczero = frac_zero_sim_cells)
ks_statistic_lib_fraczero_cells <- ks::kde.test(x1 = lib_fraczero_real_cells, x2 = lib_fraczero_sim_cells)

# ks_statistic_pearson_cells
pearson_real_cells <- reshape2::melt(cor(as.matrix(Matrix::t(real_counts)), method = "pearson"))
pearson_sim_cells <- reshape2::melt(cor(as.matrix(Matrix::t(sim_counts)), method = "pearson"))
ks_statistic_pearson_cells <- ks::kde.test(x1 = as.numeric(pearson_real_cells$value), x2 = as.numeric(pearson_sim_cells$value))

# ks_statistic_scaled_var_genes
scaled_var_real_genes <- scale(resample::colVars(real_counts))
scaled_var_sim_genes <- scale(resample::colVars(sim_counts))
ks_statistic_scaled_var_genes <- ks::kde.test(x1 = as.numeric(scaled_var_sim_genes), x2 = as.numeric(scaled_var_sim_genes))

# ks_statistic_scaled_mean_genes
scaled_mean_real_genes <- scale(colMeans(real_counts))
scaled_mean_sim_genes <- scale(colMeans(sim_counts))
ks_statistic_scaled_mean_genes <- ks::kde.test(x1 = as.numeric(scaled_mean_real_genes), x2 = as.numeric(scaled_mean_sim_genes))

# ks_statistic_pearson_genes
pearson_real_genes <- reshape2::melt(cor(as.matrix(real_counts), method = "pearson"))
pearson_sim_genes <- reshape2::melt(cor(as.matrix(sim_counts), method = "pearson"))
ks_statistic_pearson_genes <- ks::kde.test(x1 = as.numeric(pearson_real_genes$value), x2 = as.numeric(pearson_sim_genes$value))

# ks_statistic_mean_var_genes
mean_var_real_genes <- data.frame(mean = colMeans(real_counts), var = resample::colVars(real_counts))
mean_var_sim_genes <- data.frame(mean = colMeans(sim_counts), var = resample::colVars(sim_counts))
ks_statistic_mean_var_genes <- ks::kde.test(x1 = mean_var_real_genes, x2 = mean_var_sim_genes)

# ks_statistic_mean_fraczero_genes
mean_fraczero_real_genes <- data.frame(mean = colMeans(real_counts), fraczero = frac_zero_real_genes)
mean_fraczero_sim_genes <- data.frame(mean = colMeans(sim_counts), fraczero = frac_zero_sim_genes)
ks_statistic_mean_fraczero_genes <- ks::kde.test(x1 = mean_fraczero_real_genes, x2 = mean_fraczero_sim_genes)




cat("Combine metric values\n")
uns_metric_ids <- c(
  "ks_statistic_frac_zero_genes",
  "ks_statistic_frac_zero_cells",
  "ks_statistic_lib_size_cells",
  "ks_statistic_efflib_size_cells",
  "ks_statistic_scaled_var_cells",
  "ks_statistic_scaled_mean_cells",
  "ks_statistic_lib_fraczero_cells",
  "ks_statistic_pearson_cells",
  "ks_statistic_scaled_var_genes",
  "ks_statistic_scaled_mean_genes",
  "ks_statistic_pearson_genes",
  "ks_statistic_mean_var_genes",
  "ks_statistic_mean_fraczero_genes"
)
uns_metric_values <- c(
  ks_statistic_frac_zero_genes$zstat,
  ks_statistic_frac_zero_cells$zstat,
  ks_statistic_lib_size_cells$zstat,
  ks_statistic_efflib_size_cells$zstat,
  ks_statistic_scaled_var_cells$zstat,
  ks_statistic_scaled_mean_cells$zstat,
  ks_statistic_lib_fraczero_cells$zstat,
  ks_statistic_pearson_cells$zstat,
  ks_statistic_scaled_var_genes$zstat,
  ks_statistic_scaled_mean_genes$zstat,
  ks_statistic_pearson_genes$zstat,
  ks_statistic_mean_var_genes$zstat,
  ks_statistic_mean_fraczero_genes$zstat
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
