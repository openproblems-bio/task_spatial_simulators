requireNamespace("anndata", quietly = TRUE)
requireNamespace("edgeR", quietly = TRUE)
requireNamespace("ks", quietly = TRUE)
requireNamespace("resample", quietly = TRUE)
requireNamespace("reshape2", quietly = TRUE)
requireNamespace("scFeatures", quietly = TRUE)
requireNamespace("CARD", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)
requireNamespace("SummarizedExperiment", quietly = TRUE)


## VIASH START
par <- list(
  input_spatial_dataset = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  input_singlecell_dataset = "resources_test/datasets/MOBNEW/dataset_sc.h5ad",
  input_simulated_dataset = "resources_test/datasets/MOBNEW/simulated_dataset.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "ks_statistic",
  resources_dir = "target/executable/metrics/ks_statistics_scFeatures"
)
## VIASH END

source(paste0(meta[["resources_dir"]], "/utils.R"))

input_real_sp <- anndata::read_h5ad(par$input_spatial_dataset)
input_sc <- anndata::read_h5ad(par$input_singlecell_dataset)
input_simulated_sp <- anndata::read_h5ad(par$input_simulated_dataset)


cat("Reading input files\n")
real_sce <- scater::logNormCounts(SingleCellExperiment::SingleCellExperiment(
  list(counts = Matrix::t(input_real_sp$layers[["counts"]])),
  colData = input_real_sp$obs,
  metadata = input_real_sp$obsm
))

real_log_count <- SummarizedExperiment::assay(real_sce, "logcounts")
real_prob_matrix <- real_sce@metadata$celltype_prop
colnames(real_prob_matrix) <- paste0("ct", seq_len(ncol(real_prob_matrix)))
rownames(real_prob_matrix) <- colnames(real_log_count)

sim_sce <- scater::logNormCounts(SingleCellExperiment::SingleCellExperiment(
  list(counts = Matrix::t(input_simulated_sp$layers[["counts"]])),
  colData = input_simulated_sp$obs,
  metadata = input_simulated_sp$obsm
))

sim_log_count <- SummarizedExperiment::assay(sim_sce, "logcounts")

# build cell type deconvolution in simulated data
sim_prob_matrix <- CARD_processing(input_real_sp, input_sc)

feat_types <- c("L_stats","celltype_interaction","nn_correlation","morans_I")

real_scfeatures_result <- scFeatures::scFeatures(real_log_count,
                                sample = rep("sample1", ncol(real_log_count)),
                                spatialCoords = list(SingleCellExperiment::colData(real_sce)$row,SingleCellExperiment::colData(real_sce)$col),
                                feature_types = feat_types,
                                type = "spatial_t",
                                species = sc_species,
                                spotProbability =  t(real_prob_matrix))

sim_scfeatures_result <- scFeatures::scFeatures(sim_log_count,
                                sample = rep("sample1", ncol(sim_log_count)),
                                spatialCoords = list(SingleCellExperiment::colData(sim_sce)$row,SingleCellExperiment::colData(sim_sce)$col),
                                feature_types = feat_types,
                                type = "spatial_t",
                                species = sc_species,
                                spotProbability =  t(sim_prob_matrix))

ks_statistic_L_stats <- ks::kde.test(x1 = as.numeric(real_scfeatures_result$L_stats), x2 = as.numeric(sim_scfeatures_result$L_stats))
ks_statistic_celltype_interaction <- ks::kde.test(x1 = as.numeric(real_scfeatures_result$celltype_interaction), x2 = as.numeric(sim_scfeatures_result$celltype_interaction))
ks_statistic_nn_correlation <- ks::kde.test(x1 = as.numeric(real_scfeatures_result$nn_correlation), x2 = as.numeric(sim_scfeatures_result$nn_correlation))
ks_statistic_morans_I <- ks::kde.test(x1 = as.numeric(real_scfeatures_result$morans_I), x2 = as.numeric(sim_scfeatures_result$morans_I))

cat("Combining metric values\n")
uns_metric_ids <- c(
  "ks_statistic_L_stats",
  "ks_statistic_celltype_interaction",
  "ks_statistic_nn_correlation",
  "ks_statistic_morans_I"
)

uns_metric_values <- c(
  ks_statistic_L_stats,
  ks_statistic_celltype_interaction,
  ks_statistic_nn_correlation,
  ks_statistic_morans_I
)

cat("Writing output AnnData to file\n")
output <- anndata::AnnData(
  uns = list(
    dataset_id = input_simulated_sp$uns[["dataset_id"]],
    method_id = input_simulated_sp$uns[["method_id"]],
    metric_ids = uns_metric_ids,
    metric_values = uns_metric_values
  ),
  shape = c(0L, 0L)
)
output$write_h5ad(par[["output"]], compression = "gzip")