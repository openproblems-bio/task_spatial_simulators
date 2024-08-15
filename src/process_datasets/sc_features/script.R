suppressMessages(library(SingleCellExperiment, quietly = TRUE))

## VIASH START
par <- list(
  # inputs
  input_sp = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  sc_species = "Mus musculus",
  # outputs
  output_sp = "resources_test/datasets/MOBNEW/dataset_sp.h5ad"
)
## VIASH END

cat("Read input files\n")
input_sp <- anndata::read_h5ad(par$input_sp)

cat("Spatial dataset:\n")
print(input_sp)


cat("Run scFeatures\n")
sce <- scater::logNormCounts(SingleCellExperiment(
  list(counts = Matrix::t(input_sp$layers[["counts"]])),
  colData = input_sp$obs,
  metadata = input_sp$obsm
))

feat_types <- c("L_stats","celltype_interaction","nn_correlation","morans_I")
log_count <- assay(sce, "logcounts")

prob_matrix <- sce@metadata$celltype_prop
colnames(prob_matrix) <- paste0("ct", seq_len(ncol(prob_matrix)))
rownames(prob_matrix) <- colnames(log_count)

scfeatures_result <- scFeatures::scFeatures(log_count,
                                sample = rep("sample1", ncol(log_count)),
                                spatialCoords = list(colData(sce)$row,colData(sce)$col),
                                feature_types = c("L_stats","celltype_interaction","nn_correlation","morans_I"),
                                type = "spatial_t",
                                species = sc_species,
                                spotProbability =  t(prob_matrix))
cat("Transforming spatial into AnnData\n")

# scFeatures

output_sp <- anndata::AnnData(
  layers = list(
    counts = input_sp$layers[["counts"]],
    logcounts = input_sp$layers[["logcounts"]]
  ),
  obs = input_sp$obs,
  var = input_sp$var,
  obsm = input_sp$obsm,
  uns = list(
    dataset_id = input_sp$uns["dataset_id"],
    dataset_name = input_sp$uns["dataset_name"],
    dataset_description = input_sp$uns["dataset_description"],
    dataset_url = input_sp$uns["dataset_url"],
    dataset_reference = input_sp$uns["dataset_reference"],
    dataset_summary = input_sp$uns["dataset_summary"],
    dataset_organism = input_sp$uns["dataset_organism"],
    L_stats = scfeatures_result$L_stats,
    celltype_interaction = scfeatures_result$celltype_interaction,
    nn_correlation = scfeatures_result$nn_correlation,
    morans_I = scfeatures_result$morans_I
  )
)

cat("Write output files\n")
output_sp$write_h5ad(par$output_sp, compression = "gzip")
